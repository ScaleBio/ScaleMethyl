#!/usr/bin/env python
"""
Extract methylation calls from sorted dedup BAM file to context-specific parquet files
"""
import pysam
import pyarrow as pa
import pyarrow.parquet as pq
import numpy as np
from collections import defaultdict
import polars as pl
import re
from pathlib import Path
import duckdb
import argparse
from multiprocessing import Pool
from itertools import repeat
from pyfaidx import Fasta

BATCH_SIZE = 10_000


def write_batch(
    df: pl.DataFrame,
    chr: str,
    names: str,
    writer: pq.ParquetWriter,
    cur_pos: int = float("inf"),
):
    """
    Write methylation calls to chromosome-specific Parquet files for a batch of reads.
    Since we are grouping coverage by position, only write calls for positions less than the start position of current read.
    """
    batch = (
        df.filter(pl.col("pos") < cur_pos)
        .unique(
            # since paired-end reads have the same QNAME in the SAM file the unique function
            # will only retain one call for any position that is covered by both reads
            subset=["qname", "pos"]
        )
        .select(
            pl.col("qname").str.extract(r":([ACGT]+\+[ACGT]+\+[ACGT]+)$").alias("barcode"),
            pl.col("pos").cast(pl.UInt32),
            pl.col("strand"),
            pl.col("xb"),
        )
        .with_columns(
            pl.when(pl.col("xb").str.contains(r"[XYZ]")).then(1).otherwise(0).alias("methylated").cast(pl.UInt32),
            pl.when(pl.col("xb").str.contains(r"[xyz]")).then(1).otherwise(0).alias("unmethylated").cast(pl.UInt32),
            pl.when(pl.col("xb").str.contains(r"[xX]")).then(pl.lit("CG")).otherwise(pl.lit("CH")).alias("context"),
        )
        .group_by(["barcode", "pos"])
        .agg(
            pl.col("strand").first(),
            pl.col("context").first(),
            pl.lit(chr).alias("chr"),
            pl.col("methylated").sum(),
            pl.col("unmethylated").sum(),
        )
        .select(names)
        .sort("pos")
    )
    writer.write_table(batch.to_arrow())


def write_chr_parquet_bsbolt(bam: Path, chr: str, sample: str, threshold: float):
    """
    Write methylation calls to a Parquet file for a single chromosome
    """
    with pysam.AlignmentFile(bam, "rb") as samfile:
        high_ch_reads = defaultdict(int)
        schema = [
            ("qname", pl.String),
            ("pos", pl.UInt32),
            ("strand", pl.String),
            ("xb", pl.String),
        ]
        df = pl.DataFrame({}, schema=schema)
        calls = {"qname": [], "pos": [], "strand": [], "xb": []}
        pa_schema = pa.schema(
            [
                ("barcode", pa.large_string()),  # Polars uses large_string
                ("chr", pa.large_string()),
                ("pos", pa.uint32()),
                ("strand", pa.large_string()),
                ("context", pa.large_string()),
                ("methylated", pa.uint32()),
                ("unmethylated", pa.uint32()),
            ]
        )
        with pq.ParquetWriter(f"{sample}.RNAME_{chr}.parquet", pa_schema) as writer:
            reads = samfile.fetch(region=chr)
            num_reads = 0
            for read in reads:
                num_reads += 1
                if not read.has_tag("XB"):
                    continue  # skip reads without XB tag
                calls_idx = []
                ch_count = 0
                xb = read.get_tag("XB")
                if xb.isnumeric():
                    continue  # skip reads without methylation calls
                ch_count = xb.count("Y") + xb.count("Z")
                total_ch_count = xb.count("y") + xb.count("z") + ch_count
                if total_ch_count > 0 and (ch_count / total_ch_count) > threshold:
                    high_ch_reads[read.qname.split(":")[-1]] += 1
                    continue  # skip reads with CH methylation greater than threshold
                xb = re.findall(r"[xyzXYZ]|\d+", xb)
                aligned_pairs = dict(read.get_aligned_pairs())
                softClippedNumber = read.query_alignment_start
                # pos_offset is the offset to the aligned xb tag positions in the read, this is used to map the XB tags to the correct index in alignedpairs to be set in calls
                # Softclipped reads' xb tags don't align with the start of the read properly, you need to add the number of soft-clipped bases to the offset to get the correct position.
                pos_offset = (
                    np.array([1 if not char.isnumeric() else int(char) for char in xb]).cumsum() + softClippedNumber
                )
                calls_idx = [i for i in range(len(xb)) if not xb[i].isnumeric()]
                for idx in calls_idx:
                    calls["qname"].append(read.qname)
                    calls["pos"].append(aligned_pairs[pos_offset[idx] - 1])
                    calls["strand"].append(
                        {"C": "-", "W": "+"}[read.get_tag("YS")[0] if read.has_tag("YS") else "W"]
                    )  # YS tag identifies mapping strand
                    calls["xb"].append(xb[idx])

                if (num_reads % BATCH_SIZE) == 0:
                    df = df.extend(pl.DataFrame(calls, schema=schema))
                    write_batch(df, chr, pa_schema.names, writer, read.pos)
                    df = df.filter(pl.col("pos") >= read.pos)  # process in next batch
                    calls = {"qname": [], "pos": [], "strand": [], "xb": []}

            df = df.extend(pl.DataFrame(calls, schema=schema))
            write_batch(df, chr, pa_schema.names, writer)
        return high_ch_reads


def get_methylation_context(bases: str, meth_status: bool, reverse: bool) -> str:
    """
    Get the methylation context for a given base pair.
    Args:
        bases: Bases from the reference genome
        meth_status: Methylation status (True or False)
        reverse: Whether the read is on the reverse strand
    Returns:
        Methylation context (X, Y, Z, x, y, z)
    """
    checkBase = "G" if not reverse else "C"
    if bases[1] == checkBase:
        return "X" if meth_status else "x"
    elif bases[2] == checkBase:
        return "Y" if meth_status else "y"
    else:
        return "Z" if meth_status else "z"


def extract_methylation(
    ref_start: int,
    qname: str,
    refSeq: str,
    seq: str,
    x: int,
    y: int,
    ch_count: int,
    total_ch_count: int,
    reverse: bool,
    calls: dict,
) -> list:
    """
    Extract methylation calls from the read and reference sequence.

    Args:
        qlength: Length of the read
        ref_start: Start position of the read on the reference genome
        qname: Query name of the read
        refSeq: Reference sequence for the read
        seq: Read sequence
        x: Index for the read
        y: Index for the reference genome
        cigLen: Length of the cigar match
        ch_count: Count of CH methylation calls
        threshold: Threshold for CH methylation
        reverse: Whether the read is on the reverse strand
        calls: Dictionary to store methylation calls
        ref_base: Reference base (C or G)
        methyl_base: Methylated base (A or T)
    Returns:
        ch_count: Updated count of CH methylation calls
    """

    ref_base = "G" if reverse else "C"
    methyl_base = "A" if reverse else "T"
    # if the base on the reference genome is a C (G on reverse strand), test for methylation and contexts
    # y+2 should be equivalient to x for a cigar match, so y+i+2 is reference base we are testing
    if refSeq[y + 2] == ref_base:
        bases = ""
        # Take the base and the next two bases from the reference genome to test for context
        if not reverse:
            # Make sure we don't go out of bounds
            if y + 4 >= len(refSeq):
                return [ch_count,total_ch_count]

            # y+i+2 to y+i+4 is 3 bases to test for context on the forward strand
            bases = refSeq[y + 2 : y + 5 : 1]

        else:
            if y + 2 >= len(refSeq):
                return [ch_count,total_ch_count]
            # if y - 1 == -1, then the slicing breaks, so you have to index to the front
            if y - 1 < 0:
                bases = refSeq[y + 2 :: -1]
            else:
                bases = refSeq[y + 2 : y - 1 : -1]

        meth_state = ""
        # If the bp on the read is a C (G on reverse strand), it is methylated
        if seq[x] == ref_base:
            meth_state = get_methylation_context(bases, True, reverse)
        # If the bp on the read is an T (A on reverse strand), it is unmethylated
        elif seq[x] == methyl_base:
            meth_state = get_methylation_context(bases, False, reverse)
        if meth_state != "":
            if meth_state in ["Y", "Z"]:
                ch_count += 1
            if meth_state in ["y", "z", "Y", "Z"]:
                total_ch_count += 1
            calls["qname"].append(qname)
            calls["pos"].append(ref_start + y)
            calls["strand"].append("-" if reverse else "+")
            calls["xb"].append(meth_state)
    return [ch_count,total_ch_count]


def write_chr_parquet_bwa_meth(bam: Path, chr: str, sample: str, threshold: float, ref: Path) -> list:
    """
    Write methylation calls to a Parquet file for a single chromosome
    """
    fasta = Fasta(ref)
    with pysam.AlignmentFile(bam, "rb") as samfile:
        high_ch_reads = defaultdict(int)
        schema = [
            ("qname", pl.String),
            ("pos", pl.UInt32),
            ("strand", pl.String),
            ("xb", pl.String),
        ]
        df = pl.DataFrame({}, schema=schema)
        calls = {"qname": [], "pos": [], "strand": [], "xb": []}
        pa_schema = pa.schema(
            [
                ("barcode", pa.large_string()),  # Polars uses large_string
                ("chr", pa.large_string()),
                ("pos", pa.uint32()),
                ("strand", pa.large_string()),
                ("context", pa.large_string()),
                ("methylated", pa.uint32()),
                ("unmethylated", pa.uint32()),
            ]
        )
        with pq.ParquetWriter(f"{sample}.RNAME_{chr}.parquet", pa_schema) as writer:
            reads = samfile.fetch(region=chr)
            num_reads = 0
            for read in reads:
                temp_calls = {"qname": [], "pos": [], "strand": [], "xb": []}
                num_reads += 1
                ch_count = 0
                total_ch_count = 0
                if not read.cigarstring:
                    continue
                seq = read.query_sequence.upper()
                # reference sequence for the read with 2 extra bp at the start and 4 extra bp at the end to context testing
                start = read.reference_start - 2
                if start < 0:
                    continue
                end = read.reference_start + read.reference_length + 4

                # get the reference sequence for the read
                refSeq = str(fasta[chr][start:end].seq).upper()
                reverse = read.is_reverse
                # x is the index for the read, y is the index for the reference genome; y+2 should follow x for cigar match
                pairs = read.get_aligned_pairs(matches_only=True)
                high_ch = False
                for pair in pairs:
                    x = pair[0]
                    y = pair[1] - read.reference_start

                    if y < 0:
                        continue
                    # extract methylation calls from the read and reference sequence
                    [ch_count,total_ch_count] = extract_methylation(
                        read.reference_start,
                        read.qname,
                        refSeq,
                        seq,
                        x,
                        y,
                        ch_count,
                        total_ch_count,
                        reverse,
                        temp_calls,
                    )
                    # stop processing read if CH methylation greater than threshold
                    if total_ch_count > 0 and (ch_count / total_ch_count) > threshold:
                        high_ch_reads[read.qname.split(":")[-1]] += 1
                        high_ch = True
                        break
                if high_ch:
                    continue
                else:
                    calls["qname"].extend(temp_calls["qname"])
                    calls["pos"].extend(temp_calls["pos"])
                    calls["strand"].extend(temp_calls["strand"])
                    calls["xb"].extend(temp_calls["xb"])

                if (num_reads % BATCH_SIZE) == 0:
                    df = df.extend(pl.DataFrame(calls, schema=schema))
                    write_batch(df, chr, pa_schema.names, writer, read.pos)
                    df = df.filter(pl.col("pos") >= read.pos)  # process in next batch
                    calls = {"qname": [], "pos": [], "strand": [], "xb": []}

            df = df.extend(pl.DataFrame(calls, schema=schema))
            write_batch(df, chr, pa_schema.names, writer)
        fasta.close()
        return high_ch_reads


def generate_cell_info_db_prompt(sample, contexts):
    """
    Generate the DuckDB prompt to combine the context-specific cellInfo Parquet files
    into a single cellInfo.txt file.

    COPY (
        SELECT
            CellID,
            COALESCE(CG_Cov, 0) + COALESCE(CH_Cov, 0) AS Coverage,
            COALESCE(CG_Cov, 0) AS CG_Cov,
            COALESCE(CG_mC_Pct, 0) AS CG_mC_Pct,
            COALESCE(CH_Cov, 0) AS CH_Cov,
            COALESCE(CH_mC_Pct, 0) AS CH_mC_Pct
        FROM '{f"{sample}.cellInfo_{context1}.parquet"}'
        FULL OUTER JOIN '{f"{sample}_cellInfo_{context2}.parquet"}' USING (CellID)
    ) TO '{f"{sample}.cellInfo.txt"}' (HEADER, DELIMITER '\t');

    Args:
        sample: Sample name
        contexts: List of contexts to combine
    Returns:
        cellInfoDBPrompt: DuckDB prompt to combine the context-specific cellInfo Parquet files
    """
    sample2 = sample.replace(".", "_")
    cellInfoDBPrompt = "COPY (\n    SELECT\n        CellID,"

    # Add Coverage calculation
    coverage_parts = [f"COALESCE({context}_Cov, 0)" for context in contexts]
    cellInfoDBPrompt += f"\n        {' + '.join(coverage_parts)} AS Coverage,"

    # Add individual context columns
    for context in contexts:
        cellInfoDBPrompt += f"\n        COALESCE({context}_Cov, 0) AS {context}_Cov,"
        cellInfoDBPrompt += f"\n        COALESCE({context}_mC_Pct, 0) AS {context}_mC_Pct,"

    # Remove the trailing comma from the last context column
    cellInfoDBPrompt = cellInfoDBPrompt.rstrip(",")

    # Add FROM and JOIN clauses
    from_clause = f"\n    FROM '{sample2}_cellInfo_{contexts[0]}.parquet'"
    join_clauses = [
        f"\n    FULL OUTER JOIN '{sample2}_cellInfo_{context}.parquet' USING (CellID)" for context in contexts[1:]
    ]

    cellInfoDBPrompt += from_clause + "".join(join_clauses)
    cellInfoDBPrompt += f"\n) TO '{sample}.cellInfo.txt' (HEADER, DELIMITER '\t');"

    return cellInfoDBPrompt


def met_extract(bam: Path, sample: str, threshold: float, nproc: int, contexts: list, aligner: str, ref: Path):
    """
    Process dedup BAM file aligned by bsbolt to extract methylation calls
    to context-specific Parquet files.

    Take XB tag from bsbolt alignment output and extract the methylation
    calls creating a record in the Parquet file for each barcode, chromosome,
    position, and count of reads with methylated and unmethylated calls.

    Args:
        bam: Path to BAM file
        sample: Sample name
        threshold: Discard reads with CH methylation greater than threshold
        nproc: Number of subprocesses to use for parallel processing
        contexts: List of contexts to extract
        aligner: Aligner used to generate the BAM file
        ref: Path to reference genome
    """
    pysam.index(bam.as_posix())  # str argument required by pysam
    with pysam.AlignmentFile(bam, "rb") as samfile:
        # get all contigs with reads
        chrs = [stat.contig for stat in samfile.get_index_statistics() if stat.total > 0]

    with Pool(nproc) as p:
        if aligner == "bsbolt":
            chr_high_ch_reads = p.starmap(
                write_chr_parquet_bsbolt, zip(repeat(bam), chrs, repeat(sample), repeat(threshold))
            )
        elif aligner == "bwa-meth" or aligner == "parabricks":
            chr_high_ch_reads = p.starmap(
                write_chr_parquet_bwa_meth, zip(repeat(bam), chrs, repeat(sample), repeat(threshold), repeat(ref))
            )

    # create CG and CH context met calls parquet file
    # by combining all chromosome-specific parquet files
    # NB: Switched from polars.scan_parquet to duckdb to avoid issue with some records being re-ordered
    chr_pqs = Path(".").glob(f"{sample}.RNAME_*.parquet")
    # sort by chromosome name
    chr_pqs = sorted([file.as_posix() for file in chr_pqs])
    cell_infos = {}
    for context in contexts:
        duckdb.sql(
            f"""
        COPY (
            SELECT
                barcode,
                chr,
                pos,
                strand,
                methylated,
                unmethylated
            FROM
                read_parquet({chr_pqs})
            WHERE
                context = '{context}'
                
        )
        TO '{f"{sample}.met_{context}.parquet"}' (COMPRESSION ZSTD);
        """
        )

        # create cell info summary stats
        sample2 = sample.replace(".", "_")

        cell_infos[context] = duckdb.sql(
            f"""
        COPY(SELECT
            barcode AS CellID,
            SUM(methylated) + SUM(unmethylated) AS {context}_Cov,
            ROUND(SUM(methylated) / (SUM(methylated) + SUM(unmethylated)) * 100, 2) AS {context}_mC_Pct,
        FROM
            '{f"{sample}.met_{context}.parquet"}'
        GROUP BY
            barcode) 
        TO '{f"{sample2}_cellInfo_{context}.parquet"}';
        """
        )

    cellInfoDBPrompt = generate_cell_info_db_prompt(sample, contexts)
    duckdb.sql(cellInfoDBPrompt)

    # Add high CH reads stats to cellInfo.txt
    cell_info_df = pl.read_csv(f"{sample}.cellInfo.txt", separator="\t")
    high_ch_reads = {}
    for high_ch in chr_high_ch_reads:
        high_ch_reads.update(high_ch)
    if high_ch_reads:
        high_ch_df = pl.DataFrame({"CellID": high_ch_reads.keys(), "CH_high": high_ch_reads.values()})
        cell_info_df = cell_info_df.join(high_ch_df, on="CellID", how="left")
        cell_info_df = cell_info_df.with_columns(pl.col("CH_high").fill_null(0).alias("CH_high"))
    else:
        # No high CH reads
        cell_info_df = cell_info_df.with_columns(pl.lit(0).alias("CH_high"))
    cell_info_df.write_csv(f"{sample}.cellInfo.txt", separator="\t", include_header=True)


def main():
    parser = argparse.ArgumentParser("Create methylation calls parquet files from dedup bam file")
    parser.add_argument("bam", type=Path, help="Coordinate-sorted dedup bam file")
    parser.add_argument("--sample", type=str, help="Sample name")
    parser.add_argument(
        "--threshold",
        type=float,
        help="Discard reads with CH methylation greater than threshold",
    )
    parser.add_argument(
        "--subprocesses",
        type=int,
        help="How many subprocesses to use for parallel processing",
        default=4,
    )
    parser.add_argument(
        "--contexts",
        type=str,
        help="Comma-separated list of contexts to extract",
        default="CG,CH",
    )
    parser.add_argument("--aligner", type=str, help="Aligner used to generate the BAM file", default="bwa-meth")
    parser.add_argument("--ref", type=Path, help="Path to reference genome")
    args = parser.parse_args()
    met_extract(
        args.bam, args.sample, args.threshold, args.subprocesses, args.contexts.split(","), args.aligner, args.ref
    )


if __name__ == "__main__":
    main()
