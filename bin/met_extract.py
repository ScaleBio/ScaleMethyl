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
BATCH_SIZE = 10_000

def write_batch(df:pl.DataFrame, chr:str, names:str, writer:pq.ParquetWriter, cur_pos:int=float("inf")):
    """
    Write methylation calls to chromosome-specific Parquet files for a batch of reads.
    Since we are grouping coverage by position, only write calls for positions less than the start position of current read.
    """
    batch = df.filter(pl.col('pos') < cur_pos).unique(
        # since paired-end reads have the same QNAME in the SAM file the unique function
        # will only retain one call for any position that is covered by both reads
        subset=["qname", "pos"]
        ).select(
            pl.col("qname").str.extract(r":([ACGT]+\+[ACGT]+\+[ACGT]+)$").alias("barcode"),
            pl.col("pos").cast(pl.UInt32),
            pl.col("strand"),
            pl.col("xb")
        ).with_columns(
            pl.when(pl.col("xb").str.contains(r"[XYZ]")).then(1).otherwise(0).alias("methylated").cast(pl.UInt32),
            pl.when(pl.col("xb").str.contains(r"[xyz]")).then(1).otherwise(0).alias("unmethylated").cast(pl.UInt32),
            pl.when(pl.col("xb").str.contains(r"[xX]")).then(pl.lit("CG")).otherwise(pl.lit("CH")).alias("context")
        ).group_by(["barcode", "pos"]).agg(
            pl.col("strand").first(),
            pl.col("context").first(),
            pl.lit(chr).alias("chr"),
            pl.col("methylated").sum(),
            pl.col("unmethylated").sum(),
        ).select(names).sort("pos")
    writer.write_table(batch.to_arrow())


def write_chr_parquet(bam:Path, chr:str, sample:str, threshold:float):
    """
    Write methylation calls to a Parquet file for a single chromosome
    """
    with pysam.AlignmentFile(bam, "rb") as samfile:
        high_ch_reads = defaultdict(int)
        schema = [
            ("qname", pl.String),
            ("pos", pl.UInt32),
            ("strand", pl.String),
            ("xb", pl.String)
        ]
        df = pl.DataFrame({}, schema=schema)
        calls = {'qname': [], 'pos': [], 'strand': [], 'xb': []}
        pa_schema = pa.schema([
            ("barcode", pa.large_string()), # Polars uses large_string
            ("chr", pa.large_string()),
            ("pos", pa.uint32()),
            ("strand", pa.large_string()),
            ("context", pa.large_string()),
            ("methylated", pa.uint32()),
            ("unmethylated", pa.uint32())
        ])
        with pq.ParquetWriter(f"{sample}.RNAME_{chr}.parquet", pa_schema) as writer:
            reads = samfile.fetch(region=chr)
            num_reads = 0
            for read in reads:
                num_reads += 1
                if not read.has_tag('XB'):
                    continue # skip reads without XB tag
                xb = read.get_tag('XB')
                if xb.isnumeric():
                    continue # skip reads without methylation calls
                ch_count = xb.count('Y') + xb.count('Z')
                if (ch_count / read.query_length) > threshold:
                    high_ch_reads[read.qname.split(':')[-1]] += 1
                    continue # skip reads with CH methylation greater than threshold
                xb = re.findall(r"[xyzXYZ]|\d+", xb)
                pos_offset = np.array([1 if not char.isnumeric() else int(char) for char in xb]).cumsum()
                calls_idx = [i for i in range(len(xb)) if not xb[i].isnumeric()]
                aligned_pairs = dict(read.get_aligned_pairs())
                for idx in calls_idx:
                    calls['qname'].append(read.qname)
                    calls['pos'].append(aligned_pairs[pos_offset[idx]-1])
                    calls['strand'].append({'C': '-', 'W': '+'}[read.get_tag('YS')[0] if read.has_tag('YS') else 'W']) # YS tag identifies mapping strand
                    calls['xb'].append(xb[idx])
                if (num_reads % BATCH_SIZE) == 0:
                    df = df.extend(pl.DataFrame(calls, schema=schema))
                    write_batch(df, chr, pa_schema.names, writer, read.pos)
                    df = df.filter(pl.col('pos') >= read.pos) # process in next batch
                    calls = {'qname': [], 'pos': [], 'strand': [], 'xb': []}
            df = df.extend(pl.DataFrame(calls, schema=schema))
            write_batch(df, chr, pa_schema.names, writer)
        return high_ch_reads
    

def met_extract(bam:Path, sample:str, threshold:float, nproc:int):
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
    """
    pysam.index(bam.as_posix()) # str argument required by pysam
    with pysam.AlignmentFile(bam, "rb") as samfile:
        # get all contigs with reads
        chrs = [stat.contig for stat in samfile.get_index_statistics() if stat.total > 0]
    with Pool(nproc) as p:
        chr_high_ch_reads = p.starmap(write_chr_parquet, zip(repeat(bam), chrs, repeat(sample), repeat(threshold)))
    
    # create CG and CH context met calls parquet file
    # by combining all chromosome-specific parquet files
    # NB: Switched from polars.scan_parquet to duckdb to avoid issue with some records being re-ordered
    chr_pqs = Path(".").glob(f"{sample}.RNAME_*.parquet")
    # sort by chromosome name
    chr_pqs = sorted([file.as_posix() for file in chr_pqs])

    # CG context
    duckdb.sql(f"""
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
            context = 'CG'
               
    )
    TO '{f"{sample}.met_CG.parquet"}' (COMPRESSION ZSTD);
    """)

    # CH context
    duckdb.sql(f"""
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
            context = 'CH'
               
    )
    TO '{f"{sample}.met_CH.parquet"}' (COMPRESSION ZSTD);
    """)

    # create cell info summary stats
    cg_info = duckdb.sql(f"""
    SELECT
        barcode AS CellID,
        SUM(methylated) + SUM(unmethylated) AS CG_Cov,
        ROUND(SUM(methylated) / (SUM(methylated) + SUM(unmethylated)) * 100, 2) AS CG_mC_Pct,
    FROM
        '{f"{sample}.met_CG.parquet"}'
    GROUP BY
        barcode
    """)
    ch_info = duckdb.sql(f"""
    SELECT
        barcode AS CellID,
        SUM(methylated) + SUM(unmethylated) AS CH_Cov,
        ROUND(SUM(methylated) / (SUM(methylated) + SUM(unmethylated)) * 100, 2) AS CH_mC_Pct,
    FROM
        '{f"{sample}.met_CH.parquet"}'
    GROUP BY
        barcode
    """)
    duckdb.sql(f"""
    COPY (
        SELECT
            CellID,
            COALESCE(CG_Cov, 0) + COALESCE(CH_Cov, 0) AS Coverage,
            COALESCE(CG_Cov, 0) AS CG_Cov,
            COALESCE(CG_mC_Pct, 0) AS CG_mC_Pct,
            COALESCE(CH_Cov, 0) AS CH_Cov,
            COALESCE(CH_mC_Pct, 0) AS CH_mC_Pct
        FROM cg_info
        FULL OUTER JOIN ch_info
        USING (CellID)
    ) TO '{f"{sample}.cellInfo.txt"}' (HEADER, DELIMITER '\t');
    """)
    
    # Add high CH reads stats to cellInfo.txt
    cell_info_df = pl.read_csv(f"{sample}.cellInfo.txt", separator='\t')
    high_ch_reads = {}
    for high_ch in chr_high_ch_reads:
        high_ch_reads.update(high_ch)
    if high_ch_reads:
        high_ch_df = pl.DataFrame({'CellID': high_ch_reads.keys(), 'CH_high': high_ch_reads.values()})
        cell_info_df = cell_info_df.join(high_ch_df, on='CellID', how='left')
        cell_info_df = cell_info_df.with_columns(
            pl.col("CH_high").fill_null(0).alias("CH_high")
        )
    else:
        # No high CH reads
        cell_info_df = cell_info_df.with_columns(pl.lit(0).alias("CH_high"))
    cell_info_df.write_csv(f"{sample}.cellInfo.txt", separator='\t', include_header=True)
    

def main():
    parser = argparse.ArgumentParser("Create methylation calls parquet files from dedup bam file")
    parser.add_argument("bam", type=Path, help="Coordinate-sorted dedup bam file aligned by bsbolt")
    parser.add_argument("--sample", type=str, help="Sample name")
    parser.add_argument("--threshold", type=float, help="Discard reads with CH methylation greater than threshold")
    parser.add_argument("--subprocesses", type=int, help="How many subprocesses to use for parallel processing", default=4)
    args = parser.parse_args()
    met_extract(args.bam, args.sample, args.threshold, args.subprocesses)

if __name__ == "__main__":
    main()