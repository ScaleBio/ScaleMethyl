#!/usr/bin/env python
"""
Process tsv files created by split_bam and extract methylation calls to context-specific parquet files
"""
import polars as pl
import re
from pathlib import Path
import duckdb
import argparse

def met_extract(input_dir:Path):
    """
    Processes TSV files created by split_bam and extracts methylation calls
    to context-specific Parquet files.

    Takes XB tag from bsbolt alignment output and extracts the methylation
    calls creating a record in the Parquet file for each barcode, chromosome,
    position, and met context identified in each read.

    Args:
        input_dir: Directory containing the TSV files to be processed.
    """
    # pickup batch files created by split_bam
    batch = input_dir.glob("*.RNAME_*.tsv")
    for chr in batch:
        # read tsv file lazily for potential memory savings
        df = pl.scan_csv(chr,
                         has_header=False,
                         separator="\t",
                         new_columns=["qname", "chr", "pos", "strand", "xb"],
                         dtypes=[pl.String, pl.String, pl.UInt32, pl.String, pl.String]
                         )

        # extract met calls from XB tag
        df = df.with_columns(
            pl.col("xb").str.extract_all(r"[xyzXYZ]|\d+"),
            pl.col("xb").str.extract_all(r"[xyzXYZ]|\d+")
                .cast(pl.List(pl.UInt8))
                .alias("lengths")
        ).with_row_index()

        # create a row for each met call
        df = df.explode(["xb", "lengths"]).with_columns(
            pl.col("lengths").fill_null(1).cum_sum().over(pl.col("index")).alias("pos_offset"),
        ).filter(
            # discard rows with digits in xb
            pl.col("xb").cast(pl.UInt8, strict=False).is_null()
        ).select(
            pl.col("qname"),
            pl.col("chr"),
            pl.col("pos").add(pl.col("pos_offset")).alias("pos"),
            pl.col("strand"),
            pl.col("xb").alias("met"),
        # since paired-end reads have the same QNAME in the SAM file the unique function
        # below will only retain one call for any position that is covered by both reads
        ).unique(subset=["qname", "chr", "pos"])

        df = df.select(
            pl.col("qname").str.extract(r":([ACGT]+\+[ACGT]+\+[ACGT]+)$").alias("barcode"),
            pl.col("chr"),
            pl.col("pos"),
            pl.col("strand"),
            pl.col("met")
        )
        df.collect().write_parquet(input_dir / f"{chr.stem}.parquet")
        
    sample = re.sub(r"\.RNAME_.*.tsv$", "", chr.name)

    # create CG context met calls parquet file
    batch = input_dir.glob(f"{sample}.RNAME_*.parquet")
    pl.scan_parquet(
        [file for file in batch]
    ).filter(
        pl.col("met").str.contains(r"[xX]")
    ).sink_parquet(input_dir / f"{sample}.met_CG.parquet")

    # create CH context met calls parquet file
    batch = input_dir.glob(f"{sample}.RNAME_*.parquet")
    pl.scan_parquet(
        [file for file in batch]
    ).filter(
        pl.col("met").str.contains(r"[yzYZ]")
    ).sink_parquet(input_dir / f"{sample}.met_CH.parquet")

    # create cell info summary stats
    duckdb.sql(f"""
    COPY (
        SELECT
            barcode AS CellID,
            COUNT(met) AS Coverage,
            SUM(CASE WHEN met ILIKE 'X' THEN 1 ELSE 0 END) AS CG_Cov,
            ROUND(SUM(CASE WHEN met = 'X' THEN 1 ELSE 0 END) / SUM(CASE WHEN met ILIKE 'X' THEN 1 ELSE 0 END) * 100, 2) AS CG_mC_Pct,
            SUM(CASE WHEN met SIMILAR TO '[yYzZ]' THEN 1 ELSE 0 END) AS CH_Cov,
            ROUND(SUM(CASE WHEN met SIMILAR TO '[YZ]' THEN 1 ELSE 0 END) / SUM(CASE WHEN met SIMILAR TO '[yYzZ]' THEN 1 ELSE 0 END) * 100, 2) AS CH_mC_Pct
        FROM
            '{input_dir / f"{sample}.met_*.parquet"}'
        GROUP BY
            barcode
    ) TO '{input_dir / f"{sample}.cellInfo.txt"}' (HEADER, DELIMITER '\t');
    """)
    
    # Add high CH reads stats to cellInfo.txt
    cell_info_df = pl.read_csv(input_dir / f"{sample}.cellInfo.txt", separator='\t')
    high_ch_reads = input_dir / f"{sample}.high_ch_reads.tsv"
    if high_ch_reads.exists():
        high_ch_df = pl.read_csv(high_ch_reads, separator='\t', has_header=False, new_columns=['CellID', 'CH_high'])
        cell_info_df = cell_info_df.join(high_ch_df, on='CellID', how='left')
        cell_info_df = cell_info_df.with_columns(
            pl.col("CH_high").fill_null(0).alias("CH_high")
        )
    else:
        # No high CH reads
        cell_info_df = cell_info_df.with_columns(pl.lit(0).alias("CH_high"))
        
    cell_info_df.write_csv(input_dir / f"{sample}.cellInfo.txt", separator='\t', include_header=True)

    

def main():
    parser = argparse.ArgumentParser("Create methylation calls parquet files from split bam files")
    parser.add_argument("input_dir", nargs='?', type=Path, default=Path.cwd(), help="Directory for input files")
    args = parser.parse_args()
    met_extract(args.input_dir)

if __name__ == "__main__":
    main()