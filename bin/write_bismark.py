#!/usr/bin/env python
"""
Write Bismark coverage files, one file per passing cell
https://github.com/anders-biostat/MethSCAn/blob/master/docs/tutorial.md#what-you-will-need
"""
import duckdb
import argparse
import subprocess
from pathlib import Path
from reporting.reporting import Utils

def write_bismark(met_calls:Path, barcodes:list, context:str):
    """
    Writes Bismark coverage format files, one file per passing cell.

    Args:
        met_calls: Parquet file containing methylation calls.
        barcodes: Passing barcodes for this well coordinate.
        context: CG or CH
    """
    # get chromosomes present in methylation extraction calls
    chrs = sorted([row[0] for row in duckdb.sql(f"SELECT DISTINCT(chr) FROM read_parquet('{met_calls}');").fetchall()])
    outdir = Path.cwd() / context
    outdir.mkdir(exist_ok=True)
    for bc in barcodes:
        for chr in chrs:
            duckdb.sql(f"""
            COPY (
                SELECT
                    any_value(chr),
                    pos,
                    ROUND(SUM(CASE WHEN met SIMILAR TO '[XYZ]' THEN 1 ELSE 0 END) / COUNT(met), 2) AS met_perc,
                    SUM(CASE WHEN met SIMILAR TO '[XYZ]' THEN 1 ELSE 0 END) as met_count,
                    SUM(CASE WHEN met SIMILAR TO '[xyz]' THEN 1 ELSE 0 END) as unmet_count
                FROM
                    read_parquet('{met_calls}')
                WHERE
                    barcode = '{bc}' AND chr = '{chr}'
                GROUP BY
                    pos
                ORDER BY
                    pos
            ) TO '{outdir / f"{bc}.{chr}.{context}.cov"}' (HEADER false, DELIMITER '\t');
            """)
        # concatenate all chromosome files into one file
        chr_files = sorted([file for file in outdir.glob(f"{bc}.*.{context}.cov")])
        cat_process = subprocess.Popen(["cat"] + chr_files, stdout=subprocess.PIPE)
        with open(outdir / f"{bc}.{context}.cov.gz", "wb") as f:
            gzip_process = subprocess.Popen(["bgzip", "-c"], stdin=cat_process.stdout, stdout=f)
            gzip_process.wait()
        # delete the chromosome files
        for file in chr_files:
            file.unlink()


def main():
    parser = argparse.ArgumentParser("Write per-cell Bismark coverage files")
    parser.add_argument("--met_calls", type=Path, help="Parquet file with met calls")
    parser.add_argument("--all_cells", type=Path, help="File with all cell barcode information")
    parser.add_argument("--sample", required=True, help="Sample name with well coordinate")
    args = parser.parse_args()
    write_bismark(
        met_calls=args.met_calls,
        barcodes=Utils.get_passing_cells(args.all_cells, args.sample),
        context = "CG" if ".met_CG" in args.met_calls.suffixes else "CH"
        )
    
if __name__ == "__main__":
    main()
