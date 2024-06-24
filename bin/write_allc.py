#!/usr/bin/env python
"""
Write ALLC format files, one file per passing cell
https://lhqing.github.io/ALLCools/start/input_files.html#allc-file
"""
import duckdb
import argparse
import subprocess
from pathlib import Path
from reporting.reporting import Utils

def write_allc(met_calls:Path, barcodes:list):
    """
    Writes ALLC format files, one file per passing cell.

    Args:
        met_calls: Parquet file containing methylation calls.
        barcodes: Passing barcodes for this well coordinate.
    """
    # get chromosomes present in methylation extraction calls
    chrs = sorted([row[0] for row in duckdb.sql(f"SELECT DISTINCT(chr) FROM read_parquet('{met_calls}');").fetchall()])

    outdir = Path.cwd() / "CG" if ".met_CG" in met_calls.suffixes else Path.cwd() / "CH"
    outdir.mkdir(exist_ok=True)
    for bc in barcodes:
        # batch by chromosome to reduce memory usage
        for chr in chrs:
            duckdb.sql(f"""
            COPY (
                SELECT
                    any_value(chr),
                    pos,
                    any_value(strand),
                    CASE WHEN any_value(met) SIMILAR TO '[xX]' THEN 'CG' WHEN any_value(met) SIMILAR TO '[yY]' THEN 'CHG' ELSE 'CHH' END as context,
                    SUM(CASE WHEN met SIMILAR TO '[XYZ]' THEN 1 ELSE 0 END) AS met,
                    COUNT(met) as cov,
                    1 as col7
                FROM
                    read_parquet('{met_calls}')
                WHERE
                    barcode = '{bc}' AND chr = '{chr}'
                GROUP BY
                    pos
                ORDER BY
                    pos
            ) TO '{outdir / f"{bc}.{chr}.allc.tsv"}' (HEADER false, DELIMITER '\t');
            """)
        # concatenate all chromosome files into one file
        chr_files = sorted([file for file in outdir.glob(f"{bc}.*.allc.tsv")])
        cat_process = subprocess.Popen(["cat"] + chr_files, stdout=subprocess.PIPE)
        with open(outdir / f"{bc}.allc.tsv.gz", "wb") as f:
            gzip_process = subprocess.Popen(["bgzip", "-c"], stdin=cat_process.stdout, stdout=f)
            gzip_process.wait()
        # delete the chromosome files
        for file in chr_files:
            file.unlink()
    
    
def main():
    parser = argparse.ArgumentParser("Write per-cell ALLC files")
    parser.add_argument("--sample", required=True, help="Sample name with well coordinate")
    parser.add_argument("--met_calls", type=Path, help="Parquet file with met calls")
    parser.add_argument("--all_cells", type=Path, help="File with all cell barcode information")
    args = parser.parse_args()
    write_allc(
        met_calls=args.met_calls,
        barcodes=Utils.get_passing_cells(args.all_cells, args.sample)
    )

    
if __name__ == "__main__":
    main()
