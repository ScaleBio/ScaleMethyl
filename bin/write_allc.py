#!/usr/bin/env python
"""
Write ALLC format files, one file per passing cell
https://lhqing.github.io/ALLCools/start/input_files.html#allc-file
"""
import duckdb
import argparse
import subprocess
import threading
import os
from pathlib import Path
from reporting.reporting import Utils

def write_to_pipe(met_calls:Path, bc:str, context:str):
    """
    Send duckdb query to named pipe
    """
    duckdb.query(f"""
        COPY(
            SELECT
                chr,
                pos,
                strand,
                '{context}' as context,
                methylated,
                (methylated + unmethylated) as cov,
                1 as col7
            FROM read_parquet('{met_calls}')
            WHERE barcode = '{bc}'
        ) TO './duckdb.pipe' (FORMAT CSV, HEADER false, DELIMITER '\t');
    """)
    

def write_allc(met_calls:Path, barcodes:list):
    """
    Writes ALLC format files, one file per passing cell.

    Args:
        met_calls: Parquet file containing methylation calls.
        barcodes: Passing barcodes for this well coordinate.
    """
    context = "CG" if ".met_CG" in met_calls.suffixes else "CH"
    outdir = Path.cwd() / context
    outdir.mkdir(exist_ok=True)
    if not Path("duckdb.pipe").exists():
        os.mkfifo("duckdb.pipe")
    for bc in barcodes:
        # send query to thread which will write to named pipe (FIFO) that is input to bgzip
        query_thread = threading.Thread(target=write_to_pipe, args=(met_calls, bc, context))
        query_thread.start()
        with open(outdir / f"{bc}.allc.tsv.gz", "wb") as f:
            gzip_process = subprocess.Popen(["bgzip", "-c", "./duckdb.pipe"], stdout=f)
            gzip_process.communicate()
    os.unlink("duckdb.pipe")
    
    
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
