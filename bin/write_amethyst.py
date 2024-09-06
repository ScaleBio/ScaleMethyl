#!/usr/bin/env python
"""
Write Amethyst hdf5 coverage files, one file per TN5 well
https://github.com/lrylaarsdam/amethyst
"""
import duckdb
import argparse
import h5py
import numpy as np
from pathlib import Path
from reporting.reporting import Utils


def write_amethyst(met_cg:Path, met_ch:Path, barcodes:list, sample:str):
    """
    Writes Amethyst hdf5 files

    Args:
        met_cg: Parquet file containing CG context methylation calls.
        met_ch: Parquet file containing CH context methylation calls.
        barcodes: Passing barcodes for this well coordinate.
        sample: Sample name with well coordinate for output file name.
    """
    with h5py.File(f"{sample}_cov.h5", 'w') as f:
        cg_group = f.create_group("CG")
        ch_group = f.create_group("CH")
        for bc in barcodes:
            for group, file in zip([cg_group, ch_group], [met_cg, met_ch]):
                cov = duckdb.query(f"""
                    SELECT
                        chr,
                        pos,
                        ROUND(methylated / (methylated + unmethylated), 2) as pct,
                        unmethylated as t,
                        methylated as c
                    FROM read_parquet('{file}')
                    WHERE barcode = '{bc}';
                """).fetchnumpy()
                # https://github.com/adeylab/premethyst/blob/main/premethyst_commands/calls2h5.py
                arr = np.zeros(len(cov['chr']), dtype = [('chr', 'S10'), ('pos', int), ('pct', float), ('t', int), ('c', int)])
                for field in cov.keys():
                    arr[field] = cov[field]
                group.create_dataset(bc, data=arr, compression='gzip', compression_opts=9)


def main():
    parser = argparse.ArgumentParser("Write Amethyst h5 file")
    parser.add_argument("--met_cg", type=Path, help="Parquet file with CG met calls")
    parser.add_argument("--met_ch", type=Path, help="Parquet file with CH met calls")
    parser.add_argument("--all_cells", type=Path, help="File with all cell barcode information")
    parser.add_argument("--sample", required=True, help="Sample name with well coordinate")
    args = parser.parse_args()
    write_amethyst(
        met_cg=args.met_cg,
        met_ch=args.met_ch,
        barcodes=Utils.get_passing_cells(args.all_cells, args.sample),
        sample = args.sample
        )
    
if __name__ == "__main__":
    main()
