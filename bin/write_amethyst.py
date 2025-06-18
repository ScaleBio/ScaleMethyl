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


def write_amethyst(met_cg: list[Path], met_ch: list[Path], barcodes: list, sample: str):
    """
    Writes Amethyst hdf5 files with this file structure: https://github.com/lrylaarsdam/amethyst/raw/main/images/h5structure.png?raw=true

    Args:
        met_cg: List of Parquet files containing CG context methylation calls.
        met_ch: List of Parquet files containing CH context methylation calls.
        barcodes: Passing barcodes for this well coordinate.
        sample: Sample name with well coordinate for output file name.
    """
    with h5py.File(f"{sample}_cov.h5", "w") as f:
        mets_list = []
        mets_list.append(met_cg)
        context_groups = []
        context_groups.append(f.create_group("CG"))
        if met_ch != []:
            mets_list.append(met_ch)
            context_groups.append(f.create_group("CH"))
        
        for bc in barcodes:
            for group, files in zip(context_groups, mets_list):
                file_names = [file.name for file in files]
                cov = duckdb.query(
                    f"""
                    SELECT
                        chr,
                        pos,
                        ROUND(methylated / (methylated + unmethylated), 2) as pct,
                        unmethylated as t,
                        methylated as c
                    FROM read_parquet({file_names})
                    WHERE barcode = '{bc}';
                """
                ).fetchnumpy()
                arr = np.zeros(
                    len(cov["chr"]),
                    dtype=[
                        ("chr", "S10"),
                        ("pos", int),
                        ("pct", float),
                        ("t", int),
                        ("c", int),
                    ],
                )
                for field in cov.keys():
                    arr[field] = cov[field]
                post = group.create_group(bc)
                # base-level values are given the "1" dataset name
                post.create_dataset("1", data=arr, compression="gzip", compression_opts=9)


def main():
    parser = argparse.ArgumentParser("Write Amethyst h5 file")
    parser.add_argument("--met_cg", type=Path, nargs="+", help="Parquet file with CG met calls")
    parser.add_argument(
        "--met_ch",
        type=Path,
        nargs="+",
        help="Parquet file with CH met calls",
        default=[],
    )
    parser.add_argument(
        "--all_cells",
        type=Path,
        nargs="+",
        help="File with all cell barcode information",
    )
    parser.add_argument("--sample", required=True, help="Sample name with well coordinate")
    args = parser.parse_args()
    write_amethyst(
        met_cg=args.met_cg,
        met_ch=args.met_ch,
        barcodes=Utils.get_passing_cells(args.all_cells, args.sample),
        sample=args.sample,
    )


if __name__ == "__main__":
    main()
