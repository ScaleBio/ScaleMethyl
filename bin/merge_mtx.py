#!/usr/bin/env python
"""
Merge matrices from individual processes
"""
import argparse
import gzip
from pathlib import Path
from scipy.io import mmread, mmwrite
from scipy.sparse import csc_array, hstack
from reporting.reporting import Utils
import numpy as np


def merge_mtx(
    sample: str,
    barcode_files: list[Path],
    features_file: Path,
    mtx_files: list[Path],
    passing_bc: list[str],
):
    """
    Merge matrices

    Args:
        sample: Sample name
        barcode_files: Barcode files
        features_file: Features tsv file
        mtx_files: mtx files to merge
        passing_bc: Passing barcodes
    """
    context = "CG" if ".CG" in mtx_files[0].suffixes else "CH"

    # Read features to get the number of features
    with open(features_file, "r") as f:
        feature_count = sum(1 for _ in f)

    # Initialize an empty sparse matrix with the correct shape
    mtx = csc_array((feature_count, 0), dtype=np.float32)  # Start with an empty matrix with float32 type
    column_names = []
    idx_dict = {val: idx for idx, val in enumerate(passing_bc)}
    for mtx_file, barcode_file in zip(mtx_files, barcode_files):
        with open(mtx_file, "rb") as f, open(barcode_file, "r") as bar:
            matrix = csc_array(mmread(f), dtype=np.float32)  # Ensure the matrix is of type float32
            bcs = [line.strip() for line in bar.readlines()]
            filtered_bc = [bc for bc in bcs if bc in passing_bc]
            column_names.extend(filtered_bc)
            indices = [bidx for bidx in range(len(bcs)) if bcs[bidx] in idx_dict]
            matrix = matrix[:, indices]
            mtx = hstack([mtx, matrix], format="csc")

    # get indices to re-sort matrix to match barcodes order in all_cells
    score_str = ".score" if context == "CG" else ""
    comment = "\nCG context score matrix\n" if context == "CG" else "\nCH context methylation rate matrix\n"
    with gzip.open(f"{sample}.{context}{score_str}.mtx.gz", "wb", compresslevel=6) as f:
        mmwrite(f, mtx, comment=comment, precision=3)

    with open(f"{sample}.barcodes.tsv", "w") as f:
        for bc in column_names:
            f.write(f"{bc}\n")

    # write out features file for this context
    features_file.rename(f"{sample}.{context}.features.tsv")


def main():
    parser = argparse.ArgumentParser("Merge matrices into sparse mtx format")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--barcodes", nargs="+", type=Path, required=True, help="Barcode tsv files")
    parser.add_argument("--features", type=Path, required=True, help="Features tsv file")
    parser.add_argument("--mtx", nargs="+", type=Path, required=True, help=".mtx files to merge")
    parser.add_argument("--all_cells", type=Path, nargs="+", help="File with all cell info")
    args = parser.parse_args()
    merge_mtx(
        sample=args.sample,
        barcode_files=args.barcodes,
        features_file=args.features,
        mtx_files=args.mtx,
        passing_bc=Utils.get_passing_cells(args.all_cells, args.sample),
    )


if __name__ == "__main__":
    main()
