#!/usr/bin/env python
"""
Merge matrices from individual processes
"""
import argparse
import gzip
from pathlib import Path
from scipy.io import mmread, mmwrite
from scipy.sparse import hstack, csc_array
from reporting.reporting import Utils


def merge_mtx(sample:str, barcode_files:list[Path], features_file:Path, mtx_files:list[Path], passing_bc:list[str]):
    """
    Merge matrices
    
    Args:
        sample: Sample name
        barcodes: Barcode tsv files
        features: Features tsv file
        mtx_files: mtx files to merge
    """
    context = "CG" if ".CG" in mtx_files[0].suffixes else "CH"
    barcodes = []
    for file in barcode_files:
        with open(file, 'r') as f:
            barcodes.extend([line.strip() for line in f.readlines()])
    
    # get number of features
    feature_count = 0
    with open(features_file, 'r') as f:
        for line in f:
            feature_count += 1
    
    mtx = hstack([csc_array(mmread(file)) for file in mtx_files])
    idx_dict = {val: idx for idx, val in enumerate(barcodes)}
    # get indices to re-sort matrix to match barcodes order in all_cells
    sorted_indices = [idx_dict[bc] for bc in passing_bc]
    mtx = mtx[:, sorted_indices]

    score_str = ".score" if context == "CG" else ""
    comment = '\nCG context score matrix\n' if context == "CG" else '\nCH context methylation rate matrix\n'
    with gzip.open(f"{sample}.{context}{score_str}.mtx.gz", 'wb', compresslevel=6) as f:
        mmwrite(f, mtx, comment=comment, precision=3)
    
    with open(f"{sample}.barcodes.tsv", 'w') as f:
        for bc in passing_bc:
            f.write(f"{bc}\n")
            
    # write out features file for this context
    features_file.rename(f"{sample}.{context}.features.tsv")


def main():
    parser = argparse.ArgumentParser("Merge matrices into sparse mtx format")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--barcodes", nargs="+", type=Path, required=True, help="Barcode tsv files")
    parser.add_argument("--features", type=Path, required=True, help="Features tsv file")
    parser.add_argument("--mtx", nargs="+", type=Path, required=True, help=".mtx files to merge")
    parser.add_argument("--all_cells", type=Path, help="File with all cell info")
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
