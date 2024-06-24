#!/usr/bin/env python
"""
Merge matrices from individual processes
"""
import argparse
import gzip
from pathlib import Path
from scipy.io import mmread, mmwrite
from scipy.sparse import hstack


def merge_mtx(sample:str, barcode_files:list[Path], features_file:Path, mtx_files:list[Path]):
    """
    Merge score matrices
    
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
    
    mtx = [mmread(file) for file in mtx_files]
    score_str = ".score" if context == "CG" else ""
    with gzip.open(f"{sample}.{context}{score_str}.mtx.gz", 'wb') as f:
        mmwrite(f, hstack(mtx))
    
    with open(f"{sample}.barcodes.tsv", 'w') as f:
        for bc in barcodes:
            f.write(f"{bc}\n")
            
    # write out features file for this context
    features_file.rename(f"{sample}.{context}.features.tsv")


def main():
    parser = argparse.ArgumentParser("Merge met and cov matrices into dense mtx format")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--barcodes", nargs="+", type=Path, required=True, help="Barcode tsv files")
    parser.add_argument("--features", type=Path, required=True, help="Features tsv file")
    parser.add_argument("--scoreMtx", nargs="+", type=Path, required=True, help="score.mtx files to merge")
    args = parser.parse_args()

    merge_mtx(
        sample=args.sample,
        barcode_files=args.barcodes,
        features_file=args.features,
        mtx_files=args.scoreMtx
    )
    

if __name__ == "__main__":
    main()
