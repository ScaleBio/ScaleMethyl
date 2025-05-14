#!/usr/bin/env python
"""
Pre-process samples.csv to simplify parsing downstream.
Rename deprecriated column names and fill in defaults
Basic validation

Writes the normalized samples.csv to stdout
"""
import argparse
import csv
import sys
import json
import pandas as pd
from pathlib import Path
from typing import Dict


def validateName(name):
    """
    Check sample and library names for invalid characters
    Print error and exit for invalid names

    Args:
        name: sample/library name to check
    """
    for n in name:
        if not (n.isalnum() or n in "-"):
            print(
                f"Name should only contain [a-z],[A-Z],[0-9], dash (-): {name}",
                file=sys.stderr,
            )
            sys.exit(1)


def get_barcodes_range(libraryStruct: Path) -> str:
    """
    Get upper and lower end of permissible range for barcodes
    Args:
        libraryStruct: Path to library structure json file
    """
    # Parse lib json
    with open(libraryStruct) as f:
        libJson = json.load(f)

    if "sample_barcode" not in libJson:
        raise ValueError("sample_barcode not found in library structure json")

    # Get file corresponding to sample barcode so maximum well coordinate string can be computed
    for entry in libJson["barcodes"]:
        if entry["name"] == libJson["sample_barcode"]:
            start, end = get_first_and_last_entry(libraryStruct.parent / f"{entry['sequences']}")

    return f"{start}-{end}"


def get_first_and_last_entry(file: Path, sep: str = "\t") -> tuple[str, str]:
    """
    Get the first and last entry of a file
    Args:
        file: Path to file
        sep: Separator used in the file
    Returns:
        (first_entry, last_entry): Tuple of first and last entry
    """
    with open(file) as f:
        lines = f.readlines()
        first_entry = lines[0].strip().split(sep)[1]
        last_entry = lines[-1].strip().split(sep)[1]

    return first_entry, last_entry


def check_whether_barcode_ranges_overlap(all_samples_barcode_range: Dict, libraryStruct: Path):
    """
    Check whether the user supplied barcode ranges overlap amongst samples
    Throw exception if there is an overlap

    Args:
        all_samples_barcode_range: Dictionary of barcode ranges for all samples with libName as key
        libraryStruct: Path to the library structure json file
    """
    lib_struct_dir = libraryStruct.parent
    lib_struct = json.load(open(libraryStruct))
    # Get the txt file that corresponds to the sample_barcode in the lib json
    for barcode_info in lib_struct["barcodes"]:
        if lib_struct["sample_barcode"] == barcode_info["name"]:
            fname = barcode_info["sequences"]
            break
    barcode_whitelist = pd.read_csv(lib_struct_dir / fname, sep="\t", names=["barcode", "well"])
    for libName in all_samples_barcode_range:
        # Needs to be reset for each libName
        verbose_barcodes_list = []
        for sample_barcodes in all_samples_barcode_range[libName]:
            # String can be a single barcode or a range of barcodes separated by semi-colon
            for semi_colon_separated in sample_barcodes.split(";"):
                if "-" in semi_colon_separated:
                    # Get well coordinate that corresponds to starting of the barcodes range for a sample
                    starting = barcode_whitelist.index[barcode_whitelist["well"] == semi_colon_separated.split("-")[0]][
                        0
                    ]
                    # Get well coordinate that corresponds to end of the barcodes range for a sample
                    end = barcode_whitelist.index[barcode_whitelist["well"] == semi_colon_separated.split("-")[1]][0]
                    # Retrieve all the well coordinates that correspond to the barcodes range for a sample
                    all_barcodes = barcode_whitelist.loc[starting:end, "well"].tolist()
                    # extend because all_barcodes is a list
                    verbose_barcodes_list.extend(all_barcodes)
                else:
                    verbose_barcodes_list.append(semi_colon_separated)
        # Check whether the barcode ranges overlap amongst individual samples
        if len(verbose_barcodes_list) != len(set(verbose_barcodes_list)):
            print(
                "The barcodes range mentioned for each sample overlap amongst individual samples",
                file=sys.stderr,
            )
            sys.exit(1)


def main(samplesCsv: Path, splitFastq: bool, libraryStruct: Path, merge: bool):
    """
    Writes normalized samples.csv to stdout
    Args:
        samplesCsv: Path to samples.csv file
        splitFastq: Flag that denotes whether bcParser will split per TN5 well to ensure parallelization
        libraryStruct: Path to the library structure json file
        merge: samples.csv is for merging
    """
    barcodes_range = get_barcodes_range(libraryStruct)
    rows = {}  # sample -> samples.csv row
    with open(samplesCsv) as csvfile:
        samples = csv.reader(csvfile)
        cols = next(samples)
        # Trim empty columns
        while not cols[-1]:
            cols = cols[:-1]

        if "sample" not in cols:
            print("'sample' (the sample name) is required", file=sys.stderr)
            sys.exit(1)
        for row in samples:
            if not row or len(row[0].strip()) == 0:
                continue
            # Trim empty columns
            while not row[-1] and len(row) > len(cols):
                row = row[:-1]

            if len(row) != len(cols):
                print("Unexpected number of columns:", row, sep="\n", file=sys.stderr)
                sys.exit(1)
            sample = row[cols.index("sample")].strip()
            libName = row[cols.index("libName")].strip()
            id = sample + "-" + libName
            validateName(id)
            if merge:
                id = id + "-" + row[cols.index("resultDir")].strip()
            if id in rows:
                if merge:
                    print(f"Duplicate sample/lib name/resultDir combination: {id}", file=sys.stderr)
                else:
                    print(f"Duplicate sample/lib name combination: {id}", file=sys.stderr)
                sys.exit(1)
            rows[id] = row
    if "barcodes" not in cols:
        cols.insert(len(cols), "barcodes")
        for r in rows.values():
            r.insert(len(cols), "1A01-3H12")
    if "libName" not in cols and "fastqName" in cols:
        cols[cols.index("fastqName")] = "libName"
    if "libIndex" not in cols and "index" in cols:
        cols[cols.index("index")] = "libIndex"
    if "libIndex" not in cols and "fastqIndex" in cols:
        cols[cols.index("fastqIndex")] = "libIndex"
    # Normalize index sequences
    all_samples_barcode_range = {}
    if "libIndex" in cols:
        libIndexInd = cols.index("libIndex")
        for s, r in rows.items():
            index = "".join(r[libIndexInd].split())  # Remove whitespace
            index = index.upper()
            index = ";".join(sorted(index.split(";")))  # sort sequences
            r[libIndexInd] = index
        # Default libName based on index
        if "libName" not in cols:
            cols.insert(libIndexInd, "libName")
            for r in rows.values():
                name = r[libIndexInd].split(";")[0]
                r.insert(libIndexInd, name)
    elif "libName" not in cols:  # Only one (unnamed) library; Use default name
        cols.insert(1, "libName")
        for r in rows.values():
            r.insert(1, "ScaleMethyl")

    libNameIndex = cols.index("libName")
    for r in rows:
        all_samples_barcode_range[r[libNameIndex]] = []
    if "barcodes" not in cols:
        cols.insert(len(cols), "barcodes")
        for r in rows.values():
            r.insert(len(cols), barcodes_range)
            all_samples_barcode_range[r[libNameIndex]].append(barcodes_range)
    else:  # check whether each entry for barcodes is empty or not
        barcodesIndex = cols.index("barcodes")
        for r in rows.values():
            if r[barcodesIndex].strip() == "":
                r[barcodesIndex] = barcodes_range
                all_samples_barcode_range[r[libNameIndex]].append(barcodes_range)
    check_whether_barcode_ranges_overlap(all_samples_barcode_range, libraryStruct)
    if splitFastq:
        if "split" not in cols:
            cols.insert(len(cols), "split")
            for r in rows.values():
                r.insert(len(cols), "true")

    w = csv.writer(sys.stdout)
    w.writerow(cols)
    w.writerows(rows.values())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Standardize samples.csv for the workflow (defaults, column names, etc."
    )
    parser.add_argument(
        "samples",
        metavar="SAMPLES.csv",
        type=Path,
        help="CSV with sample names and information for the workflow run",
    )
    parser.add_argument(
        "--splitFastq",
        help="Flag that denotes whether bcParser will split per TN5 well to ensure parallelization",
        action="store_true",
        default=False,
    )
    parser.add_argument("--libraryStruct", help="Library structure json file", type=Path)
    parser.add_argument("--merge", help="Merge samples.csv files", action="store_true", default=False)
    args = parser.parse_args()
    main(args.samples, args.splitFastq, args.libraryStruct, args.merge)
