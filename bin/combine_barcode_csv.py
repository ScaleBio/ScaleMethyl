#!/usr/bin/env python
"""
Combine mulitple CSV files.
"""

import argparse
import pandas as pd
from pathlib import Path
from shutil import copyfileobj


def combine_csv(file_list: list, header: bool, outFile: Path):
    """
    Combine multiple CSV files into a single output file.
    Args:
        file_list: list of CSV files to combine
    Returns:
        DataFrame containing the combined csv files
    """

    with open(outFile, "w") as out:
        first_file = True
        for file_location in file_list:
            first_line = True
            with open(file_location, "r") as file:
                for line in file:
                    if first_line and not first_file:
                        first_line = False
                        if header:
                            continue
                    out.write(line.strip() + "\n")
            first_file = False


def main():
    parser = argparse.ArgumentParser(description="Combine CSV files.")
    parser.add_argument("--inFiles", nargs="+", type=Path, help="List of CSV files to combine.")
    parser.add_argument("--outFile", type=Path, help="Output path for combined CSV.")
    parser.add_argument("--header", action=argparse.BooleanOptionalAction, help="Is there a header in the CSV files?")

    args = parser.parse_args()
    combine_csv(args.inFiles, args.header, args.outFile)


if __name__ == "__main__":
    main()
