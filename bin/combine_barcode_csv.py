#!/usr/bin/env python
"""
Combine mulitple CSV files.
"""

import argparse
import pandas as pd
from pathlib import Path


def combine_csv(file_list: list,header:bool) -> pd.DataFrame:
    """
    Combine multiple CSV files into a single DataFrame.
    Args:
        file_list: list of CSV files to combine
    Returns:
        DataFrame containing the combined csv files
    """
    combined_df = []
    files = []
    print(header)
    for file in file_list:
        if header:
            read = pd.read_csv(file)
        else:
            read = pd.read_csv(file,header=None,names=["bc","data1","data2","data3","data4"],dtype={"bc":str,"data1":float,"data2":float,"data3":float,"data4":float})
        files.append(read)
    combined_df = pd.concat(files)
    return combined_df


def main():
    parser = argparse.ArgumentParser(description="Combine CSV files.")
    parser.add_argument(
        "--inFiles", nargs="+", type=Path, help="List of CSV files to combine."
    )
    parser.add_argument(
        "--outFile", type=Path, help="Output path for combined CSV."
    )
    parser.add_argument(
        "--header", action=argparse.BooleanOptionalAction, help="Is there a header in the CSV files?"
    )
    
    args = parser.parse_args()
    combined = combine_csv(args.inFiles, args.header)
    if args.header:
        combined.to_csv(args.outFile, index=False)
    else:
        combined.to_csv(args.outFile, index=False, header=False)


if __name__ == "__main__":
    main()
