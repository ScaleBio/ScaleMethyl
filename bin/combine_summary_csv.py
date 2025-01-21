#!/usr/bin/env python
"""
Combine mulitple dedup_stats.tsv and fragment_hist.tsv files. Some rows are summed and some are given a weighted average (% rows).
"""
import argparse


# Calculate the weighted average of the values in the two % rows, based on the number of reads in each sample
def calculate_weighted_average(values: dict, key: str, row: str) -> float:
    """
    Calculate the weighted average of the values in the two % rows, based on the number of reads in each sample
    Args:
        values: dictionary of values from file
        key: the key to calculate the weighted average for
        row: the key of the row to use for the weighting (usually "Passing Reads" or "Unique Reads")
    Returns:
        The weighted average of the values in a % row
    """
    weightedSum = 0.0
    total = 0.0
    for index, value in enumerate(values[key]):
        weightedSum += float(value) * float(values[row][index])
        total += float(values[row][index])
    if total != 0:
        weightedSum = weightedSum / total
    else:
        weightedSum = 0.0
    return weightedSum


# Read the values from the input files and store them in a dictionary and return the header
def input_values(files: list) -> tuple:
    """
    Read the values from the input files and store them in a dictionary and return the header
    Args:
        files: list of files to read

    Returns:
        A tuple containing the dictionary of values and the header
    """
    values = {}
    header = ""
    for file in files:
        header = file.readline().strip()
        for line in file:
            line = line.strip()
            if line == "":
                continue
            parts = line.split(",")
            if parts[0] not in values:
                values[parts[0]] = []
            values[parts[0]].append(float(parts[1]))
        file.close()
    return values, header


def main():
    parser = argparse.ArgumentParser("Combine summary csv files")
    parser.add_argument(
        "--inDedupFiles",
        type=argparse.FileType("r"),
        nargs="+",
        help="dedup_stats.tsv files to combine",
    )
    parser.add_argument(
        "--outDedupFile",
        type=argparse.FileType("w"),
        help="dedup_stats.tsv file to write",
    )
    parser.add_argument(
        "--inFragFiles",
        type=argparse.FileType("r"),
        nargs="+",
        help="fragment_hist.tsv files to combine",
    )
    parser.add_argument(
        "--outFragFile",
        type=argparse.FileType("w"),
        help="fragment_hist.tsv file to write",
    )
    args = parser.parse_args()

    # Perform calculations on the deduplication files
    values, header = input_values(args.inDedupFiles)

    args.outDedupFile.write(header + "\n")
    for key in values.keys():
        # The % rows are calculated differently, they are a weighted average based on the number of reads in each sample (a straight average would be misleading)
        if key == "Passing Reads %":
            weightedAverage = calculate_weighted_average(values, key, "Passing Reads")
            # write the weighted average to one decimal place to the output file
            args.outDedupFile.write(f"{key},{weightedAverage:.1f}\n")
        elif key == "Unique Reads %":
            weightedAverage = calculate_weighted_average(values, key, "Unique Reads")
            # write the weighted average to one decimal place to the output file
            args.outDedupFile.write(f"{key},{weightedAverage:.1f}\n")
        else:
            # write the sum to one decimal place to the output file
            args.outDedupFile.write(f"{key},{sum(values[key]):.1f}\n")

    args.outDedupFile.close()

    # Perform calculations on the fragment files
    values, header = input_values(args.inFragFiles)

    args.outFragFile.write(header + "\n")
    for key in values.keys():
        args.outFragFile.write(f"{key},{sum(values[key])}\n")
    args.outFragFile.close()


if __name__ == "__main__":
    main()
