#!/usr/bin/env python
"""
Generate bed file with genomic windows
"""
import argparse
from pathlib import Path


def parse_args():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(description="Create BED file with genomic windows")
    parser.add_argument("--genomeTiles", type=int, required=True, help="Size of the windows")
    parser.add_argument("--chromSizes", type=Path, required=True, help="Path to the chrom.sizes file")
    parser.add_argument("--blacklist", type=Path, required=False, help="Path to the blacklist file")
    parser.add_argument("--output", type=Path, required=True, help="Path to the output BED file")
    parser.add_argument(
        "--minimumWindowSize",
        type=int,
        required=False,
        default=1,
        help="Minimum size of the windows",
    )
    return parser.parse_args()


def create_windows(
    chrom_sizes_file: Path,
    blacklist: Path,
    window_size: int,
    output_file: Path,
    minimum_window_size: int = 1,
):
    """
    Create a BED file with genomic windows from a window size and a chrom.sizes file

    Args:
        chrom_sizes_file: path to the chrom.sizes file
        blacklist: path to the blacklist file
        window_size: size of the windows
        output_file: path to the output BED file
    """
    # Read blacklist regions into a set of tuples
    blacklist_regions = {}
    if blacklist:
        with open(blacklist, "r") as bl:
            for line in bl:
                chrom, start, end = line.strip().split()
                if chrom not in blacklist_regions:
                    blacklist_regions[chrom] = set()
                blacklist_regions[chrom].add((int(start), int(end)))

    # Create windows using the chrom.sizes file
    with open(chrom_sizes_file, "r") as f, open(output_file, "w") as out:
        # for each chromosome, create windows of size window_size starting at 0
        for line in f:
            chrom, size = line.strip().split()
            size = int(size)
            for start in range(0, size, window_size):
                end = min(start + window_size, size)
                # Check if the window overlaps with any blacklist region
                if chrom in blacklist_regions:
                    for bl_start, bl_end in blacklist_regions[chrom]:
                        if start <= bl_start <= end:
                            # blacklist region is in the middle of the window, so split the window
                            if start <= bl_end <= end:
                                if bl_end - start >= minimum_window_size:
                                    out.write(f"{chrom}\t{start}\t{bl_end}\n")
                                start = bl_end
                            # blacklist region overlaps the end of the window, so end the window sooner
                            else:
                                end = bl_start
                        # blacklist region overlaps the start of the window, so start the window later
                        elif start <= bl_end <= end:
                            start = bl_end
                        elif bl_start <= start and bl_end >= end:
                            continue
                # Write the window to the output file
                if end - start >= minimum_window_size:
                    out.write(f"{chrom}\t{start}\t{end}\n")


def main():
    args = parse_args()
    create_windows(
        args.chromSizes,
        args.blacklist,
        args.genomeTiles,
        args.output,
        args.minimumWindowSize,
    )


if __name__ == "__main__":
    main()
