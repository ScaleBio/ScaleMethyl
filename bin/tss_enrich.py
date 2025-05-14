#!/usr/bin/env python
"""
Create csv file with number of reads mapping to tss and background bedfiles along with percentage enriched
"""
import pybedtools
from collections import defaultdict
from pathlib import Path
import argparse
import csv
import sys


def bam_intersect(id: str, bam_file: Path, tss_bed: Path, bg_bed: Path, cells: Path):
    passing_cells = []
    with open(cells) as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            if row["pass"] == "pass":
                passing_cells.append(row["cell_id"])

    bg_counts = defaultdict(int)
    background = pybedtools.BedTool(bg_bed)
    if len(background.intersect(background)) == 0:
        print("No background elements found. Background bed file may be empty or invalid.")
        sys.exit(1)
    bam = pybedtools.BedTool(bam_file)
    bg_overlaps = bam.intersect(background, u=True, bed=True, stream=True)

    for overlap in bg_overlaps:
        # intersect tool adds a slash plus 1 or 2 to qname to indicate which mate overlapped
        bc = overlap.name.split(":")[-1].split("/")[0]
        if bc in passing_cells:
            bg_counts[bc] += 1
    tss_counts = {key: 0 for key in bg_counts}
    tss = pybedtools.BedTool(tss_bed)
    if len(tss.intersect(tss)) == 0:
        print("No TSS elements found. TSS bed file may be empty or invalid.")
        sys.exit(1)
    tss_overlaps = bam.intersect(tss, u=True, bed=True, stream=True)

    for overlap in tss_overlaps:
        # intersect tool adds a slash plus 1 or 2 to qname to indicate which mate overlapped
        bc = overlap.name.split(":")[-1].split("/")[0]
        if bc in tss_counts:
            tss_counts[bc] += 1
    # write output csv file
    cols = ["cell_id", "tss_enrich", "tss_counts", "background_counts"]
    with open(f"{id}.tss_enrich.csv", mode="w") as output_file:
        csv_writer = csv.writer(output_file, delimiter=",")
        csv_writer.writerow(cols)
        for k, count in bg_counts.items():
            tss_count = tss_counts[k]
            enrich = tss_count / count
            row = [k, round(enrich, 4), tss_count, count]
            csv_writer.writerow(row)


def main():
    parser = argparse.ArgumentParser("Count overlaps of BAM to bed files")
    parser.add_argument("--id", required=True, help="For output file name")
    parser.add_argument("--bam", type=Path, required=True, help="Dedup BAM file")
    parser.add_argument("--tss_bed", type=Path, required=True, help="tss BED file")
    parser.add_argument("--bg_bed", type=Path, required=True, help="Background BED file")
    parser.add_argument("--cells", type=Path, required=True, help="Cells file")
    args = parser.parse_args()

    bam_intersect(
        id=args.id,
        bam_file=args.bam,
        tss_bed=args.tss_bed,
        bg_bed=args.bg_bed,
        cells=args.cells,
    )


if __name__ == "__main__":
    main()
