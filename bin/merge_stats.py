#!/usr/bin/env python
import argparse
import glob
import sys
import os
import re
import pandas as pd
from pathlib import Path


def accept_args() -> argparse.Namespace:
    """
    Accept command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Concatenate per well stats")
    parser.add_argument("--cellInfoPath", default=".", type=Path,
                        help="Path to directory with files containing methylation calls for each split file")
    parser.add_argument("--cellStatsPath", default=".", type=Path,
                        help="Path to directory containing split cell_stats.tsv files")
    parser.add_argument("--fragmentHistPath", default=".", type=Path,
                        help="Path to directory containing split fragment_hist.tsv files")
    parser.add_argument("--sampleName", type=str, required=True,
                        help="Identifier for this sample")
    parser.add_argument("--mappingStatsPath", default=".", type=Path,
                        help="Path to directory with log files containing bsbolt output")
    parser.add_argument("--trimmingStatsPath", default=".", type=Path,
                        help="Path to directory with log files containing trim_galore output")
    parser.add_argument("--dedupStatsPath", default=".", type=Path,
                        help="Path to directory containing dedup_stats.tsv files")
    parser.add_argument("--outDir", default=".", type=Path,
                        help="Path to write concatenated files to")
    args = parser.parse_args()
    return args


def concat_cell_stats(sampleName: str, outDir: Path, cellStatsPath: Path):
    """
    Concatenate cell_stats.tsv files for a sample
    """
    tsv_files = glob.glob(f"{cellStatsPath}/{sampleName}*.cell_stats.tsv")
    columns = ["BC", "total", "passing", "uniq", "MitoReads"]
    if tsv_files == []:
        print("No cell_stats.tsv file found", file=sys.stderr)
    if "MouseUnique" in next(open(tsv_files[0])):
        columns.extend(['HumanUnique', 'MouseUnique'])
    with open(f"{outDir}/{sampleName}.cell_stats.csv", "w") as outfile:
        outfile.write(f"{','.join(columns)}\n")
        for file in tsv_files:
            with open(file, "r") as infile:
                lines = infile.readlines()
                for line in lines:
                    if "Barcode" not in line:
                        outfile.write(line.replace("\t", ","))

def concat_dedup_stats(sampleName: str, outDir: Path, dedupStatsPath: Path):
    """
    Concatenate dedup_stats.tsv files for a sample
    """
    final_df = pd.DataFrame()
    file_list = glob.glob(f"{dedupStatsPath}/{sampleName}*.dedup_stats.tsv")
    if file_list == []:
        print("No dedup_stats.tsv file found", file=sys.stderr)
    for idx, fname in enumerate(file_list):
        df = pd.read_csv(fname, sep="\t", header=None, names=["Metric", "Value", "Value2"])
        df = df.drop("Value2",axis=1)
        if idx == 0:
            final_df = df
        else:
            final_df["Value"] = df["Value"].add(final_df["Value"])

    input_reads = final_df[final_df["Metric"] == "Input Reads"]["Value"].to_list()[0]
    passing_reads = final_df[final_df["Metric"] == "Passing Reads"]["Value"].to_list()[0]
    unique_reads = final_df[final_df["Metric"] == "Unique Reads"]["Value"].to_list()[0]
    final_df.loc[len(final_df)] = ['Passing Reads %', round(passing_reads/input_reads*100, 1)]
    final_df.loc[len(final_df)] = ['Unique Reads %', round(unique_reads/passing_reads*100, 1)]
    final_df.to_csv(f"{outDir}/{sampleName}.dedup_stats.csv", index=False)

def concat_fragment_histogram(sampleName: str, outDir: Path, fragmentHistPath: Path):
    """
    Concatenate fragment_hist.tsv files for a sample
    """
    file_list = glob.glob(f"{fragmentHistPath}/{sampleName}*.fragment_hist.tsv")
    if file_list == []:
        print("No fragment_hist.tsv file found", file=sys.stderr)
    data = {}
    for file_path in file_list:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip().split('\t')
                key = line[0]
                value = int(line[1])
                data[key] = data.get(key, 0) + value
    sorted_data = sorted(data.items(), key=lambda x: int(x[0]))
    with open(f"{outDir}/{sampleName}.fragment_hist.csv", 'w') as output:
        for key, value in sorted_data:
            output.write(f"{key},{value}\n")


def concat_cellinfo(sampleName: str, outDir: Path, cellInfoPath: Path):
    """
    Concatenate cellInfo.txt files for a sample
    """
    file_list = glob.glob(f"{cellInfoPath}/{sampleName}*.cellInfo.txt")
    if file_list == []:
        print("No cellInfo.txt file found", file=sys.stderr)
    with open(f"{outDir}/{sampleName}.cellInfo.csv", "w") as outfile:
        for file_name in file_list:
            with open(file_name, "r") as infile:
                contents = infile.read()
                # Exclude lines containing "CellID"
                filtered_contents = [line for line in contents.split("\n") if "CellID" not in line]
                outfile.write(("\n".join(filtered_contents)).replace("\t", ","))


def build_mapping_stats(sampleName: str, outDir: Path, mappingStatsPath: Path):
    """
    Build mapping stats from the bsbolt log files
    """
    suffix = "*bsbolt.log"
    fullPath = os.path.join(mappingStatsPath, sampleName + suffix)
    filenamesList = glob.glob(fullPath)
    if filenamesList == []:
        print("No logs found for bsbolt(mapping)", file=sys.stderr)
        return
    total_reads = []
    unaligned_reads = []
    for file in filenamesList:
        with open(file) as fh:
            for line in fh:
                if line.startswith("TotalAlignments:"):
                    count = int(re.sub('TotalAlignments: |\n', "", line))
                    total_reads.append(count)
                if line.startswith("Unaligned:"):
                    count = int(re.sub('Unaligned: |\n', "", line))
                    unaligned_reads.append(count)
    mapping_report = pd.DataFrame.from_dict({
        "Metric": ["total_reads", "unaligned_reads", "mapped_reads", "mapped_percent"],
        "Value": [int(sum(total_reads)), int(sum(unaligned_reads)), int(sum(total_reads) - sum(unaligned_reads)),
                  round(((sum(total_reads) - sum(unaligned_reads)) / sum(total_reads) * 100), 2)]
    }, dtype="object")
    mapping_report.to_csv(f'{outDir}/{sampleName}.mapping_stats.csv', index=False)


def build_trimming_stats(sampleName: str, outDir: Path, trimmingStatsPath: Path):
    """
    Build trimming stats from trimgalore report files
    """
    suffix = "*.fq.gz_trimming_report_R2.txt"
    fullPath = os.path.join(trimmingStatsPath, sampleName + suffix)
    filenamesList = glob.glob(fullPath)
    if filenamesList == []:
        print("No logs found for trim galore(trimming)", file=sys.stderr)
        return
    total_pairs = []
    adapter_pairs = []
    removed_pairs = []
    for file in filenamesList:
        with open(file) as fh:
            for line in fh:
                if line.startswith("Total reads processed:"):
                    count = int(re.sub('Total reads processed:|\s+|\n|,', "", line))
                    total_pairs.append(count)
                if line.startswith("Reads with adapters:"):
                    count = int(re.sub('Reads with adapters:|\s+|\n|,|\(.*$', "", line))
                    adapter_pairs.append(count)
                if line.startswith("Number of sequence pairs removed because at least one read"):
                    count = int(re.sub('^(.*:)|\s+|\n|,|\(.*$', "", line))
                    removed_pairs.append(count)
    trimming_report = pd.DataFrame.from_dict({
        "Metric": ["total_pairs", "adapter_pairs", "removed_pairs", "passing_pairs", "percent_passing", "total_reads", "passing_reads"],
        "Value": [int(sum(total_pairs)), int(sum(adapter_pairs)), int(sum(removed_pairs)), int(sum(total_pairs) - sum(removed_pairs)),
                  round(((sum(total_pairs) - sum(removed_pairs)) / sum(total_pairs) * 100), 2), int(sum(total_pairs)*2), int((sum(total_pairs) - sum(removed_pairs))*2)]
    }, dtype="object")
    trimming_report.to_csv(f'{outDir}/{sampleName}.trimming_stats.csv', index=False)


def main():
    args = accept_args()
    args.outDir.mkdir()
    concat_dedup_stats(args.sampleName, args.outDir, args.dedupStatsPath)
    concat_cell_stats(args.sampleName, args.outDir, args.cellStatsPath)
    concat_fragment_histogram(args.sampleName, args.outDir, args.fragmentHistPath)
    concat_cellinfo(args.sampleName, args.outDir, args.cellInfoPath)
    build_mapping_stats(args.sampleName, args.outDir, args.mappingStatsPath)
    build_trimming_stats(args.sampleName, args.outDir, args.trimmingStatsPath)


if __name__ == "__main__":
    main()
