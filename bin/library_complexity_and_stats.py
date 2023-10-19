#!/usr/bin/env python
"""
Generate sample report and associated metrics
"""
import argparse
import pandas as pd
import plotly.express as px
import os
import re
import glob
import sys
import csv
import numpy as np
import datapane as dp
import shutil
from pathlib import Path
from reporting.reporting import BuildReadsPage, BuildMethylPage, build_plate_plot


def accept_args() -> argparse.Namespace:
    """
    Accept command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Generate sample report")
    parser.add_argument("--tgmt_barcodes", type=Path,
                        help="Path to file containing tgmt barcodes")
    parser.add_argument("--i5_barcodes", type=Path,
                        help="Path to file containing i5 barcodes")
    parser.add_argument("--i7_barcodes", type=Path,
                        help="Path to file containing i7 barcodes")
    parser.add_argument("--cellInfoPath", default=".", type=Path,
                        help="Path to directory with files containing methylation calls for each split file")
    parser.add_argument("--cellStatsPath", default=".", type=Path,
                        help="Path to directory containing split cell_stats.tsv files")
    parser.add_argument("--fragmentHistPath", default=".", type=Path,
                        help="Path to directory containing split fragment_hist.tsv files")
    parser.add_argument("--sampleName", type=str, required=True,
                        help="Identifier for this sample")
    parser.add_argument("--topCellPercentage", default=99.0, type=float,
                        help="Percentage of cells over minReads to use as 'robust max'")
    parser.add_argument("--minCellRatio", default=200.0, type=float,
                        help="Ratio between transcript counts of top cells and the lower cell threshold")
    parser.add_argument("--minReads", default=1000.0, type=float,
                        help="Minimum counts to consider a barcode as a potential cell")
    parser.add_argument("--minUniqTotal", default=0.0, type=float,
                        help="Minimum uniqCounts/Total percentage cutoff")
    parser.add_argument("--maxUniq", type=float,
                        help="Maximum unique counts cutoff (optional; default = max of tsv)")
    parser.add_argument("--maxUniqTotal", default=100.0, type=float,
                        help="Maximum uniqCounts/Total percentage cutoff")
    parser.add_argument("--mappingStatsPath", default=".", type=Path,
                        help="Path to directory with log files containing bsbolt output")
    parser.add_argument("--trimmingStatsPath", default=".", type=Path,
                        help="Path to directory with log files containing trim_galore output")
    parser.add_argument("--threshold", type=int, help="User provided custom threshold to determine passing cells")
    parser.add_argument("--outDir", default=".", type=Path,
                        help="Path to output directory for plots and report")
    
    args = parser.parse_args()
    return args


def construct_cell_stats_df(maxUniq: int, maxUniqTotal: int, sampleName: str, minUniqTotal: int, i5_barcodes: Path, i7_barcodes: Path,
                            outDir: Path, topCellPercentage: float, minCellRatio: float, minReads: float, tgmt_barcodes: Path, threshold: float) -> pd.DataFrame:
    """
    Construct cell stats for this sample from the cell_stats.tsv file
    Also determine threshold for passing cells based on logic from nf-rna
    """
    cell_stats_complexity_df = pd.read_csv(f"{outDir}/tmp/{sampleName}.cell_stats.csv", header=None,
                                           names=["BC", "total", "passing", "uniq", "MitoReads"])
    cell_stats_complexity_df['percent'] = round(
        (cell_stats_complexity_df['uniq'] / cell_stats_complexity_df['total'] * 100), 3)
    cell_stats_complexity_df['pct_uniq_pass'] = round(
        (cell_stats_complexity_df['uniq'] / cell_stats_complexity_df['passing'] * 100), 3)
    cell_stats_complexity_df['pct_pass_total'] = round(
        (cell_stats_complexity_df['passing'] / cell_stats_complexity_df['total'] * 100), 3)
    cell_stats_complexity_df['pct_pass_cell'] = round(
        (cell_stats_complexity_df['passing'] / cell_stats_complexity_df['passing'].sum() * 100), 3)

    if not maxUniq:
        maxUniq = cell_stats_complexity_df['uniq'].max()

    cell_stats_complexity_df['sampleName'] = sampleName
    cell_stats_complexity_df['pass_filter'] = "fail"
    if not threshold or threshold==0:
        expectedCells = (cell_stats_complexity_df["uniq"] >= minReads).sum()
        threshold = float(np.percentile(cell_stats_complexity_df["uniq"][:expectedCells], (topCellPercentage))/minCellRatio)
    print(f"Threshold: {threshold}", file=sys.stderr)
    cell_stats_complexity_df.loc[
        (cell_stats_complexity_df['percent'] >= minUniqTotal) & (cell_stats_complexity_df['uniq'] >= threshold) &
        (cell_stats_complexity_df['percent'] <= maxUniqTotal) & (cell_stats_complexity_df['uniq'] <= maxUniq),
        "pass_filter"
    ] = "pass"
    cell_stats_complexity_df['threshold'] = threshold
    
    tgmt_barcodes_df = pd.read_csv(tgmt_barcodes, names=["tgmt", "tgmt_well"], sep="\t")
    cell_stats_complexity_df["tgmt"] = cell_stats_complexity_df["BC"].str.split("+").str[2]
    cell_stats_complexity_df = pd.merge(cell_stats_complexity_df, tgmt_barcodes_df, left_on="tgmt", right_on="tgmt", how="left")

    i5_barcodes_df = pd.read_csv(i5_barcodes, names=["i5", "i5_well"], sep="\t")
    cell_stats_complexity_df["i5"] = cell_stats_complexity_df["BC"].str.split("+").str[1]
    cell_stats_complexity_df = pd.merge(cell_stats_complexity_df, i5_barcodes_df, left_on="i5", right_on="i5", how="left")

    i7_barcodes_df = pd.read_csv(i7_barcodes, names=["i7", "i7_well"], sep="\t")
    cell_stats_complexity_df["i7"] = cell_stats_complexity_df["BC"].str.split("+").str[0]
    cell_stats_complexity_df = pd.merge(cell_stats_complexity_df, i7_barcodes_df, left_on="i7", right_on="i7", how="left")
    
    print(f"Number of rows in allCells dataframe: {len(cell_stats_complexity_df.index)}", file=sys.stderr)
    cell_stats_complexity_df.to_csv(f"{outDir}/csv/{sampleName}.allCells.csv", index=False)
    return cell_stats_complexity_df


def construct_met_dataframe(sampleName: str, cell_stats_complexity_df: pd.DataFrame, outDir: Path) -> pd.DataFrame:
    """
    Construct methylation stats for passing cells from the cellInfo.txt file and the cell_stats_complexity_df
    """
    met_df = pd.read_csv(f"{outDir}/tmp/{sampleName}.cellInfo.csv",
                         names=["BC", "Coverage", "CG_Cov", "CG_mC_Pct", "CH_Cov", "CH_mC_Pct"])
    met_df['sampleName'] = sampleName
    passing_complexity_df = cell_stats_complexity_df[cell_stats_complexity_df['pass_filter'] == 'pass']
    met_passing = met_df.merge(passing_complexity_df, on=["BC", "sampleName"])
    met_passing['cg_total_ratio'] = round((met_passing['uniq'] / met_passing['total']), 2)
    print(f"Number of rows in passingCellsMapMethylStats dataframe: {len(met_passing)}", file=sys.stderr)
    met_passing.to_csv(f"{outDir}/csv/{sampleName}.passingCellsMapMethylStats.csv", index=False)
    return met_passing


def filter_cells(sampleName: str, barcodes: Path, outDir: Path):
    """
    Construct barcodes for passing cells on a per well basis
    """
    data = []
    with open(f"{outDir}/csv/{sampleName}.allCells.csv", "r") as infile:
        reader = csv.reader(infile)
        for row in reader:
            if row[10] == "pass":
                cell = row[0].split("+")
                data.append((row[0], cell[2]))

    PV = {}
    with open(barcodes, "r") as bcfile:
        reader = csv.reader(bcfile, delimiter="\t")
        for row in reader:
            PV[row[0]] = row[1]

    Path("filter_cells").mkdir(parents=True, exist_ok=True)
    for row in data:
        if row[1] in PV:
            outfile_name = f"filter_cells/{sampleName}.{PV[row[1]]}.cells"
            with open(outfile_name, "a") as outfile:
                outfile.write(row[0] + "\n")


def build_mapping_stats(mappingStatsPath: Path, sampleName: str, outDir: Path) -> pd.DataFrame | bool:
    """
    Build mapping stats from the bsbolt log files
    """
    suffix = "*bsbolt.log"
    fullPath = os.path.join(mappingStatsPath, sampleName + suffix)
    filenamesList = glob.glob(fullPath)
    if filenamesList == []:
        print("No logs found for bsbolt(mapping)", file=sys.stderr)
        return False
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
    mapping_report.to_csv(f'{outDir}/csv/{sampleName}.mapping_stats.csv', index=False)
    return mapping_report


def build_trimming_stats(trimmingStatsPath: Path, sampleName: str, outDir: Path) -> pd.DataFrame | bool:
    """
    Build trimming stats from trimgalore report files
    """
    suffix = "*.fq.gz_trimming_report_R2.txt"
    fullPath = os.path.join(trimmingStatsPath, sampleName + suffix)
    filenamesList = glob.glob(fullPath)
    if filenamesList == []:
        print("No logs found for trim galore(trimming)", file=sys.stderr)
        return False
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
    trimming_report.to_csv(f'{outDir}/csv/{sampleName}.trimming_stats.csv', index=False)
    return trimming_report


class ConcatPerBarcodeStats:
    """
    Contains function to concatenate files that were generated on a per well basis
    """
    def __init__(self, sampleName: str, writeDir: Path):
        self.sampleName = sampleName
        self.writeDir = writeDir

    def concat_cell_stats(self, cellStatsPath: Path):
        """
        Concatenate cell_stats.tsv files for a sample
        """
        tsv_files = glob.glob(f"{cellStatsPath}/{self.sampleName}*.cell_stats.tsv")
        if tsv_files == []:
            print("No cell_stats.tsv file found", file=sys.stderr)
        with open(f"{self.writeDir}/tmp/{self.sampleName}.cell_stats.csv", "w") as outfile:
            for file in tsv_files:
                with open(file, "r") as infile:
                    lines = infile.readlines()
                    for line in lines:
                        if "Barcode" not in line:
                            outfile.write(line.replace("\t", ","))

    def concat_fragment_histogram(self, fragmentHistPath: Path):
        """
        Concatenate fragment_hist.tsv files for a sample
        """
        file_list = glob.glob(f"{fragmentHistPath}/{self.sampleName}*.fragment_hist.tsv")
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
        with open(f"{self.writeDir}/csv/{self.sampleName}.fragment_hist.csv", 'w') as output:
            for key, value in sorted_data:
                output.write(f"{key},{value}\n")

    def concat_cellinfo(self, cellInfoPath: Path):
        """
        Concatenate cellInfo.txt files for a sample
        """
        file_list = glob.glob(f"{cellInfoPath}/{self.sampleName}*.cellInfo.txt")
        if file_list == []:
            print("No cellInfo.txt file found", file=sys.stderr)
        with open(f"{self.writeDir}/tmp/{self.sampleName}.cellInfo.csv", "w") as outfile:
            for file_name in file_list:
                with open(file_name, "r") as infile:
                    contents = infile.read()
                    # Exclude lines containing "CellID"
                    filtered_contents = [line for line in contents.split("\n") if "CellID" not in line]
                    outfile.write(("\n".join(filtered_contents)).replace("\t", ","))


def main():
    args = accept_args()
    writeDir = args.outDir 
    writeDir.mkdir(parents=True, exist_ok=True)
    Path(writeDir / "csv").mkdir(parents=True, exist_ok=True)
    Path(writeDir / "png").mkdir(parents=True, exist_ok=True)
    Path(writeDir / "tmp").mkdir(parents=True, exist_ok=True)

    concat_stats_obj = ConcatPerBarcodeStats(args.sampleName, writeDir)
    concat_stats_obj.concat_cell_stats(args.cellStatsPath)
    concat_stats_obj.concat_fragment_histogram(args.fragmentHistPath)
    concat_stats_obj.concat_cellinfo(args.cellInfoPath)

    cell_stats_complexity_df = construct_cell_stats_df(
        args.maxUniq, args.maxUniqTotal, args.sampleName, args.minUniqTotal, args.i5_barcodes, args.i7_barcodes,
        args.outDir, args.topCellPercentage, args.minCellRatio, args.minReads, args.tgmt_barcodes, args.threshold)
    met_passing = construct_met_dataframe(args.sampleName, cell_stats_complexity_df, args.outDir)
    mapping_stats = build_mapping_stats(args.mappingStatsPath, args.sampleName, writeDir)
    trimming_stats = build_trimming_stats(args.trimmingStatsPath, args.sampleName, writeDir)

    dp_list_read_qc = []
    dp_list_methyl_qc = []
    dp_page = []
    fragment_df = pd.read_csv(f"{writeDir}/csv/{args.sampleName}.fragment_hist.csv", header=None, names=["bin", "count"])

    # Initialize object for building reads page of the report
    reads_page_obj = BuildReadsPage(args.sampleName, writeDir, cell_stats_complexity_df, False, fragment_df)
    dp_list_read_qc.append(reads_page_obj.build_knee_plot())
    dp_list_read_qc.append(reads_page_obj.build_fragment_length_histogram())
    dp_list_read_qc.append(reads_page_obj.create_total_complexity_plot())
    dp_list_read_qc.append(reads_page_obj.construct_passing_cell_stats(met_passing))
    dp_list_read_qc.append(reads_page_obj.build_summary_stats(mapping_stats, trimming_stats))
    dp_list_read_qc.append(reads_page_obj.build_summary_stats_table(mapping_stats, trimming_stats))
    
    # Initialize object for building methyl page of report
    methyl_page_obj = BuildMethylPage(args.sampleName, writeDir, met_passing, False)
    dp_list_methyl_qc.append(methyl_page_obj.build_cell_covered_box())
    dp_list_methyl_qc.append(methyl_page_obj.build_cg_cell_methyl_percent_box())
    dp_list_methyl_qc.append(methyl_page_obj.build_ch_cell_methyl_percent_box())
    dp_list_methyl_qc.append(methyl_page_obj.build_cell_cg_per_total())
    dp_list_methyl_qc.append(methyl_page_obj.build_total_and_unique_reads_box())
    dp_list_methyl_qc.append(methyl_page_obj.build_uniq_over_total_percent_box())

    dp_page.append(dp.Page(dp.Group(blocks=dp_list_read_qc, columns=2), title="Read QC"))
    dp_page.append(dp.Page(dp.Group(blocks=dp_list_methyl_qc, columns=2), title="Methylation QC"))
    # Build barcodes page of the report
    dp_page.append(dp.Page(dp.Group(blocks=build_plate_plot(cell_stats_complexity_df[cell_stats_complexity_df["pass_filter"]=="pass"], args.tgmt_barcodes, writeDir, args.sampleName), columns=2), title="Barcodes"))
    dp.save_report(dp_page, path=f"{writeDir}/{args.sampleName}.report.html")

    filter_cells(args.sampleName, args.tgmt_barcodes, args.outDir)
    shutil.rmtree(writeDir / "tmp")


if __name__ == "__main__":
    main()
