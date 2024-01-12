#!/usr/bin/env python
"""
Generate sample report and associated metrics
"""
import argparse
import pandas as pd
import plotly.express as px
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
    parser.add_argument("--cellInfo", type=Path,
                        help="Path to directory with files containing methylation calls for each split file")
    parser.add_argument("--cellStats", type=Path,
                        help="Path to directory containing split cell_stats.tsv files")
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
    parser.add_argument("--threshold", type=int, help="User provided custom threshold to determine passing cells")
    parser.add_argument("--outDir", default=".", type=Path,
                        help="Path to output directory for plots and report")
    
    args = parser.parse_args()
    return args


def construct_cell_stats_df(maxUniq: int, maxUniqTotal: int, sampleName: str, minUniqTotal: int, i5_barcodes: Path, i7_barcodes: Path, cellStats: Path,
                            outDir: Path, topCellPercentage: float, minCellRatio: float, minReads: float, tgmt_barcodes: Path, threshold: float) -> pd.DataFrame:
    """
    Construct cell stats for this sample from the cell_stats.tsv file
    Also determine threshold for passing cells based on logic from nf-rna
    """
    cell_stats_complexity_df = pd.read_csv(cellStats, header=None, names=["BC", "total", "passing", "uniq", "MitoReads"])
    cell_stats_complexity_df['percent'] = round(
        (cell_stats_complexity_df['uniq'] / cell_stats_complexity_df['total'] * 100), 3)
    cell_stats_complexity_df['pct_uniq_pass'] = round(
        (cell_stats_complexity_df['uniq'] / cell_stats_complexity_df['passing'] * 100), 3)
    cell_stats_complexity_df['pct_pass_total'] = round(
        (cell_stats_complexity_df['passing'] / cell_stats_complexity_df['total'] * 100), 3)
    cell_stats_complexity_df['pct_pass_cell'] = round(
        (cell_stats_complexity_df['passing'] / cell_stats_complexity_df['passing'].sum() * 100), 3)
    cell_stats_complexity_df['pct_mito'] = round(
        (cell_stats_complexity_df["MitoReads"] / cell_stats_complexity_df['total'] * 100), 3)
    cell_stats_complexity_df.sort_values(by='uniq', inplace=True, ascending=False)
    if not maxUniq:
        maxUniq = cell_stats_complexity_df['uniq'].max()

    cell_stats_complexity_df['sampleName'] = sampleName
    cell_stats_complexity_df['pass_filter'] = "fail"
    if not threshold or threshold==0:
        expectedCells = (cell_stats_complexity_df["uniq"] >= minReads).sum()
        if expectedCells == 0:
            threshold = minReads
        else:
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
    cell_stats_complexity_df.to_csv(f"{outDir}/{sampleName}.allCells.csv", index=False)
    return cell_stats_complexity_df


def construct_met_dataframe(sampleName: str, cell_stats_complexity_df: pd.DataFrame, outDir: Path, cellInfo: Path) -> pd.DataFrame:
    """
    Construct methylation stats for passing cells from the cellInfo.txt file and the cell_stats_complexity_df
    """
    met_df = pd.read_csv(cellInfo, names=["BC", "Coverage", "CG_Cov", "CG_mC_Pct", "CH_Cov", "CH_mC_Pct"])
    met_df['sampleName'] = sampleName
    passing_complexity_df = cell_stats_complexity_df[cell_stats_complexity_df['pass_filter'] == 'pass']
    met_passing = met_df.merge(passing_complexity_df, on=["BC", "sampleName"])
    met_passing['cg_total_ratio'] = round((met_passing['uniq'] / met_passing['total']), 2)
    print(f"Number of rows in passingCellsMapMethylStats dataframe: {len(met_passing)}", file=sys.stderr)
    met_passing.to_csv(f"{outDir}/{sampleName}.passingCellsMapMethylStats.csv", index=False)


def filter_cells(sampleName: str, barcodes: Path, outDir: Path):
    """
    Construct barcodes for passing cells on a per well basis
    """
    data = []
    with open(f"{outDir}/{sampleName}.allCells.csv", "r") as infile:
        reader = csv.reader(infile)
        for row in reader:
            if row[11] == "pass":
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


def main():
    args = accept_args()
    writeDir = args.outDir
    if not writeDir == ".":
        writeDir.mkdir(parents=True, exist_ok=True)
    
    cell_stats_complexity_df = construct_cell_stats_df(
        args.maxUniq, args.maxUniqTotal, args.sampleName, args.minUniqTotal, args.i5_barcodes, args.i7_barcodes, args.cellStats,
        args.outDir, args.topCellPercentage, args.minCellRatio, args.minReads, args.tgmt_barcodes, args.threshold)
    construct_met_dataframe(args.sampleName, cell_stats_complexity_df, args.outDir, args.cellInfo)

    filter_cells(args.sampleName, args.tgmt_barcodes, args.outDir)

if __name__ == "__main__":
    main()