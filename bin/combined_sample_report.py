#!/usr/bin/env python
"""
Generate combined sample report for all samples in a library
"""
import argparse
import pandas as pd
import datapane as dp
import sys
from pathlib import Path
from reporting.reporting import build_combined_plate_plot
from reporting.reporting import BuildReadsPage
from reporting.reporting import BuildMethylPage


def accept_args() -> argparse.Namespace:
    """
    Accept command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Generate combined sample report")
    parser.add_argument("--libraryName", type=str, required=True,
                        help="Identifier for this library")
    parser.add_argument("--outDir", default=".", type=Path,
                        help="Path to output directory for plots and report")
    parser.add_argument("--complexity_stats", nargs="+", type=Path,
                        help="Path to complexity stats files")
    parser.add_argument("--methyl_stats", nargs="+", type=Path, required=True,
                        help="Path to methyl stats files")
    parser.add_argument("--passing_met_stats", nargs="+", type=Path, required=True,
                        help="Path to passing methyl stats files")
    parser.add_argument("--bcParserMetrics", default="demuxMetrics.json", type=Path,
                        help="Path to json file containing bcParser metrics")
    parser.add_argument("--i5_barcodes", type=Path,
                        help="Path to file containing i5 barcodes")
    parser.add_argument("--i7_barcodes", type=Path,
                        help="Path to file containing i7 barcodes")
    parser.add_argument("--tgmt_barcodes", type=Path,
                        help="Path to file containing tagmentation barcodes")
    args = parser.parse_args()
    return args


def main():
    args = accept_args()
    writeDir = args.outDir
    writeDir.mkdir(parents=True, exist_ok=True)
    Path(writeDir / "png").mkdir()
    Path(writeDir / "csv").mkdir()

    dp_list_methyl_qc = []
    dp_list_read_qc = []
    dp_page = []

    complexity_dfs = [pd.read_csv(complexity_stat) for complexity_stat in args.complexity_stats]
    library_complexity_df = pd.concat(complexity_dfs, ignore_index=True, axis=0)
    print(f"Number of rows in combined allCells dataframe: {len(library_complexity_df.index)}", file=sys.stderr)
    
    met_stats_df = [pd.read_csv(methyl_stat) for methyl_stat in args.passing_met_stats]
    met_stats_df = [df.rename({"Value": df[df["Metric"] == "sample_name"]["Value"].to_list()[0]}, axis=1) for df in met_stats_df]
    met_stats_df = [df.set_index("Metric") for df in met_stats_df]
    library_met_stats_df = pd.concat(met_stats_df, axis=1)
    library_met_stats_df.reset_index(inplace=True)
    print(f"Number of rows in combined passingCellsMapMethylStats dataframe: {len(library_met_stats_df)}", file=sys.stderr)

    # Initialize object for building reads page of the report
    reads_page_obj = BuildReadsPage(f"library.{args.libraryName}", writeDir, library_complexity_df, True)
    dp_list_read_qc.append(reads_page_obj.create_total_complexity_plot())
    dp_list_read_qc.append(reads_page_obj.build_knee_plot())
    if args.bcParserMetrics.is_file():
        dp_list_read_qc.append(reads_page_obj.build_bc_parser_stats(args.bcParserMetrics))
    else:
        print("bc_parser metrics not available", file=sys.stderr)
        dp_list_read_qc.append(reads_page_obj.build_bc_parser_stats(False))
    dp_list_read_qc.append(reads_page_obj.build_threshold_table())
    reads_page_obj.build_combined_passing_met_stats_table(library_met_stats_df)

    methyl_dfs = [pd.read_csv(methyl_stat) for methyl_stat in args.methyl_stats]
    library_methyl_df = pd.concat(methyl_dfs, ignore_index=True, axis=0)
    
    # Initialize object for building methyl page of report
    methyl_page_obj = BuildMethylPage(f"library.{args.libraryName}", writeDir, library_methyl_df, True)
    dp_list_methyl_qc.append(methyl_page_obj.build_cell_covered_box())
    dp_list_methyl_qc.append(methyl_page_obj.build_cg_cell_methyl_percent_box())
    dp_list_methyl_qc.append(methyl_page_obj.build_ch_cell_methyl_percent_box())
    dp_list_methyl_qc.append(methyl_page_obj.build_cell_cg_per_total())
    dp_list_methyl_qc.append(methyl_page_obj.build_total_and_unique_reads_box())
    dp_list_methyl_qc.append(methyl_page_obj.build_uniq_over_total_percent_box())

    dp_list_barcodes = []
    dp_list_barcodes.append(build_combined_plate_plot(library_complexity_df[library_complexity_df["pass_filter"]=="pass"], args.i5_barcodes, "i5", args.tgmt_barcodes, args.libraryName, writeDir))
    dp_list_barcodes.append(build_combined_plate_plot(library_complexity_df[library_complexity_df["pass_filter"]=="pass"], args.i7_barcodes, "i7", args.tgmt_barcodes, args.libraryName, writeDir))
    
    dp_page.append(dp.Page(dp.Group(blocks=dp_list_read_qc, columns=2), title="Read QC"))
    dp_page.append(dp.Page(dp.Group(blocks=dp_list_methyl_qc, columns=2), title="Methylation QC"))
    dp_page.append(dp.Page(dp.Group(blocks=[y for x in dp_list_barcodes for y in x], columns=2), title="Barcodes"))
    dp.save_report(dp_page, path=f"{writeDir}/library.{args.libraryName}.report.html")


if __name__ == "__main__":
    main()
