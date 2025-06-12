#!/usr/bin/env python
"""
Generate combined sample report for all samples in a library
"""
import argparse
import pandas as pd
from pathlib import Path
from reporting.reporting import (
    WhitelistAndLibJsonHelper,
    GenerateMetrics,
)


def accept_args() -> argparse.Namespace:
    """
    Accept command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Generate combined sample report",
    )
    parser.add_argument("--library_name", type=str, required=True, help="Identifier for this library")
    parser.add_argument(
        "--out_dir",
        default=".",
        type=Path,
        help="Path to output directory for plots and report",
    )
    parser.add_argument("--all_cells", nargs="+", type=Path, help="Path to complexity stats files")
    parser.add_argument(
        "--passing_cell_stats",
        nargs="+",
        type=Path,
        required=True,
        help="Path to passing methyl stats files",
    )
    parser.add_argument(
        "--library_structure_json",
        type=Path,
        help="Path to library structure json file",
    )
    parser.add_argument("--workflow_version", type=str, help="Version of the workflow")
    args = parser.parse_args()
    return args


def main():
    args = accept_args()
    write_dir = args.out_dir
    write_dir.mkdir(parents=True, exist_ok=True)
    Path(write_dir / "csv").mkdir(exist_ok=True)

    all_cells = [pd.read_csv(all_cell, dtype={"tgmt_well": "str"}) for all_cell in args.all_cells]
    library_all_cells = pd.concat(all_cells, ignore_index=True, axis=0)
    library_methyl_df = library_all_cells[library_all_cells["pass"] == "pass"]
    library_methyl_df = library_methyl_df.sort_values(by="sample")

    lib_json = WhitelistAndLibJsonHelper.get_lib_json(args.library_structure_json)
    barcode_list = lib_json["barcodes"]
    barcode_name_list = []
    for barcodes in barcode_list:
        if barcodes["name"] != "umi" and barcodes["name"] != lib_json["sample_barcode"]:
            barcode_name_list.append(barcodes["name"])

    GenerateMetrics.build_combined_plate_plot_csv(
        library_methyl_df,
        barcode_name_list,
        args.library_name,
        write_dir,
        args.library_structure_json,
    )

    passing_cell_summary_stat = []
    for file in args.passing_cell_stats:
        passing_cell_file = pd.read_csv(file)
        passing_cell_value = passing_cell_file.rename(
            {"Value": passing_cell_file[passing_cell_file["Metric"] == "sample_name"]["Value"].to_list()[0]}, axis=1
        )
        passing_cell_summary_stat.append(passing_cell_value)

    library_passing_cell_summary_stats = pd.concat(passing_cell_summary_stat, axis=1)
    library_passing_cell_summary_stats.reset_index(inplace=True)
    GenerateMetrics.build_combined_passing_met_stats_table_csv(
        write_dir, args.library_name, library_passing_cell_summary_stats
    )


if __name__ == "__main__":
    main()
