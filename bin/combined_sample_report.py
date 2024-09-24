#!/usr/bin/env python
"""
Generate combined sample report for all samples in a library
"""
import argparse
import pandas as pd
import datapane as dp
import sys
from pathlib import Path
from reporting.reporting import BuildDatapane, WhitelistAndLibJsonHelper, PlatePlotBuilder


def accept_args() -> argparse.Namespace:
    """
    Accept command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Generate combined sample report")
    parser.add_argument("--library_name", type=str, required=True,
                        help="Identifier for this library")
    parser.add_argument("--out_dir", default=".", type=Path,
                        help="Path to output directory for plots and report")
    parser.add_argument("--all_cells", nargs="+", type=Path,
                        help="Path to complexity stats files")
    parser.add_argument("--passing_cell_stats", nargs="+", type=Path, required=True,
                        help="Path to passing methyl stats files")
    parser.add_argument("--bc_parser_metrics", default="demuxMetrics.json", type=Path,
                        help="Path to json file containing bcParser metrics")
    parser.add_argument("--library_structure_json", type=Path,
                        help="Path to library structure json file")
    args = parser.parse_args()
    return args


def main():
    args = accept_args()
    write_dir = args.out_dir
    write_dir.mkdir(parents=True, exist_ok=True)
    Path(write_dir / "png").mkdir(exist_ok=True)
    Path(write_dir / "csv").mkdir(exist_ok=True)

    dp_list_methyl_qc = []
    dp_list_read_qc = []
    dp_page = []

    all_cells = [pd.read_csv(all_cell, dtype={'tgmt_well': 'str'}) for all_cell in args.all_cells]
    library_all_cells = pd.concat(all_cells, ignore_index=True, axis=0)
    print(f"Number of rows in combined allCells dataframe: {len(library_all_cells.index)}", file=sys.stderr)
    
    passing_cell_summary_stat = [pd.read_csv(methyl_stat) for methyl_stat in args.passing_cell_stats]
    passing_cell_summary_stat = [df.rename({"Value": df[df["Metric"] == "sample_name"]["Value"].to_list()[0]}, axis=1) for df in passing_cell_summary_stat]
    passing_cell_summary_stat = [df.set_index("Metric") for df in passing_cell_summary_stat]
    library_passing_cell_summary_stats = pd.concat(passing_cell_summary_stat, axis=1)
    library_passing_cell_summary_stats.reset_index(inplace=True)

    library_all_cells = library_all_cells.sort_values(by='sample')

    datapane_obj = BuildDatapane(f"library.{args.library_name}", write_dir, library_all_cells, True)
    dp_list_read_qc.append(datapane_obj.create_total_complexity_plot())
    dp_list_read_qc.append(datapane_obj.build_knee_plot())
    if args.bc_parser_metrics.is_file():
        dp_list_read_qc.append(datapane_obj.build_bc_parser_stats(args.bc_parser_metrics))
    else:
        print("bc_parser metrics not available", file=sys.stderr)
        dp_list_read_qc.append(datapane_obj.build_bc_parser_stats(False))
    dp_list_read_qc.append(datapane_obj.build_cell_table())
    dp_list_read_qc.append(datapane_obj.build_total_and_unique_reads_box())
    dp_list_read_qc.append(datapane_obj.build_saturation_box())
    datapane_obj.build_combined_passing_met_stats_table(library_passing_cell_summary_stats)
    dp_page.append(dp.Page(dp.Group(blocks=dp_list_read_qc, columns=2), title="Read QC"))

    library_methyl_df = library_all_cells[library_all_cells["pass"] == "pass"]
    library_methyl_df = library_methyl_df.sort_values(by='sample')
    
    if not library_methyl_df.empty:
        dp_list_methyl_qc.append(datapane_obj.build_cg_cell_methyl_percent_box())
        dp_list_methyl_qc.append(datapane_obj.build_ch_cell_methyl_percent_box())
        dp_list_methyl_qc.append(datapane_obj.build_cell_covered_box())
        dp_page.append(dp.Page(dp.Group(blocks=dp_list_methyl_qc, columns=2), title="Methylation QC"))

    dp_list_barcodes = []

    lib_json = WhitelistAndLibJsonHelper.get_lib_json(args.library_structure_json)
    for barcodes in lib_json["barcodes"]:
        if barcodes["name"] != "umi" and barcodes["name"] != lib_json["sample_barcode"]:
            dp_list_barcodes.append(
                PlatePlotBuilder.build_combined_plate_plot(
                    library_methyl_df, barcodes['name'], args.library_name, write_dir, args.library_structure_json
                )
            )
    dp_page.append(dp.Page(dp.Group(blocks=[y for x in dp_list_barcodes for y in x], columns=2), title="Barcodes"))
    
    dp.save_report(dp_page, path=f"{write_dir}/library.{args.library_name}.report.html")


if __name__ == "__main__":
    main()
