#!/usr/bin/env python
"""
Generate sample report and associated metrics
"""
import argparse
import pandas as pd
import datapane as dp
import os
import numpy as np
import shutil
from pathlib import Path

from reporting.reporting import BuildDatapane, PlatePlotBuilder


def accept_args() -> argparse.Namespace:
    """
    Accept command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Generate sample report")
    parser.add_argument("--library_structure_json", type=Path,
                        help="Path to library structure json file")
    parser.add_argument("--all_cells", type=Path,
                        help="Path to all_cells csv file")
    parser.add_argument("--mapping_stats", type=Path,
                        help="Path to concatenated bsbolt log file")
    parser.add_argument("--trimming_stats", type=Path,
                        help="Path to concatenated trim_galore log file")
    parser.add_argument("--fragment_hist", type=Path,
                        help="Path to concatenated fragment_hist.tsv file")
    parser.add_argument("--dedup_stats", type=Path,
                        help="Path to concatenated dedup_stats.tsv file")
    parser.add_argument("--out_dir", default=".", type=Path,
                        help="Path to output directory for plots and report")
    parser.add_argument("--sample_name", type=str, required=True,
                        help="Identifier for this sample")
    parser.add_argument("--tss_enrich", type=str, nargs="*",
                        help="Path to files containing tss enrichment stats")
    parser.add_argument("--reporting_only", action="store_true", default=False,
                        help="Whether this a reporting only run")
    
    args = parser.parse_args()
    return args


def main():
    args = accept_args()
    writeDir = args.out_dir 
    writeDir.mkdir(parents=True, exist_ok=True)
    Path(writeDir / "csv").mkdir(parents=True, exist_ok=True)
    Path(writeDir / "png").mkdir(parents=True, exist_ok=True)

    # Copy over files that were generated in either generateMetrics or merge_stats process, but we want to put them in the report folder
    shutil.copy(args.dedup_stats, f"{args.sample_name}/csv/{args.dedup_stats.name}")

    if args.mapping_stats and args.trimming_stats:
        mapping_stats = pd.read_csv(args.mapping_stats)
        trimming_stats = pd.read_csv(args.trimming_stats)
        shutil.copy(args.mapping_stats, f"{args.sample_name}/csv/{args.mapping_stats.name}")
        shutil.copy(args.trimming_stats, f"{args.sample_name}/csv/{args.trimming_stats.name}")
    else:
        mapping_stats = False
        trimming_stats = False
    all_cells = pd.read_csv(args.all_cells, dtype={'tgmt_well': 'str'})
    fragment_df = pd.read_csv(args.fragment_hist, header=None, names=["bin", "count"])

    concat_tss_enrich_df = pd.DataFrame()
    if args.reporting_only:
        if args.tss_enrich:
            concat_tss_enrich_df = pd.read_csv(args.tss_enrich[0])
    else:
        if args.tss_enrich:
            tss_enrich_dfs = []
            for tss_enrich in args.tss_enrich:
                # Some files may be empty, so we need to check for that
                if os.stat(tss_enrich).st_size > 1:
                    tss_enrich_dfs.append(pd.read_csv(tss_enrich, sep="\t"))
            if len(tss_enrich_dfs) > 0:
                concat_tss_enrich_df = pd.concat(tss_enrich_dfs, ignore_index=True, axis=0)
                concat_tss_enrich_df = concat_tss_enrich_df.drop(columns=["tss_counts", "background_counts"])
                concat_tss_enrich_df.to_csv(f"{args.sample_name}.tss_enrich.csv", index=False)
    
    if not concat_tss_enrich_df.empty:
        all_cells = all_cells.merge(concat_tss_enrich_df, on=["cell_id"], how="left")
        print(f"Number of NaN values in tss enrichment column: {all_cells['tss_enrich'].isna().sum()}")
        print(f"Number of Inf values in tss enrichment column: {all_cells['tss_enrich'].isin([np.inf]).sum()}")
    met_passing = all_cells[all_cells["pass"]=="pass"].reset_index(drop=True)
    all_cells.to_csv(f"{args.sample_name}.allCells.csv", index=False)

    shutil.copy(args.fragment_hist, f"{args.sample_name}/csv/{args.fragment_hist.name}")
    
    # Sort to keep plotting order deterministic
    # Helps ensure same samples are given same colors every time this script is run
    all_cells = all_cells.sort_values(by='sample')
    met_passing = met_passing.sort_values(by='sample')

    dp_list_read_qc = []
    dp_list_methyl_qc = []
    dp_list_summary = []
    dp_page = []

    datapane_obj = BuildDatapane(args.sample_name, writeDir, all_cells, False, fragment_df)
    datapane_obj.build_summary_stats(mapping_stats, trimming_stats)
    dp_list_summary.append(datapane_obj.build_knee_plot())
    dp_list_summary.append(datapane_obj.create_total_complexity_plot())
    if not met_passing.empty:
        dp_list_summary.append(datapane_obj.construct_passing_cell_stats())
    dp_list_summary.append(datapane_obj.build_summary_stats_table(mapping_stats, trimming_stats))
    dp_page.append(dp.Page(dp.Group(blocks=dp_list_summary, columns=2), title="Summary"))

    dp_list_read_qc.append(datapane_obj.build_total_and_unique_reads_box())
    dp_list_read_qc.append(datapane_obj.build_fragment_length_histogram())
    if not concat_tss_enrich_df.empty:
        dp_list_read_qc.append(datapane_obj.build_tss_enrich_box())
    dp_list_read_qc.append(datapane_obj.build_mito_box())
    dp_page.append(dp.Page(dp.Group(blocks=dp_list_read_qc, columns=2), title="Cell Stats"))

    if not met_passing.empty:
        dp_list_methyl_qc.append(datapane_obj.build_cg_cell_methyl_percent_box())
        dp_list_methyl_qc.append(datapane_obj.build_ch_cell_methyl_percent_box())
        dp_list_methyl_qc.append(datapane_obj.build_cell_covered_box())
        dp_list_methyl_qc.append(datapane_obj.build_mito_tss_ch_table())
        dp_page.append(dp.Page(dp.Group(blocks=dp_list_methyl_qc, columns=2), title="Methylation QC"))
        
        # Build barcodes page of the report
        dp_page.append(dp.Page(dp.Group(blocks=PlatePlotBuilder.build_plate_plot(met_passing, args.library_structure_json, writeDir, args.sample_name),
                                        columns=2), title="Barcodes"))
    dp.save_report(dp_page, path=f"{writeDir}/{args.sample_name}.report.html")


if __name__ == "__main__":
    main()
