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
from reporting.reporting import BuildReadsPage, BuildMethylPage, build_plate_plot


def accept_args() -> argparse.Namespace:
    """
    Accept command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Generate sample report")
    parser.add_argument("--tgmt_barcodes", type=Path,
                        help="Path to file containing tgmt barcodes")
    parser.add_argument("--allCells", type=Path,
                        help="Path to allCells csv file")
    parser.add_argument("--passingCellsMapMethylStats", type=Path,
                        help="Path to passing cells methyl stats csv file")
    parser.add_argument("--mappingStats", type=Path,
                        help="Path to concatenated bsbolt log file")
    parser.add_argument("--trimmingStats", type=Path,
                        help="Path to concatenated trim_galore log file")
    parser.add_argument("--fragmentHist", type=Path,
                        help="Path to concatenated fragment_hist.tsv file")
    parser.add_argument("--dedupStats", type=Path,
                        help="Path to concatenated dedup_stats.tsv file")
    parser.add_argument("--outDir", default=".", type=Path,
                        help="Path to output directory for plots and report")
    parser.add_argument("--sampleName", type=str, required=True,
                        help="Identifier for this sample")
    parser.add_argument("--tssEnrich", type=str, nargs="*",
                        help="Path to files containing tss enrichment stats")
    parser.add_argument("--reporting_only", action="store_true", default=False,
                        help="Whether this a reporting only run")
    
    args = parser.parse_args()
    return args


def main():
    args = accept_args()
    writeDir = args.outDir 
    writeDir.mkdir(parents=True, exist_ok=True)
    Path(writeDir / "csv").mkdir(parents=True, exist_ok=True)
    Path(writeDir / "png").mkdir(parents=True, exist_ok=True)
    shutil.copy(args.dedupStats, f"{args.sampleName}/csv/{args.dedupStats.name}")

    if args.mappingStats and args.trimmingStats:
        mapping_stats = pd.read_csv(args.mappingStats)
        trimming_stats = pd.read_csv(args.trimmingStats)
        shutil.copy(args.mappingStats, f"{args.sampleName}/csv/{args.mappingStats.name}")
        shutil.copy(args.trimmingStats, f"{args.sampleName}/csv/{args.trimmingStats.name}")
    else:
        mapping_stats = False
        trimming_stats = False
    cell_stats_complexity_df = pd.read_csv(args.allCells, dtype={'tgmt_well': 'str'})
    met_passing = pd.read_csv(args.passingCellsMapMethylStats)
    fragment_df = pd.read_csv(args.fragmentHist, header=None, names=["bin", "count"])

    concat_tss_enrich_df = pd.DataFrame()
    if args.reporting_only:
        if args.tssEnrich:
            concat_tss_enrich_df = pd.read_csv(args.tssEnrich[0])
    else:
        if args.tssEnrich:
            tss_enrich_dfs = []
            for tss_enrich in args.tssEnrich:
                if os.stat(tss_enrich).st_size > 1:
                    tss_enrich_dfs.append(pd.read_csv(tss_enrich, sep="\t"))
            if len(tss_enrich_dfs) > 0:
                concat_tss_enrich_df = pd.concat(tss_enrich_dfs, ignore_index=True, axis=0)
                concat_tss_enrich_df.to_csv(f"{args.sampleName}.mergedTssEnrich.csv")
    if not concat_tss_enrich_df.empty:
        met_passing = pd.merge(met_passing, concat_tss_enrich_df, left_on="BC", right_on="BC", how="left")
        print(f"Number of NaN values in tss enrichment column: {met_passing['tss_enrich'].isna().sum()}")
        print(f"Number of Inf values in tss enrichment column: {met_passing['tss_enrich'].isin([np.inf]).sum()}")
    met_passing.to_csv(f"{args.sampleName}/csv/{args.passingCellsMapMethylStats.name}", index=False)
    
    # Copy over files that were generated in either generateMetrics or merge_stats process, but we want to put them in the report folder
    shutil.copy(args.allCells, f"{args.sampleName}/csv/{args.allCells.name}")
    shutil.copy(args.fragmentHist, f"{args.sampleName}/csv/{args.fragmentHist.name}")
    
    cell_stats_complexity_df = cell_stats_complexity_df.sort_values(by='sampleName')
    met_passing = met_passing.sort_values(by='sampleName')

    dp_list_read_qc = []
    dp_list_methyl_qc = []
    dp_page = []

    # Initialize object for building reads page of the report
    reads_page_obj = BuildReadsPage(args.sampleName, writeDir, cell_stats_complexity_df, False, fragment_df)
    dp_list_read_qc.append(reads_page_obj.build_knee_plot())
    dp_list_read_qc.append(reads_page_obj.create_total_complexity_plot())
    dp_list_read_qc.append(reads_page_obj.create_unique_over_total_plot())
    dp_list_read_qc.append(reads_page_obj.create_unique_over_passing_plot())
    dp_list_read_qc.append(reads_page_obj.build_fragment_length_histogram())
    dp_list_read_qc.append(reads_page_obj.build_summary_stats(mapping_stats, trimming_stats))
    dp_list_read_qc.append(reads_page_obj.construct_passing_cell_stats(met_passing))
    dp_list_read_qc.append(reads_page_obj.build_summary_stats_table(mapping_stats, trimming_stats, met_passing))
    dp_page.append(dp.Page(dp.Group(blocks=dp_list_read_qc, columns=2), title="Read QC"))

    if not met_passing.empty:
        # Initialize object for building methyl page of report
        methyl_page_obj = BuildMethylPage(args.sampleName, writeDir, met_passing, False)
        dp_list_methyl_qc.append(methyl_page_obj.build_cell_covered_box())
        dp_list_methyl_qc.append(methyl_page_obj.build_cg_cell_methyl_percent_box())
        dp_list_methyl_qc.append(methyl_page_obj.build_ch_cell_methyl_percent_box())
        dp_list_methyl_qc.append(methyl_page_obj.build_cell_cg_per_total())
        dp_list_methyl_qc.append(methyl_page_obj.build_total_and_unique_reads_box())
        dp_list_methyl_qc.append(methyl_page_obj.build_uniq_over_total_percent_box())
        if not concat_tss_enrich_df.empty:
            dp_list_methyl_qc.append(methyl_page_obj.build_tss_enrich_box())
        dp_list_methyl_qc.append(methyl_page_obj.build_mito_box(cell_stats_complexity_df))
        dp_list_methyl_qc.append(methyl_page_obj.build_mito_table(cell_stats_complexity_df))
        dp_page.append(dp.Page(dp.Group(blocks=dp_list_methyl_qc, columns=2), title="Methylation QC"))

    if not cell_stats_complexity_df[cell_stats_complexity_df["pass_filter"]=="pass"].empty:
        # Build barcodes page of the report
        dp_page.append(dp.Page(dp.Group(blocks=build_plate_plot(cell_stats_complexity_df[cell_stats_complexity_df["pass_filter"]=="pass"],
                                                                args.tgmt_barcodes, writeDir, args.sampleName), columns=2), title="Barcodes"))
    dp.save_report(dp_page, path=f"{writeDir}/{args.sampleName}.report.html")


if __name__ == "__main__":
    main()
