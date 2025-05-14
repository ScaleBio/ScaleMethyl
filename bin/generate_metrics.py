#!/usr/bin/env python
"""
Generate sample report and associated metrics
"""
import argparse
import pandas as pd
import sys
import json
import numpy as np
import shutil
from pathlib import Path
from reporting.reporting import GenerateMetrics


def accept_args() -> argparse.Namespace:
    """
    Accept command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Generate sample report",
    )
    parser.add_argument(
        "--library_structure_json",
        type=Path,
        help="Path to library structure json file",
    )
    parser.add_argument(
        "--cell_info",
        type=Path,
        help="Path to directory with files containing methylation calls for each split file",
    )
    parser.add_argument(
        "--cell_stats",
        type=Path,
        help="Path to directory containing split cell_stats.tsv files",
    )
    parser.add_argument("--sample_name", type=str, required=True, help="Identifier for this sample")
    parser.add_argument(
        "--top_cell_percentage",
        default=99.0,
        type=float,
        help="Percentage of cells over min_reads to use as 'robust max'",
    )
    parser.add_argument(
        "--min_cell_ratio",
        default=200.0,
        type=float,
        help="Ratio between transcript counts of top cells and the lower cell threshold",
    )
    parser.add_argument(
        "--min_reads",
        default=1000.0,
        type=float,
        help="Minimum counts to consider a barcode as a potential cell",
    )
    parser.add_argument(
        "--min_uniq_total",
        default=0.0,
        type=float,
        help="Minimum uniqCounts/Total percentage cutoff",
    )
    parser.add_argument(
        "--max_uniq",
        type=float,
        help="Maximum unique counts cutoff (optional; default = max of tsv)",
    )
    parser.add_argument(
        "--max_uniq_total",
        default=100.0,
        type=float,
        help="Maximum uniqCounts/Total percentage cutoff",
    )
    parser.add_argument(
        "--threshold",
        type=int,
        help="User provided custom threshold to determine passing cells",
    )
    parser.add_argument(
        "--out_dir",
        default=".",
        type=Path,
        help="Path to output directory for plots and report",
    )
    parser.add_argument(
        "--mapping_stats",
        type=Path,
        help="Path to mapping stats file",
    ),
    parser.add_argument(
        "--trimming_stats",
        type=Path,
        help="Path to trimming stats file",
    )

    args = parser.parse_args()
    return args


def construct_all_cells(
    max_uniq: int,
    max_uniq_total: int,
    sample_name: str,
    min_uniq_total: int,
    cell_stats: Path,
    out_dir: Path,
    top_cell_percentage: float,
    min_cell_ratio: float,
    min_reads: float,
    threshold: float,
    cell_info: Path,
    library_structure_json: Path,
):
    """
    Construct cell stats for this sample from the cell_stats.tsv file
    Also determine threshold for passing cells based on logic from nf-rna
    """
    met_df = pd.read_csv(
        cell_info,
        names=[
            "cell_id",
            "cov",
            "cg_cov",
            "mcg_pct",
            "ch_cov",
            "mch_pct",
            "ch_high_reads",
        ],
    )
    met_df["sample"] = sample_name
    all_cells = pd.read_csv(cell_stats)
    if len(all_cells) == 0:
        print(f"allCells is empty. No cells found in {cell_stats}", file=sys.stderr)
        sys.exit(1)
    all_cells.rename(
        columns={
            "BC": "cell_id",
            "total": "total_reads",
            "uniq": "unique_reads",
            "passing": "passing_reads",
        },
        inplace=True,
    )
    all_cells["percent"] = round((all_cells["unique_reads"] / all_cells["total_reads"] * 100), 3)
    all_cells["mito_reads"] = round((all_cells["MitoReads"] / all_cells["total_reads"] * 100), 3)
    all_cells.sort_values(by="unique_reads", inplace=True, ascending=False)
    if not max_uniq:
        max_uniq = all_cells["unique_reads"].max()
    all_cells["sample"] = sample_name
    all_cells["pass"] = "fail"
    if not threshold or threshold == 0:
        expected_cells = (all_cells["unique_reads"] >= min_reads).sum()
        if expected_cells == 0:
            threshold = min_reads
        else:
            calc_threshold = float(
                np.percentile(all_cells["unique_reads"][:expected_cells], (top_cell_percentage)) / min_cell_ratio
            )
            # Respect min_reads threshold if calculated threshold is lower
            threshold = max(min_reads, calc_threshold)

    print(f"Threshold: {threshold}", file=sys.stderr)
    all_cells.loc[
        (all_cells["percent"] >= min_uniq_total)
        & (all_cells["unique_reads"] >= threshold)
        & (all_cells["percent"] <= max_uniq_total)
        & (all_cells["unique_reads"] <= max_uniq),
        "pass",
    ] = "pass"

    for barcode_order, barcode_info in enumerate(json.load(open(library_structure_json))["barcodes"]):
        barcodes_df = pd.read_csv(
            library_structure_json.parent / barcode_info["sequences"],
            names=[f"{barcode_info['name']}", f"{barcode_info['name']}_well"],
            sep="\t",
        )
        all_cells[f"{barcode_info['name']}"] = all_cells["cell_id"].str.split("+").str[barcode_order]
        all_cells = pd.merge(
            all_cells,
            barcodes_df,
            left_on=f"{barcode_info['name']}",
            right_on=f"{barcode_info['name']}",
            how="left",
        )
        all_cells = all_cells.drop(columns=[f"{barcode_info['name']}"])

    all_cells = all_cells.merge(met_df, on=["cell_id", "sample"], how="left")
    all_cells["ch_high_reads_percent"] = round((all_cells["ch_high_reads"] / all_cells["total_reads"] * 100), 3)
    all_cells = all_cells.drop(columns=["MitoReads", "percent", "ch_high_reads"])
    all_cells.to_csv(f"{out_dir}/{sample_name}.allCells.csv", index=False)

    print(
        f"Number of passing rows in all_cells: {all_cells[all_cells['pass'] == 'pass'].shape[0]}",
        file=sys.stderr,
    )
    return all_cells


def construct_csvs(
    all_cells: pd.DataFrame,
    library_structure_json: Path,
    write_dir: Path,
    sample: str,
    mapping_stats: Path,
    trimming_stats: Path,
):
    """
    Construct csvs for this sample
    """
    if mapping_stats and trimming_stats:
        mapping_stats_file = pd.read_csv(mapping_stats)
        trimming_stats_file = pd.read_csv(trimming_stats)
        shutil.copy(mapping_stats, f"{write_dir}/csv/{mapping_stats.name}")
        shutil.copy(trimming_stats, f"{write_dir}/csv/{trimming_stats.name}")
    else:
        mapping_stats_file = False
        trimming_stats_file = False
    GenerateMetrics.build_plate_plot_csv(all_cells, library_structure_json, write_dir, sample)
    GenerateMetrics.build_summary_stats_csv(write_dir, sample, mapping_stats_file, trimming_stats_file)
    GenerateMetrics.construct_passing_cell_stats_csv(all_cells, write_dir, sample)


def main():
    args = accept_args()
    write_dir = Path(args.sample_name)
    write_dir.mkdir(parents=True, exist_ok=True)
    Path(write_dir / "csv").mkdir(exist_ok=True)

    all_cells = construct_all_cells(
        args.max_uniq,
        args.max_uniq_total,
        args.sample_name,
        args.min_uniq_total,
        args.cell_stats,
        args.out_dir,
        args.top_cell_percentage,
        args.min_cell_ratio,
        args.min_reads,
        args.threshold,
        args.cell_info,
        args.library_structure_json,
    )
    construct_csvs(
        all_cells,
        args.library_structure_json,
        write_dir,
        args.sample_name,
        args.mapping_stats,
        args.trimming_stats,
    )


if __name__ == "__main__":
    main()
