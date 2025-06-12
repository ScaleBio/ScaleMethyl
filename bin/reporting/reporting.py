import pandas as pd
import plotly.express as px
import datapane as dp
import warnings
import json
import sys
import csv
import numpy as np
import plotly.graph_objects as go
import matplotlib.colors as colors
import seaborn as sns
from pathlib import Path
from matplotlib import pyplot as plt
import re
import gzip
import scipy.io
import scanpy as sc

# Suppress SettingWithCopyWarning
pd.options.mode.chained_assignment = None

CELL_STAT_COL_NAMES = [
    "Number of Passing Cells",
    "Total Reads in Cells",
    "Total Passing Reads",
    "Median Passing Reads per Cell",
    "Total Unique Reads",
    "Median Unique Reads per Cell",
    "Percent Reads in Passing Cells",
    "Overall Saturation",
    "Median Mitochondrial Percentage per Cell",
    "Median CG's Covered per Cell",
    "Median CH's Covered per Cell",
    "Median CG Methylation per Cell",
    "Median CH Methylation per Cell",
]


class GenerateMetrics:
    # builds the combined passing met stats table csv
    def build_combined_passing_met_stats_table_csv(out_dir, sample, combined_passing_met: pd.DataFrame):
        for idx, col_name in enumerate(CELL_STAT_COL_NAMES):
            combined_passing_met.iloc[idx, combined_passing_met.columns.get_loc("Metric")] = col_name
        if "median_tss_enrich" in combined_passing_met["Metric"].values:
            combined_passing_met.iloc[
                combined_passing_met["Metric"].values == "median_tss_enrich",
                combined_passing_met.columns.get_loc("Metric"),
            ] = "Median TSS Enrichment"
        combined_passing_met = combined_passing_met.iloc[combined_passing_met["Metric"].values != "sample_name",]
        combined_passing_met.to_csv(f"{out_dir}/{sample}.combinedPassingCellStats.csv", index=False)

    # Update well plate plot matrix with values from df
    def update_well_plate(df, well_plate, value_col):
        for _, row in df.iterrows():
            match = re.match(r"([A-Z])(\d{1,2})([A-Z])(\d{1,2})", row["combo"])
            if match:
                letter_1, number_1, letter_2, number_2 = match.groups()
                combo_1 = letter_1 + number_1
                combo_2 = letter_2 + number_2
                well_plate.at[combo_2, combo_1] = row[value_col]

    # builds the combined plate plot csvs
    def build_combined_plate_plot_csv(
        passing_cells: pd.DataFrame,
        indexing_type_list: list,
        lib_name: str,
        write_dir: Path,
        lib_json: Path,
    ):
        barcode_entry_cols = WhitelistAndLibJsonHelper.get_barcode_entry_from_libjson(lib_json, indexing_type_list[0])
        barcode_entry_index = WhitelistAndLibJsonHelper.get_barcode_entry_from_libjson(lib_json, indexing_type_list[1])
        sequence_path_cols = lib_json.parent / barcode_entry_cols["sequences"]
        sequence_path_index = lib_json.parent / barcode_entry_index["sequences"]
        max_letter_cols, max_number_cols = WhitelistAndLibJsonHelper.get_max_number_and_letter(
            sequence_path_cols, False
        )
        max_letter_index, max_number_index = WhitelistAndLibJsonHelper.get_max_number_and_letter(
            sequence_path_index, False
        )
        all_combos_cols = WhitelistAndLibJsonHelper.get_all_number_letter_combinations(
            max_letter_cols, max_number_cols, True
        )
        all_combos_index = WhitelistAndLibJsonHelper.get_all_number_letter_combinations(
            max_letter_index, max_number_index, False
        )

        well_plate_total_reads = pd.DataFrame(0, index=all_combos_index, columns=all_combos_cols)
        well_plate_total_cells = well_plate_total_reads.copy()
        well_plate_unique_reads = well_plate_total_reads.copy()

        df_by_plate = passing_cells[
            [
                "total_reads",
                "unique_reads",
                f"{indexing_type_list[1]}_well",
                f"{indexing_type_list[0]}_well",
            ]
        ]
        df_by_plate["combo"] = (
            df_by_plate[f"{indexing_type_list[0]}_well"] + df_by_plate[f"{indexing_type_list[1]}_well"]
        )
        num_cells_dict = df_by_plate["combo"].value_counts().to_dict()
        df_unique_reads = df_by_plate.groupby("combo", as_index=False)["unique_reads"].median()

        GenerateMetrics.update_well_plate(df_unique_reads, well_plate_unique_reads, "unique_reads")
        GenerateMetrics.update_well_plate(
            pd.DataFrame.from_dict(num_cells_dict, orient="index", columns=["count"])
            .reset_index()
            .rename(columns={"index": "combo"}),
            well_plate_total_cells,
            "count",
        )

        well_plate_total_cells.to_csv(
            write_dir / "csv" / f"{lib_name}.total_passing_cells_{indexing_type_list[1]}_{indexing_type_list[0]}.csv"
        )
        well_plate_unique_reads.to_csv(
            write_dir / "csv" / f"{lib_name}.median_unique_reads_{indexing_type_list[1]}_{indexing_type_list[0]}.csv"
        )

    # builds the plate plot csvs
    def build_plate_plot_csv(all_cells: pd.DataFrame, lib_json: Path, write_dir: Path, sample: str):
        sample_barcode = WhitelistAndLibJsonHelper.get_sample_barcode_name(lib_json)
        barcode_entry = WhitelistAndLibJsonHelper.get_barcode_entry_from_libjson(lib_json, sample_barcode)
        sequence_path = lib_json.parent / barcode_entry["sequences"]
        max_letter, max_number = WhitelistAndLibJsonHelper.get_max_number_and_letter(sequence_path, True)
        for i in range(1, WhitelistAndLibJsonHelper.get_number_of_plates(sequence_path) + 1):
            well_plate_total_reads = pd.DataFrame(
                0,
                columns=range(1, max_number + 1),
                index=[chr(j) for j in range(65, ord(max_letter) + 1)],
            )
            well_plate_median_unique_reads = well_plate_total_reads.copy()
            well_plate_total_cells = well_plate_total_reads.copy()
            passing_cells = all_cells[all_cells["pass"] == "pass"]
            df_by_plate = passing_cells[passing_cells[f"{sample_barcode}_well"].str[0] == str(i)]
            df_by_plate = df_by_plate[["unique_reads", f"{sample_barcode}_well"]]
            df_by_plate[f"{sample_barcode}_well"] = df_by_plate[f"{sample_barcode}_well"].str[1:]
            num_cells_dict = df_by_plate[f"{sample_barcode}_well"].value_counts().to_dict()
            df_unique_reads = df_by_plate.groupby(f"{sample_barcode}_well", as_index=False)["unique_reads"].median()

            for _, row in df_unique_reads.iterrows():
                letter = row[f"{sample_barcode}_well"][0]
                number = int(row[f"{sample_barcode}_well"][1:])
                well_plate_median_unique_reads.at[letter, number] = row["unique_reads"]

            for well in num_cells_dict:
                letter = well[0]
                number = int(well[1:])
                well_plate_total_cells.at[letter, number] = num_cells_dict[well]

            well_plate_total_cells.to_csv(
                write_dir / "csv" / f"{sample}.total_passing_cells_{sample_barcode}_plate_{i}.csv"
            )
            well_plate_median_unique_reads.to_csv(
                write_dir / "csv" / f"{sample}.median_unique_reads_{sample_barcode}_plate_{i}.csv"
            )

    # builds the summary stats csv
    def build_summary_stats_csv(
        out_dir, sample, mapping_stats: pd.DataFrame | bool, trimming_stats: pd.DataFrame | bool
    ):
        total_read_value = np.nan
        reads_passing_trimming_value = np.nan
        reads_passing_mapping_value = np.nan
        if not (isinstance(mapping_stats, bool) and isinstance(trimming_stats, bool)):
            total_read_value = int(trimming_stats[trimming_stats["Metric"] == "total_reads"]["Value"].to_list()[0])
            reads_passing_trimming_value = trimming_stats[trimming_stats["Metric"] == "passing_reads"][
                "Value"
            ].to_list()[0]
            reads_passing_mapping_value = mapping_stats[mapping_stats["Metric"] == "mapped_reads"]["Value"].to_list()[0]

        summary_stats = pd.DataFrame.from_dict(
            {
                "Metric": [
                    "total_reads",
                    "reads_passing_trimming",
                    "reads_passing_mapping",
                ],
                "Value": [
                    total_read_value,
                    reads_passing_trimming_value,
                    reads_passing_mapping_value,
                ],
            },
            dtype="object",
        )
        summary_stats.to_csv(f"{out_dir}/csv/{sample}.summaryStats.csv", index=False)

    # builds the passing cell stats csv
    def construct_passing_cell_stats_csv(all_cells, out_dir, sample):
        metric_name = [
            "number_passing_cells",
            "total_reads",
            "total_passing_reads",
            "median_passing_reads_per_passing_cell",
            "total_unique_reads",
            "median_unique_reads_per_passing_cell",
            "percent_reads_in_passing_cells",
            "saturation",
            "median_mito_pct",
            "CGcov_median",
            "CHcov_median",
            "CG_mC_Pct_median",
            "CH_mC_Pct_median",
            "sample_name",
        ]
        met_df = all_cells[all_cells["pass"] == "pass"].reset_index(drop=True)
        metric_values = []
        ch_is_nan = pd.isna(met_df["ch_cov"].median())
        passing_cells = all_cells[all_cells["pass"] == "pass"]
        metric_values = [
            len(met_df),
            all_cells["total_reads"].sum(),
            all_cells["passing_reads"].sum(),
            round(passing_cells["passing_reads"].median()),
            all_cells["unique_reads"].sum(),
            round(passing_cells["unique_reads"].median()),
            str(
                round(
                    100 * all_cells[all_cells["pass"] == "pass"]["total_reads"].sum() / all_cells["total_reads"].sum(),
                    2,
                )
            )
            + "%",
            str(
                round(
                    1 - (all_cells["unique_reads"].sum() / all_cells["total_reads"].sum()),
                    3,
                )
            ),
            str(round(met_df["mito_reads"].median() * 100, 2)) + "%",
            round(met_df["cg_cov"].median()),
            round(met_df["ch_cov"].median()) if not ch_is_nan else np.nan,
            str(round(met_df["mcg_pct"].median(), 2)) + "%",
            (str(round(met_df["mch_pct"].median(), 2)) + "%" if not ch_is_nan else np.nan),
            sample,
        ]
        if "tss_enrich" in met_df.columns:
            metric_name.append("median_tss_enrich")
            metric_values.append(str(round(met_df["tss_enrich"].median(), 2)))

        if met_df.empty:
            met_stats = pd.DataFrame.from_dict(
                {"Metric": metric_name, "Value": [np.nan] * len(metric_name)},
                dtype="object",
            )
        else:
            met_stats = pd.DataFrame.from_dict({"Metric": metric_name, "Value": metric_values}, dtype="object")
        met_stats.to_csv(f"{out_dir}/csv/{sample}.passingCellSummaryStats.csv", index=False)


class Utils:
    def get_passing_cells(all_cells: list, sample: str) -> list[str]:
        """
        Get a list of cell barcodes that passed filters

        Args:
            all_cells: Path to the file containing all cell metrics
            sample: Used for getting the well coordinate
        """
        passing_cells = []
        for all_cells_file in all_cells:
            with open(all_cells_file) as csv_file:
                csv_reader = csv.DictReader(csv_file)
                for row in csv_reader:
                    if "." in sample:
                        if row["pass"] == "pass" and row["tgmt_well"] == sample.split(".")[1]:
                            passing_cells.append(row["cell_id"])
                    else:
                        if row["pass"] == "pass" and row["sample"] == sample:
                            passing_cells.append(row["cell_id"])

        if len(passing_cells) == 0:
            print(f"No passing cells found for {sample}")
            sys.exit(0)
        return passing_cells


class WhitelistAndLibJsonHelper:
    @staticmethod
    def get_barcode_entry_from_libjson(library_structure_json: Path, name: str) -> dict:
        with open(library_structure_json) as f:
            lib_json = json.load(f)
        for entry in lib_json["barcodes"]:
            if entry["name"] == name:
                return entry

    @staticmethod
    def get_sample_barcode_name(library_structure_json: Path) -> str:
        with open(library_structure_json) as f:
            lib_json = json.load(f)
        return lib_json["sample_barcode"]

    @staticmethod
    def get_barcode_whitelist_dataframe_from_file(barcodes: Path) -> pd.DataFrame:
        return pd.read_csv(barcodes, sep="\t", names=["barcode", "well"])

    @staticmethod
    def get_barcode_whitelist_dataframe_from_libjson(library_structure_json: Path, name: str) -> pd.DataFrame:
        entry = WhitelistAndLibJsonHelper.get_barcode_entry_from_libjson(library_structure_json, name)
        return pd.read_csv(entry["sequences"], sep="\t", names=["barcode", "well"])

    @staticmethod
    def get_all_barcode_names(library_structure_json: Path) -> list:
        with open(library_structure_json) as f:
            lib_json = json.load(f)
        return [entry["name"] for entry in lib_json["barcodes"]]

    @staticmethod
    def get_all_number_letter_combinations(max_letter: str, max_number: int, order_by_numbers: bool) -> list:

        all_combos = []
        if order_by_numbers:
            for letter in range(65, ord(max_letter) + 1):
                for number in range(1, max_number + 1):
                    all_combos.append(chr(letter) + f"{number:02d}")
        else:
            for number in range(1, max_number + 1):
                for letter in range(65, ord(max_letter) + 1):
                    all_combos.append(chr(letter) + f"{number:02d}")
        return all_combos

    @staticmethod
    def get_max_number_and_letter(barcodes: Path, tgmt: bool = False):
        barcodes_df = pd.read_csv(barcodes, sep="\t", names=["barcode", "well"])
        # Last line corresponds to maximum well coordinate
        # First character denotes plate number so removing that
        max_well_coordinate = ""
        if tgmt:
            max_well_coordinate = barcodes_df.iloc[-1, 1][1:]
        else:
            max_well_coordinate = barcodes_df.iloc[-1, 1]
        max_letter = max_well_coordinate[0]
        max_number = int(max_well_coordinate[1:])
        return max_letter, max_number

    @staticmethod
    def get_number_of_plates(barcodes: Path):
        barcodes_df = pd.read_csv(barcodes, sep="\t", names=["barcode", "well"])
        # Last line corresponds to maximum well coordinate
        # First character denotes plate number
        number_of_plates = barcodes_df.iloc[-1]["well"][0]
        return int(number_of_plates)

    @staticmethod
    def get_lib_json(library_structure_json: Path):
        with open(library_structure_json) as f:
            lib_json = json.load(f)
        return lib_json


class PlatePlotBuilder:
    def build_combined_plate_plot(
        passing_cells: pd.DataFrame,
        indexing_type_list: list,
        csv_folder: Path,
        lib_name: str,
    ) -> list:

        datapane_list = []
        well_plate_total_cells = pd.read_csv(
            f"{csv_folder}/{lib_name}.total_passing_cells_{indexing_type_list[1]}_{indexing_type_list[0]}.csv",
            index_col=0,
            header=0,
        )
        well_plate_unique_reads = pd.read_csv(
            f"{csv_folder}/{lib_name}.median_unique_reads_{indexing_type_list[1]}_{indexing_type_list[0]}.csv",
            index_col=0,
            header=0,
        )
        datapane_list.append(
            PlatePlotBuilder.plot_plate_plot(
                well_plate_total_cells,
                f"Total Passing Cells({indexing_type_list[1]},{indexing_type_list[0]})",
                False,
            )
        )
        datapane_list.append(
            PlatePlotBuilder.plot_plate_plot(
                well_plate_unique_reads,
                f"Median Unique Reads per Passing Cell({indexing_type_list[1]},{indexing_type_list[0]})",
                100,
            )
        )

        return datapane_list

    def build_plate_plot(all_cells: pd.DataFrame, lib_json: Path, csv_folder: Path, sample: str) -> list:
        datapane_list = []
        sample_barcode = WhitelistAndLibJsonHelper.get_sample_barcode_name(lib_json)
        barcode_entry = WhitelistAndLibJsonHelper.get_barcode_entry_from_libjson(lib_json, sample_barcode)
        sequence_path = lib_json.parent / barcode_entry["sequences"]
        for i in range(1, WhitelistAndLibJsonHelper.get_number_of_plates(sequence_path) + 1):
            well_plate_total_cells = pd.read_csv(
                f"{csv_folder}/{sample}.total_passing_cells_tgmt_plate_{str(i)}.csv", index_col=0, header=0
            )
            well_plate_median_unique_reads = pd.read_csv(
                f"{csv_folder}/{sample}.median_unique_reads_tgmt_plate_{str(i)}.csv", index_col=0, header=0
            )

            datapane_list.append(
                PlatePlotBuilder.plot_plate_plot(
                    well_plate_total_cells,
                    f"Total Passing Cells({sample_barcode}): Plate {i}",
                    False,
                )
            )
            datapane_list.append(
                PlatePlotBuilder.plot_plate_plot(
                    well_plate_median_unique_reads,
                    f"Median Unique Reads per Passing Cell({sample_barcode}): Plate {i}",
                    100,
                )
            )

        return datapane_list

    def plot_plate_plot(df: pd.DataFrame, title: str, threshold: int | bool):
        fig = plt.figure()
        if threshold:
            norm = colors.SymLogNorm(linthresh=threshold, vmin=0, vmax=max(threshold * 10, df.values.max()))
            ax = sns.heatmap(
                df,
                linewidth=0.1,
                cmap=sns.color_palette("mako", as_cmap=True),
                norm=norm,
                xticklabels=True,
                yticklabels=True,
            )
        else:
            if (df == 0).all().all():
                norm = colors.SymLogNorm(linthresh=1, vmin=0, vmax=max(1 * 10, df.values.max()))
                ax = sns.heatmap(
                    df,
                    linewidth=0.1,
                    cmap=sns.color_palette("mako", as_cmap=True),
                    norm=norm,
                    xticklabels=True,
                    yticklabels=True,
                )
            else:
                ax = sns.heatmap(
                    df,
                    linewidth=0.1,
                    cmap=sns.color_palette("mako", as_cmap=True),
                    vmin=0.0,
                    vmax=df.values.max(),
                    xticklabels=True,
                    yticklabels=True,
                )
        ax.set_title(title)
        return dp.Plot(fig)


class DatapaneUtils:
    """
    Contains util functions required for generating a datapane report. Taken from nf-rna
    """

    @staticmethod
    def createTableIgnoreWarning(styler, label=None) -> dp.Table:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")

            return dp.Table(styler, label=label)

    @staticmethod
    def reformatIfInt(val):
        if isinstance(val, int):
            return f"{val:,}"
        elif isinstance(val, float):
            return round(val, 2)
        else:
            return val

    @staticmethod
    def styleTable(styler, title: str, hideColumnHeaders=False, boldColumn=None, numericCols=None):
        if numericCols is not None:
            styler.format(DatapaneUtils.reformatIfInt, subset=numericCols)

        # Format NA values
        styler.format(na_rep="NA")

        styler.hide(axis="index")

        if hideColumnHeaders:
            styler.hide(axis="columns")
        else:
            styler.set_table_styles(
                [{"selector": "th", "props": [("border", "black solid !important")]}],
                overwrite=False,
            )

        if boldColumn is not None:
            styler.set_properties(subset=boldColumn, **{"font-weight": "bold"})

        if title != "":
            styler.set_caption(title)
            styler.set_table_styles(
                [
                    {
                        "selector": "caption",
                        "props": [
                            ("border", "black solid !important"),
                            ("font-weight", "bold"),
                        ],
                    }
                ],
                overwrite=False,
            )

        styler.set_properties(**{"border-color": "black", "border-style": "solid !important"})
        return styler

    @staticmethod
    def define_color_map(df: pd.DataFrame) -> dict:
        color_palette = px.colors.qualitative.Dark24
        return {sample: color_palette[i] for i, sample in enumerate(sorted(df["sample"].unique()))}


class BuildDatapane:
    """
    Contains functions for building the datapane html report
    """

    def __init__(
        self,
        sample: str,
        write_dir: Path,
        all_cells: pd.DataFrame,
        is_library_report: bool,
        fragment_df: pd.DataFrame = pd.DataFrame(),
    ):
        self.all_cells = all_cells
        # sample = libName if this is a combined sample report
        self.sample = sample
        self.out_dir = write_dir
        self.fragment_df = fragment_df
        self.is_library_report = is_library_report
        self.met_df = all_cells[all_cells["pass"] == "pass"].reset_index(drop=True)
        self.color_map = DatapaneUtils.define_color_map(all_cells)

    def build_umap_plot(self, mtxFile: str, met_passing: list, writeDir: str) -> dp.Plot:

        with gzip.open(mtxFile, "rb") as f:
            matrix = scipy.io.mmread(f)
        adata = sc.AnnData(matrix.transpose().todense())

        variances = np.var(adata.X, axis=0)
        top_20000_indices = np.argsort(variances)[-20000:]
        adata = adata[:, top_20000_indices]
        # Perform nearest neighbor clustering
        # sqrt(indices/2) to get the number of comps
        knn = min(100, adata.shape[0] - 1)
        sc.pp.pca(adata, n_comps=knn)
        sc.pp.neighbors(adata, n_neighbors=knn, n_pcs=10)
        sc.tl.leiden(adata)

        # Run UMAP
        sc.tl.umap(adata)

        # Save the clusters to a TSV file
        mtx_barcodes_df = met_passing["cell_id"].to_frame()
        mtx_barcodes_df["cluster"] = adata.obs["leiden"].values
        mtx_barcodes_df.to_csv(writeDir / "csv" / f"{self.sample}.report_clusters.tsv", sep="\t", index=False)

        # Plot the UMAP embedding colored by cluster
        fig = go.Figure()
        unique_clusters = adata.obs["leiden"].unique().astype(int)
        colors = []
        # loop through colors if more than 24 clusters
        # Tooltips still exist to find the difference
        while len(unique_clusters) - len(colors) > 24:
            colors.extend(px.colors.qualitative.Dark24)
        colors.extend(px.colors.qualitative.Dark24[: len(unique_clusters)])
        if len(unique_clusters) == 1:
            color_map = {unique_clusters[0]: colors[0]}
        else:
            color_map = {cluster: colors[i] for i, cluster in enumerate(range(0, max(unique_clusters)))}
        for cluster, color in color_map.items():
            vals = adata.obs["leiden"].astype(int) == cluster
            fig.add_trace(
                go.Scatter(
                    x=adata.obsm["X_umap"][vals, 0],
                    y=adata.obsm["X_umap"][vals, 1],
                    mode="markers",
                    marker=dict(color=color),
                    legendgroup=str(cluster),
                    showlegend=True,
                    name=f"Cluster {cluster}",
                    hovertemplate=f"Cluster {cluster}",
                )
            )

        fig.update_layout(
            title="UMAP projection of the methylation data colored by cluster",
            xaxis_title="UMAP1",
            yaxis_title="UMAP2",
            width=800,
            height=800,
        )
        return dp.Plot(fig)

    def build_knee_plot(self) -> dp.Plot:
        figs = []
        for sample in self.all_cells["sample"].unique():
            df = self.all_cells[self.all_cells["sample"] == sample]
            sorted_uniq = np.sort(df["unique_reads"])[::-1]
            figs.append(
                go.Scatter(
                    x=np.arange(1, len(df) + 1),
                    y=sorted_uniq,
                    mode="lines",
                    name=sample,
                    line=dict(color=self.color_map[sample]),
                )
            )
        if not self.is_library_report:
            threshold = self.all_cells[self.all_cells["pass"] == "pass"]["unique_reads"].min()
            figs.append(
                go.Scatter(
                    x=[
                        min(np.arange(1, len(self.all_cells) + 1)),
                        max(np.arange(1, len(self.all_cells) + 1)),
                    ],
                    y=[threshold, threshold],
                    mode="lines",
                    line=dict(dash="dash", color="green"),
                )
            )

        layout = go.Layout(xaxis=dict(type="log"), yaxis=dict(type="log"))
        fig = go.Figure(data=figs, layout=layout)
        update_fig_layout(
            fig,
            "Barcode Rank Plot",
            "Barcode",
            "Unique Reads",
            False,
            800,
            500,
            self.is_library_report,
            "log",
            "log",
            None,
        )
        # Start plot from 10 -> log(10) = 1
        fig.update_layout(xaxis_range=[1, np.log10(len(df))])
        fig.write_image(f"{self.out_dir}/png/{self.sample}.kneePlot.png")
        return dp.Plot(fig)

    def build_fragment_length_histogram(self) -> dp.Plot:
        # Need to use histogram instead of bar plot to get the bars to be exactly next to each other, otherwise there are tiny white gaps
        fig = px.histogram(
            self.fragment_df,
            x="bin",
            y="count",
            histfunc="sum",
            color_discrete_sequence=[self.color_map[self.sample]],
            nbins=self.fragment_df["bin"].nunique(),
        )
        update_fig_layout(
            fig,
            "Fragment Length Histogram",
            "Fragment Length",
            "Fragment Counts",
            False,
            800,
            500,
            self.is_library_report,
            "linear",
            "linear",
            None,
        )
        fig.write_image(f"{self.out_dir}/png/{self.sample}.fragLenHist.png")
        return dp.Plot(fig)

    def create_total_complexity_plot(self):
        color_palette = px.colors.qualitative.Dark24
        passing_color_map = {
            sample + "(passing)": color_palette[i] for i, sample in enumerate(sorted(self.all_cells["sample"].unique()))
        }
        failing_color_map = {
            sample + "(background)": color_palette[i]
            for i, sample in enumerate(sorted(self.all_cells["sample"].unique()))
        }
        passing_df = self.all_cells[self.all_cells["pass"] == "pass"]
        passing_df["sample"] = passing_df["sample"] + "(passing)"
        passing_df["percent"] = round((passing_df["unique_reads"] / passing_df["total_reads"] * 100), 3)
        passing_fig = px.scatter(
            passing_df,
            x="percent",
            y="unique_reads",
            color="sample",
            opacity=0.8,
            color_discrete_map=passing_color_map,
        )

        failing_df = self.all_cells[self.all_cells["pass"] == "fail"]
        failing_df["sample"] = failing_df["sample"] + "(background)"
        failing_df["percent"] = round((failing_df["unique_reads"] / failing_df["total_reads"] * 100), 3)
        failing_fig = px.scatter(
            failing_df,
            x="percent",
            y="unique_reads",
            color="sample",
            opacity=0.05,
            color_discrete_map=failing_color_map,
        )

        if self.is_library_report:
            fig = go.Figure(data=passing_fig.data)
            # Add failing cells to the plot but turn off by default
            for trace in failing_fig.data:
                fig.add_trace(
                    go.Scatter(
                        x=trace.x,
                        y=trace.y,
                        mode="markers",
                        name=trace.name,
                        marker_color=trace.marker.color,
                        opacity=0.2,
                        visible="legendonly",
                        customdata=failing_df["sample"],
                        hovertemplate="sample=%{customdata}<br>"
                        + "percent=%{x}<br>"
                        + "unique_reads=%{y}<br>"
                        + "<extra></extra>",
                    )
                )
        else:
            fig = go.Figure(data=passing_fig.data + failing_fig.data)
        update_fig_layout(
            fig,
            "Reads Per Cell",
            "(Unique/Total Reads) %",
            "Unique Reads",
            False,
            800,
            500,
            True,
            "linear",
            "log",
            None,
        )
        # Start plot from 1 -> log(10)
        fig.update_layout(yaxis_range=[1, np.log10(self.all_cells["unique_reads"].max()) * 1.05])
        fig.write_image(f"{self.out_dir}/png/{self.sample}.complexityTotal.png")
        return dp.Plot(fig)

    def build_total_and_unique_reads_box(self) -> dp.Plot:
        tmp_df = self.all_cells[self.all_cells["pass"] == "pass"].rename(
            columns={"total_reads": "Total", "unique_reads": "Unique"}
        )
        IN = tmp_df[["sample", "Total", "Unique"]].melt(id_vars=["sample"], var_name="variable", value_name="value")
        fig = px.box(
            IN,
            x="variable",
            y="value",
            color="sample",
            points="suspectedoutliers",
            color_discrete_map=self.color_map,
        )
        update_fig_layout(
            fig,
            "Reads per Cell",
            "",
            "Reads",
            False,
            800,
            500,
            self.is_library_report,
            "category",
        )
        fig.write_image(f"{self.out_dir}/png/{self.sample}.cellTotalAndUniqueReadsBox.png")
        return dp.Plot(fig)

    def build_saturation_box(self) -> dp.Plot:
        tmp_df = self.all_cells[self.all_cells["pass"] == "pass"]
        tmp_df["saturation"] = round(100 - (tmp_df["unique_reads"] / tmp_df["total_reads"] * 100), 3)
        tmp_df = tmp_df.rename(columns={"saturation": ""})
        IN = tmp_df[["sample", ""]].melt(id_vars=["sample"], var_name="variable", value_name="value")
        fig = px.box(
            IN,
            x="variable",
            y="value",
            color="sample",
            points="suspectedoutliers",
            color_discrete_map=self.color_map,
        )
        update_fig_layout(
            fig,
            "Saturation",
            "",
            "Saturation %",
            False,
            800,
            500,
            self.is_library_report,
            "category",
        )
        fig.update_layout(yaxis=dict(range=[0, 100]))
        fig.write_image(f"{self.out_dir}/png/{self.sample}.cellSaturationBox.png")
        return dp.Plot(fig)

    def build_cell_table(self) -> dp.Table:
        cells_called = {}
        for sample in self.all_cells["sample"].unique():
            ncells = ((self.all_cells["sample"] == sample) & (self.all_cells["pass"] == "pass")).sum()
            cells_called[sample] = [f"{ncells:,}"]
        cells = pd.DataFrame(cells_called)
        cells = cells.T.reset_index()
        cells.columns = ["Sample", "Cells"]
        return DatapaneUtils.createTableIgnoreWarning(
            cells[["Sample", "Cells"]].style.pipe(DatapaneUtils.styleTable, title="Cells Called")
        )

    def build_bc_parser_stats(self, bcParserMetrics: Path | bool) -> dp.Table:
        if bcParserMetrics:
            with open(bcParserMetrics) as f:
                bc_metrics_dict = json.load(f)
            read_info = bc_metrics_dict["reads"]
            bc_metrics_df = pd.DataFrame(
                {
                    "Status": ["Pass", "BarcodeError", "LinkerError", "TooShortError"],
                    "Reads": [
                        read_info["Pass"][0],
                        read_info["BarcodeError"][0],
                        read_info["LinkerError"][0],
                        read_info["TooShortError"][0],
                    ],
                    "Percent": [
                        read_info["Pass"][1],
                        read_info["BarcodeError"][1],
                        read_info["LinkerError"][1],
                        read_info["TooShortError"][1],
                    ],
                }
            )
        else:
            bc_metrics_df = pd.DataFrame(
                {
                    "Status": ["Pass", "BarcodeError", "LinkerError", "TooShortError"],
                    "Reads": [np.nan, np.nan, np.nan, np.nan],
                    "Percent": [np.nan, np.nan, np.nan, np.nan],
                }
            )
        return DatapaneUtils.createTableIgnoreWarning(
            bc_metrics_df.style.pipe(
                DatapaneUtils.styleTable,
                title="Barcode Metrics",
                hideColumnHeaders=False,
                boldColumn=["Status"],
                numericCols=["Reads", "Percent"],
            )
        )

    def build_sample_information_table(self, sample_information: pd.DataFrame | bool) -> dp.Table:
        sample_info_df = pd.DataFrame()
        if not isinstance(sample_information, bool):
            sample_info_df = pd.DataFrame.from_dict(
                {
                    "Metric": [
                        "Sample Name",
                        "Library Name",
                        "Sample Barcodes",
                        "Workflow Version",
                        "Analysis Date",
                    ],
                    "Value": [
                        sample_information["sample_name"].to_list()[0],
                        sample_information["library_name"].to_list()[0],
                        sample_information["sample_barcodes"].to_list()[0],
                        sample_information["workflow_version"].to_list()[0],
                        sample_information["analysis_date"].to_list()[0],
                    ],
                }
            )
        else:
            sample_info_df = pd.DataFrame.from_dict(
                {
                    "Metric": [
                        "Sample Name",
                        "Library Name",
                        "Barcode Name",
                        "Workflow Version",
                        "Analysis Date",
                    ],
                    "Value": [np.nan, np.nan, np.nan, np.nan, np.nan],
                },
                dtype="object",
            )
        return DatapaneUtils.createTableIgnoreWarning(
            sample_info_df.style.pipe(
                DatapaneUtils.styleTable,
                title="Sample Information",
                hideColumnHeaders=True,
                boldColumn=["Metric"],
                numericCols=["Value"],
            )
        )

    def build_sample_information_table_combined(self, sample_information: pd.DataFrame | bool) -> dp.Table:
        sample_info_df = pd.DataFrame()
        if not isinstance(sample_information, bool):
            sample_info_df = pd.DataFrame.from_dict(
                {
                    "Metric": ["Library Name", "Workflow Version", "Analysis Date"],
                    "Value": [
                        sample_information["library_name"].to_list()[0],
                        sample_information["workflow_version"].to_list()[0],
                        sample_information["analysis_date"].to_list()[0],
                    ],
                }
            )
        else:
            sample_info_df = pd.DataFrame.from_dict(
                {
                    "Metric": ["Library Name", "Workflow Version", "Analysis Date"],
                    "Value": [np.nan, np.nan, np.nan],
                },
                dtype="object",
            )
        return DatapaneUtils.createTableIgnoreWarning(
            sample_info_df.style.pipe(
                DatapaneUtils.styleTable,
                title="Sample Information",
                hideColumnHeaders=True,
                boldColumn=["Metric"],
                numericCols=["Value"],
            )
        )

    def build_summary_stats_table(
        self, mapping_stats: pd.DataFrame | bool, trimming_stats: pd.DataFrame | bool
    ) -> dp.Table:
        metric_names = [
            "total_reads",
            "total_read_pairs",
            "percent_passing_trimming",
            "percent_passing_mapping",
            "reads_per_passing_cell",
            "saturation",
        ]
        if isinstance(mapping_stats, bool) and isinstance(trimming_stats, bool):
            summary_stats = pd.DataFrame.from_dict(
                {
                    "Metric": metric_names,
                    "Value": [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                },
                dtype="object",
            )
        elif self.met_df.empty:
            summary_stats = pd.DataFrame.from_dict(
                {
                    "Metric": metric_names,
                    "Value": [
                        int(trimming_stats[trimming_stats["Metric"] == "total_reads"]["Value"].to_list()[0]),
                        int(trimming_stats[trimming_stats["Metric"] == "total_reads"]["Value"].to_list()[0]) / 2,
                        trimming_stats[trimming_stats["Metric"] == "percent_passing"]["Value"].to_list()[0],
                        (
                            mapping_stats[mapping_stats["Metric"] == "mapped_reads"]["Value"].to_list()[0]
                            / trimming_stats[trimming_stats["Metric"] == "total_reads"]["Value"].to_list()[0]
                        ),
                        np.nan,
                        np.nan,
                    ],
                },
                dtype="object",
            )
        else:
            total_reads = f"{int(trimming_stats[trimming_stats['Metric'] == 'total_reads']['Value'].to_list()[0]):,}"
            total_read_pairs = (
                f"{int(trimming_stats[trimming_stats['Metric'] == 'total_reads']['Value'].to_list()[0]/2):,}"
            )
            trimmed_reads = f"{trimming_stats[trimming_stats['Metric'] == 'percent_passing']['Value'].to_list()[0] / 100:.1%}"  # Convert percentage to fraction
            mapped_reads = f"{(mapping_stats[mapping_stats['Metric'] == 'mapped_reads']['Value'].to_list()[0]/trimming_stats[trimming_stats['Metric'] == 'total_reads']['Value'].to_list()[0]):.1%}"
            summary_stats = pd.DataFrame.from_dict(
                {
                    "Metric": metric_names,
                    "Value": [
                        total_reads,
                        total_read_pairs,
                        trimmed_reads,
                        mapped_reads,
                        f"{int(trimming_stats[trimming_stats['Metric'] == 'total_reads']['Value'].to_list()[0]/len(self.met_df)):,}",
                        f"{round(1 - (self.all_cells['unique_reads'].sum()/self.all_cells['total_reads'].sum()), 3)}",
                    ],
                },
                dtype="object",
            )
        summary_stats.iloc[0, summary_stats.columns.get_loc("Metric")] = "Total Reads"
        summary_stats.iloc[1, summary_stats.columns.get_loc("Metric")] = "Total Read Pairs"
        summary_stats.iloc[2, summary_stats.columns.get_loc("Metric")] = "Reads Passing Trimming"
        summary_stats.iloc[3, summary_stats.columns.get_loc("Metric")] = "Reads Passing Mapping"
        summary_stats.iloc[4, summary_stats.columns.get_loc("Metric")] = "Reads Per Passing Cell"
        summary_stats.iloc[5, summary_stats.columns.get_loc("Metric")] = "Saturation"
        return DatapaneUtils.createTableIgnoreWarning(
            summary_stats[["Metric", "Value"]].style.pipe(
                DatapaneUtils.styleTable,
                title="Raw Read Stats",
                hideColumnHeaders=True,
                boldColumn=["Metric"],
                numericCols=["Value"],
            )
        )

    def construct_passing_cell_stats(self, passingCellStats) -> dp.Table:

        met_stats = pd.read_csv(passingCellStats, index_col=False)
        met_stats = met_stats[met_stats["Metric"] != "sample_name"]
        met_stats = met_stats.reset_index(drop=True)
        # Reformat metric names and numeric values for HTML report (vs. .csv)
        for idx, col_name in enumerate(CELL_STAT_COL_NAMES):
            met_stats.iloc[idx, met_stats.columns.get_loc("Metric")] = col_name
            if col_name in [
                "Number of Passing Cells",
                "Total Reads in Cells",
                "Total Reads in Passing Cells",
                "Median Reads per Cell",
                "Total Passing Reads",
                "Median Passing Reads per Cell",
                "Total Unique Reads",
                "Median Unique Reads per Cell",
                "Median Unique Reads",
                "Median CG's Covered per Cell",
                "Median CH's Covered per Cell",
            ]:
                met_stats.iloc[idx, 1] = int(met_stats.iloc[idx, 1])
                met_stats.iloc[idx, 1] = f"{met_stats.iloc[idx, 1]:,}"
            if col_name in [
                "Reads in Passing Cells",
                "Median CG Methylation",
                "Median CH Methylation",
            ]:
                met_stats.iloc[idx, 1] = f"{met_stats.iloc[idx, 1]:.1%}"
        if "tss_enrich" in self.met_df.columns:
            idx = len(CELL_STAT_COL_NAMES)
            new_row = pd.DataFrame(
                {"Metric": ["Median TSS Enrichment"], "Value": [f"{self.met_df['tss_enrich'].median():.2f}"]}
            )
            met_stats = pd.concat([met_stats, new_row], ignore_index=True)

        return DatapaneUtils.createTableIgnoreWarning(
            met_stats[["Metric", "Value"]].style.pipe(
                DatapaneUtils.styleTable,
                title="Cell Statistics",
                hideColumnHeaders=True,
                boldColumn=["Metric"],
                numericCols=["Value"],
            )
        )

    def build_cell_covered_box(self) -> dp.Plot:
        tmp_df = self.met_df.rename(columns={"cg_cov": "CG", "ch_cov": "CH"})
        IN = tmp_df[["sample", "CG", "CH"]].melt(id_vars=["sample"], var_name="variable", value_name="value")
        fig = px.box(
            IN,
            x="variable",
            y="value",
            color="sample",
            points="suspectedoutliers",
            color_discrete_map=self.color_map,
        )
        update_fig_layout(
            fig,
            "Coverage per Cell",
            "",
            "Cytosines Covered",
            False,
            800,
            500,
            self.is_library_report,
            "category",
            "log",
        )
        fig.write_image(f"{self.out_dir}/png/{self.sample}.cellCoveredBox.png")
        return dp.Plot(fig)

    def build_cg_cell_methyl_percent_box(self) -> dp.Plot:
        tmp_met_df = self.met_df.rename(columns={"mcg_pct": ""})
        IN = tmp_met_df[["sample", ""]].melt(id_vars=["sample"], var_name="", value_name="value")
        fig = px.box(
            IN,
            x="",
            y="value",
            color="sample",
            points="suspectedoutliers",
            color_discrete_map=self.color_map,
        )
        update_fig_layout(
            fig,
            "CG Methylation per Cell",
            "",
            "Methylation %",
            False,
            800,
            500,
            self.is_library_report,
            "category",
        )
        fig.update_layout(yaxis=dict(range=[0, 100]))
        fig.write_image(f"{self.out_dir}/png/{self.sample}.cellMethylCgPercentBox.png")
        return dp.Plot(fig)

    def build_ch_cell_methyl_percent_box(self) -> dp.Plot:
        tmp_met_df = self.met_df.rename(columns={"mch_pct": ""})
        IN = tmp_met_df[["sample", ""]].melt(id_vars=["sample"], var_name="", value_name="value")
        fig = px.box(
            IN,
            x="",
            y="value",
            color="sample",
            points="suspectedoutliers",
            color_discrete_map=self.color_map,
        )
        update_fig_layout(
            fig,
            "CH Methylation per Cell",
            "",
            "Methylation %",
            False,
            800,
            500,
            self.is_library_report,
            "category",
        )
        fig.update_layout(yaxis=dict(range=[0, IN["value"].max() + 1]))
        fig.write_image(f"{self.out_dir}/png/{self.sample}.cellMethylChPercentBox.png")
        return dp.Plot(fig)

    def build_tss_enrich_box(self) -> dp.Plot:
        tmp_met_df = self.met_df.rename(columns={"tss_enrich": ""})
        IN = tmp_met_df[["sample", ""]].melt(id_vars=["sample"], var_name="variable", value_name="value")
        fig = px.box(
            IN,
            x="variable",
            y="value",
            color="sample",
            points="suspectedoutliers",
            color_discrete_map=self.color_map,
        )
        update_fig_layout(
            fig,
            "TSS Enrichment per Cell",
            "",
            "TSS Fold Enrichment",
            False,
            800,
            500,
            self.is_library_report,
            "category",
        )
        fig.update_layout(yaxis=dict(range=[0, IN["value"].max() + 1]))
        fig.write_image(f"{self.out_dir}/png/{self.sample}.cellTssEnrichBox.png")
        return dp.Plot(fig)

    def build_mito_box(self) -> dp.Plot:
        tmp_met_df = self.all_cells.rename(columns={"mito_reads": ""})
        IN = tmp_met_df[
            [
                "sample",
                "",
            ]
        ].melt(id_vars=["sample"], var_name="variable", value_name="value")
        fig = px.box(
            IN,
            x="variable",
            y="value",
            color="sample",
            points="suspectedoutliers",
            color_discrete_map=self.color_map,
        )
        update_fig_layout(
            fig,
            "Mitochondrial Reads",
            "",
            "% Mito Reads",
            False,
            800,
            500,
            self.is_library_report,
            "category",
        )
        fig.update_layout(yaxis=dict(range=[0, IN["value"].max() + 1]))
        fig.write_image(f"{self.out_dir}/png/{self.sample}.mitoReadsBox.png")
        return dp.Plot(fig)

    def build_mito_tss_ch_table(self) -> dp.Table:
        mito_tss_dist_df = pd.DataFrame(
            {
                "Cells with % of Mito Reads > 1": [len(self.all_cells[self.all_cells["mito_reads"] > 1])],
                "Cells with TSS Enrichment > 3": [
                    (
                        len(self.all_cells[self.all_cells["tss_enrich"] > 3])
                        if "tss_enrich" in self.all_cells.columns
                        else pd.NA
                    )
                ],
                "Cells with CH Methylation % > 1": [len(self.met_df[self.met_df["mch_pct"] > 1])],
            }
        )
        return DatapaneUtils.createTableIgnoreWarning(
            mito_tss_dist_df.style.pipe(
                DatapaneUtils.styleTable,
                title="Mitochondrial Reads and TSS Enrichment Distribution",
                hideColumnHeaders=False,
                boldColumn=None,
            )
        )


def update_fig_layout(
    fig,
    title_text="",
    xaxis_title="",
    yaxis_title="",
    autosize=False,
    width=800,
    height=500,
    showlegend=True,
    xaxis_type="linear",
    yaxis_type="linear",
    boxgroupgap=0.5,
):
    fig.update_layout(
        title={"text": title_text, "x": 0.5, "xanchor": "center"},
        xaxis_title=xaxis_title,
        xaxis_type=xaxis_type,
        yaxis_type=yaxis_type,
        yaxis_title=yaxis_title,
        autosize=autosize,
        width=width,
        height=height,
        boxgroupgap=boxgroupgap,
        showlegend=showlegend,
    )
