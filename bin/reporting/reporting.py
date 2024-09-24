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

# Suppress SettingWithCopyWarning
pd.options.mode.chained_assignment = None

CELL_STAT_COL_NAMES = ["Number of Passing Cells",
                       "Reads in Passing Cells",
                       "Median Unique Reads",
                       "Median CG's Covered",
                       "Median CH's Covered",
                       "Median CG Methylation",
                       "Median CH Methylation",
                       ]
class Utils:
    def get_passing_cells(all_cells:Path, sample:str) -> list[str]:
        """
        Get a list of cell barcodes that passed filters

        Args:
            all_cells: Path to the file containing all cell metrics
            sample: Used for getting the well coordinate
        """
        passing_cells = []
        with open(all_cells) as csv_file:
            csv_reader = csv.DictReader(csv_file)
            for row in csv_reader:
                # splitFastq = true
                if "." in sample:
                    if row['pass'] == 'pass' and row['tgmt_well'] == sample.split('.')[1]:
                        passing_cells.append(row['cell_id'])
                # splitFastq = false
                else:
                    if row['pass'] == 'pass' and row['sample'] == sample:
                        passing_cells.append(row['cell_id'])
        
        if len(passing_cells) == 0:
            print(f"No passing cells found for {sample}")
            sys.exit(0)
        return passing_cells

class WhitelistAndLibJsonHelper:
    @staticmethod
    def get_barcode_entry_from_libjson(library_structure_json: Path, name: str) -> dict:
        with open(library_structure_json) as f:
            lib_json = json.load(f)
        for entry in lib_json['barcodes']:
            if entry['name'] == name:
                return entry

    @staticmethod
    def get_sample_barcode_name(library_structure_json: Path) -> str:
        with open(library_structure_json) as f:
            lib_json = json.load(f)
        return lib_json['sample_barcode']
    
    @staticmethod
    def get_barcode_whitelist_dataframe_from_file(barcodes: Path) -> pd.DataFrame:
        return pd.read_csv(barcodes, sep="\t", names=["barcode", "well"])

    @staticmethod
    def get_barcode_whitelist_dataframe_from_libjson(library_structure_json: Path, name: str) -> pd.DataFrame:
        entry = WhitelistAndLibJsonHelper.get_barcode_entry_from_libjson(library_structure_json, name)
        return pd.read_csv(entry['sequences'], sep="\t", names=["barcode", "well"])

    @staticmethod
    def get_all_barcode_names(library_structure_json: Path) -> list:
        with open(library_structure_json) as f:
            lib_json = json.load(f)
        return [entry['name'] for entry in lib_json['barcodes']]

    @staticmethod
    def get_max_number_and_letter(barcodes: Path):
        barcodes_df = pd.read_csv(barcodes, sep="\t", names=["barcode", "well"])
        # Last line corresponds to maximum well coordinate
        # First character denotes plate number so removing that
        max_well_coordinate = barcodes_df.iloc[-1, 1][1:]
        max_letter = max_well_coordinate[0]
        max_number = int(max_well_coordinate[1:])
        return max_letter, max_number

    @staticmethod
    def get_number_of_plates(barcodes: Path):
        barcodes_df = pd.read_csv(barcodes, sep="\t", names=["barcode", "well"])
        # Last line corresponds to maximum well coordinate
        # First character denotes plate number
        number_of_plates = barcodes_df.iloc[-1]['well'][0]
        return int(number_of_plates)
    
    @staticmethod
    def get_lib_json(library_structure_json: Path):
        with open(library_structure_json) as f:
            lib_json = json.load(f)
        return lib_json
    
class PlatePlotBuilder:
    def build_combined_plate_plot(passing_cells: pd.DataFrame, identifier: str, lib_name: str, write_dir: Path, lib_json: Path) -> list:
        num_cells_dict = passing_cells[f"{identifier}_well"].value_counts().to_dict()
        datapane_list = []
        sample_barcode = WhitelistAndLibJsonHelper.get_sample_barcode_name(lib_json)
        barcode_entry = WhitelistAndLibJsonHelper.get_barcode_entry_from_libjson(lib_json, sample_barcode)
        sequence_path = lib_json.parent / barcode_entry['sequences']
        max_letter, max_number = WhitelistAndLibJsonHelper.get_max_number_and_letter(sequence_path)
        well_plate_total_reads = pd.DataFrame(0, columns=range(1, max_number+1), index=[chr(j) for j in range(65, ord(max_letter)+1)])
        well_plate_total_cells = well_plate_total_reads.copy()
        df_by_plate = passing_cells[["total_reads", f"{identifier}_well"]]
        num_cells_dict = df_by_plate[f"{identifier}_well"].value_counts().to_dict()
        df_by_plate = df_by_plate.groupby(f"{identifier}_well", as_index=False)["total_reads"].sum()

        for _, row in df_by_plate.iterrows():
            letter = row[f'{identifier}_well'][0]
            number = int(row[f'{identifier}_well'][1:])
            well_plate_total_reads.at[letter, number] = row["total_reads"]
        
        for well in num_cells_dict:
            letter = well[0]
            number = int(well[1:])
            well_plate_total_cells.at[letter, number] = num_cells_dict[well]
        well_plate_total_cells.to_csv(write_dir / "csv" / f"{lib_name}.total_passing_cells_{identifier}.csv")
        well_plate_total_reads.to_csv(write_dir / "csv" / f"{lib_name}.total_reads_{identifier}.csv")
        datapane_list.append(PlatePlotBuilder.plot_plate_plot(well_plate_total_reads, f"Total Reads({identifier})", 100))
        datapane_list.append(PlatePlotBuilder.plot_plate_plot(well_plate_total_cells, f"Total Passing Cells({identifier})", False))

        return datapane_list

    def build_plate_plot(all_cells: pd.DataFrame, lib_json: Path, write_dir: Path, sample: str) -> list:
        datapane_list = []
        sample_barcode = WhitelistAndLibJsonHelper.get_sample_barcode_name(lib_json)
        barcode_entry = WhitelistAndLibJsonHelper.get_barcode_entry_from_libjson(lib_json, sample_barcode)
        sequence_path = lib_json.parent / barcode_entry['sequences']
        max_letter, max_number = WhitelistAndLibJsonHelper.get_max_number_and_letter(sequence_path)
        for i in range(1, WhitelistAndLibJsonHelper.get_number_of_plates(sequence_path)+1):
            well_plate_total_reads = pd.DataFrame(0, columns=range(1, max_number+1), index=[chr(j) for j in range(65, ord(max_letter)+1)])
            well_plate_total_cells = well_plate_total_reads.copy()
            df_by_plate = all_cells[all_cells[f"{sample_barcode}_well"].str[0] == str(i)]
            df_by_plate = df_by_plate[["total_reads", f"{sample_barcode}_well"]]
            df_by_plate[f'{sample_barcode}_well'] = df_by_plate[f'{sample_barcode}_well'].str[1:]
            num_cells_dict = df_by_plate[f"{sample_barcode}_well"].value_counts().to_dict()
            df_by_plate = df_by_plate.groupby(f"{sample_barcode}_well", as_index=False)["total_reads"].sum()
            
            for _, row in df_by_plate.iterrows():
                letter = row[f'{sample_barcode}_well'][0]
                number = int(row[f'{sample_barcode}_well'][1:])
                well_plate_total_reads.at[letter, number] = row["total_reads"]
            
            for well in num_cells_dict:
                letter = well[0]
                number = int(well[1:])
                well_plate_total_cells.at[letter, number] = num_cells_dict[well]
            well_plate_total_cells.to_csv(write_dir / "csv" / f"{sample}.total_passing_cells_{sample_barcode}_plate_{i}.csv")
            well_plate_total_reads.to_csv(write_dir / "csv" / f"{sample}.total_reads_{sample_barcode}_plate_{i}.csv")
            datapane_list.append(PlatePlotBuilder.plot_plate_plot(well_plate_total_reads, f"Total Reads({sample_barcode}): Plate {i}", 100))
            datapane_list.append(PlatePlotBuilder.plot_plate_plot(well_plate_total_cells, f"Total Passing Cells({sample_barcode}): Plate {i}", False))
        
        return datapane_list
    
    def plot_plate_plot(df: pd.DataFrame, title: str, threshold: int | bool):
        fig = plt.figure()
        if threshold:
            norm = colors.SymLogNorm(linthresh=threshold, vmin=0, vmax=max(threshold*10, df.values.max()))
            ax = sns.heatmap(df, linewidth=0.5, cmap=sns.color_palette("dark:#80AAFF", as_cmap=True), norm=norm)
        else:
            if (df == 0).all().all():
                norm = colors.SymLogNorm(linthresh=1, vmin=0, vmax=max(1*10, df.values.max()))
                ax = sns.heatmap(df, linewidth=0.5, cmap=sns.color_palette("dark:#80AAFF", as_cmap=True), norm=norm)
            else:
                ax = sns.heatmap(df, linewidth=0.5, cmap=sns.color_palette("dark:#80AAFF", as_cmap=True), vmin=0.0, vmax=df.values.max())
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
        if (numericCols is not None):
            styler.format(DatapaneUtils.reformatIfInt, subset=numericCols)

        # Format NA values
        styler.format(na_rep='NA')

        styler.hide(axis='index')

        if hideColumnHeaders:
            styler.hide(axis='columns')
        else:
            styler.set_table_styles([{'selector': 'th', 'props': [('border', 'black solid !important')]}], overwrite=False)

        if boldColumn is not None:
            styler.set_properties(subset=boldColumn, **{'font-weight': 'bold'})

        if title != "":
            styler.set_caption(title)
            styler.set_table_styles(
                [{'selector': 'caption', 'props': [('border', 'black solid !important'), ("font-weight", "bold")]}], overwrite=False)

        styler.set_properties(**{"border-color": 'black', "border-style": 'solid !important'})
        return styler
    
    @staticmethod
    def define_color_map(df: pd.DataFrame) -> dict:
        color_palette = px.colors.qualitative.Light24
        return {sample: color_palette[i] for i, sample in enumerate(sorted(df['sample'].unique()))}

class BuildDatapane:
    """
    Contains functions for building the datapane html report
    """
    def __init__(self, sample: str, write_dir: Path, all_cells: pd.DataFrame,
                 is_library_report: bool, fragment_df: pd.DataFrame = pd.DataFrame()):
        self.all_cells = all_cells
        # sample = libName if this is a combined sample report
        self.sample = sample
        self.out_dir = write_dir
        self.fragment_df = fragment_df
        self.is_library_report = is_library_report
        self.met_df = all_cells[all_cells["pass"]=="pass"].reset_index(drop=True)
        if len(all_cells['sample'].unique()) > 1:
            self.color_map = DatapaneUtils.define_color_map(all_cells)
        else:
            self.color_map = {all_cells.iloc[0]['sample']: '#636EFA'}

    def build_knee_plot(self) -> dp.Plot:
        figs = []
        for sample in self.all_cells["sample"].unique():
            df = self.all_cells[self.all_cells["sample"]==sample]
            sorted_uniq = np.sort(df['unique_reads'])[::-1]
            figs.append(go.Scatter(x=np.arange(1, len(df)+1), y=sorted_uniq, mode='lines', name=sample, line=dict(color=self.color_map[sample])))
        if not self.is_library_report:
            threshold = self.all_cells[self.all_cells['pass'] == 'pass']['unique_reads'].min()
            figs.append(go.Scatter(x=[min(np.arange(1, len(self.all_cells)+1)),
                                   max(np.arange(1, len(self.all_cells)+1))],
                                   y=[threshold, threshold],
                                   mode='lines', line=dict(dash="dash", color="green")))

        layout = go.Layout(xaxis=dict(type='log'), yaxis=dict(type='log'))
        fig = go.Figure(data=figs, layout=layout)
        update_fig_layout(fig, 'Barcode Rank Plot', 'Barcode', 'Unique Reads', False, 800, 500, self.is_library_report, 'log', 'log', None)
        # Start plot from 10 -> log(10) = 1
        fig.update_layout(xaxis_range=[1, np.log10(len(df))])
        fig.write_image(f"{self.out_dir}/png/{self.sample}.kneePlot.png")
        return dp.Plot(fig)

    def build_fragment_length_histogram(self) -> dp.Plot:
        fig = px.bar(self.fragment_df, x='bin', y='count')
        update_fig_layout(fig, 'Fragment Length Histogram', 'Fragment Length', 'Fragment Counts', False, 800, 500, self.is_library_report, 'linear', 'linear', None)
        fig.write_image(f"{self.out_dir}/png/{self.sample}.fragLenHist.png")
        return dp.Plot(fig)
    
    def create_total_complexity_plot(self):
        if len(self.all_cells['sample'].unique()) > 1:
            color_palette = px.colors.qualitative.Light24
            passing_color_map = {sample+"(passing)": color_palette[i] for i, sample in enumerate(sorted(self.all_cells['sample'].unique()))}
            failing_color_map = {sample+"(background)": color_palette[i] for i, sample in enumerate(sorted(self.all_cells['sample'].unique()))}
        else:
            passing_color_map = {self.all_cells.iloc[0]['sample']+"(passing)": '#636EFA'} # blue
            failing_color_map = {self.all_cells.iloc[0]['sample']+"(background)": '#808080'} # gray
        passing_df = self.all_cells[self.all_cells["pass"] == "pass"]
        passing_df['sample'] = passing_df['sample'] + "(passing)"
        passing_df['percent'] = round((passing_df['unique_reads'] / passing_df['total_reads'] * 100), 3)
        passing_fig = px.scatter(passing_df, x='percent', y='unique_reads', color='sample', opacity=0.8, color_discrete_map=passing_color_map)

        failing_df = self.all_cells[self.all_cells["pass"]=="fail"]
        failing_df['sample'] = failing_df['sample'] + "(background)"
        failing_df['percent'] = round((failing_df['unique_reads'] / failing_df['total_reads'] * 100), 3)
        failing_fig = px.scatter(failing_df, x='percent', y='unique_reads', color='sample', opacity=0.2, color_discrete_map=failing_color_map)

        if self.is_library_report:
            fig = go.Figure(data=passing_fig.data)
            # Add failing cells to the plot but turn off by default
            for trace in failing_fig.data:
                fig.add_trace(
                    go.Scatter(
                        x=trace.x,
                        y=trace.y,
                        mode='markers',
                        name=trace.name,
                        marker_color=trace.marker.color,
                        opacity=0.2,
                        visible='legendonly',
                        customdata=failing_df['sample'],
                        hovertemplate=
                        'sample=%{customdata}<br>' +
                        'percent=%{x}<br>' +
                        'unique_reads=%{y}<br>' +
                        '<extra></extra>'
                    )
                )
        else:
            fig = go.Figure(data=passing_fig.data + failing_fig.data)
        update_fig_layout(fig, 'Reads Per Cell', '(Unique/Total Reads) %', 'Unique Reads', False, 800, 500, True, 'linear', 'log', None)
        # Start plot from 1 -> log(10)
        fig.update_layout(yaxis_range=[1, np.log10(self.all_cells['unique_reads'].max())*1.05])
        fig.write_image(f"{self.out_dir}/png/{self.sample}.complexityTotal.png")
        return dp.Plot(fig)
    
    def build_total_and_unique_reads_box(self) -> dp.Plot:
        tmp_df = self.all_cells[self.all_cells['pass'] == 'pass'].rename(columns={"total_reads": "Total", "unique_reads": "Unique"})
        IN = tmp_df[["sample", "Total", "Unique"]].melt(id_vars=["sample"], var_name="variable", value_name="value")
        fig = px.box(IN, x='variable', y='value', color='sample', points='suspectedoutliers', color_discrete_map=self.color_map)
        update_fig_layout(fig, 'Reads per Cell', '', 'Reads', False, 800, 500, self.is_library_report, 'category')
        fig.write_image(f"{self.out_dir}/png/{self.sample}.cellTotalAndUniqueReadsBox.png")
        return dp.Plot(fig)

    def build_saturation_box(self) -> dp.Plot:
        tmp_df = self.all_cells[self.all_cells['pass'] == 'pass']
        tmp_df['saturation'] = round(100 - (tmp_df['unique_reads'] / tmp_df['total_reads'] * 100), 3)
        tmp_df = tmp_df.rename(columns={"saturation": ""})
        IN = tmp_df[["sample", ""]].melt(id_vars=["sample"], var_name="variable", value_name="value")
        fig = px.box(IN, x='variable', y='value', color='sample', points='suspectedoutliers', color_discrete_map=self.color_map)
        update_fig_layout(fig, 'Saturation', '', 'Saturation %', False, 800, 500, self.is_library_report, 'category')
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
        cells.columns = ['Sample', 'Cells']
        return DatapaneUtils.createTableIgnoreWarning(cells[["Sample", "Cells"]].style.pipe(
            DatapaneUtils.styleTable, title="Cells Called"))

    def build_combined_passing_met_stats_table(self, combined_passing_met: pd.DataFrame) -> dp.Table:
        for idx, col_name in enumerate(CELL_STAT_COL_NAMES):
            combined_passing_met.iloc[idx, combined_passing_met.columns.get_loc('Metric')] = col_name
        if 'tss_enrich' in combined_passing_met.columns:
            combined_passing_met.iloc[len(CELL_STAT_COL_NAMES), combined_passing_met.columns.get_loc('Metric')] = "Median TSS Enrichment"
        combined_passing_met.to_csv(f"{self.out_dir}/{self.sample}.combinedPassingCellStats.csv", index=False)

    def build_bc_parser_stats(self, bcParserMetrics: Path | bool) -> dp.Table:
        if bcParserMetrics:
            with open(bcParserMetrics) as f:
                bc_metrics_dict = json.load(f)
            read_info = bc_metrics_dict['reads']
            bc_metrics_df = pd.DataFrame({
                'Status': ['Pass', 'BarcodeError', 'LinkerError', 'TooShortError'],
                'Reads': [read_info['Pass'][0], read_info['BarcodeError'][0], read_info['LinkerError'][0], read_info['TooShortError'][0]],
                'Percent': [read_info['Pass'][1], read_info['BarcodeError'][1], read_info['LinkerError'][1], read_info['TooShortError'][1]]
            })
        else:
            bc_metrics_df = pd.DataFrame({
                'Status': ['Pass', 'BarcodeError', 'LinkerError', 'TooShortError'],
                'Reads': [np.nan, np.nan, np.nan, np.nan],
                'Percent': [np.nan, np.nan, np.nan, np.nan]
            })
        return (DatapaneUtils.createTableIgnoreWarning(bc_metrics_df.style.pipe(
            DatapaneUtils.styleTable, title="Barcode Metrics", hideColumnHeaders=False, boldColumn=['Status'], numericCols=['Reads', 'Percent'])))

    def build_summary_stats(self, mapping_stats: pd.DataFrame | bool, trimming_stats: pd.DataFrame | bool):
        if isinstance(mapping_stats, bool) and isinstance(trimming_stats, bool):
            summary_stats = pd.DataFrame.from_dict({
                "Metric": ["total_reads", "reads_passing_trimming", "reads_passing_mapping"],
                "Value": [np.nan, np.nan, np.nan]
            }, dtype="object")
        else:
            summary_stats = pd.DataFrame.from_dict({
                "Metric": ["total_reads", "reads_passing_trimming", "reads_passing_mapping"],
                "Value": [int(trimming_stats[trimming_stats["Metric"] == "total_reads"]["Value"].to_list()[0]), trimming_stats[trimming_stats["Metric"] == "passing_reads"]["Value"].to_list()[0],
                          (mapping_stats[mapping_stats["Metric"] == "mapped_reads"]["Value"].to_list()[0])]
            }, dtype="object")
        summary_stats.to_csv(f"{self.out_dir}/csv/{self.sample}.summaryStats.csv", index=False)

    def build_summary_stats_table(self, mapping_stats: pd.DataFrame | bool, trimming_stats: pd.DataFrame | bool) -> dp.Table:
        metric_names = ["total_reads", "percent_passing_trimming", "percent_passing_mapping", "reads_per_passing_cell", "saturation"]
        if isinstance(mapping_stats, bool) and isinstance(trimming_stats, bool):
            summary_stats = pd.DataFrame.from_dict({
                "Metric": metric_names,
                "Value": [np.nan, np.nan, np.nan, np.nan, np.nan]
            }, dtype="object")
        else:
            total_reads = f"{int(trimming_stats[trimming_stats['Metric'] == 'total_reads']['Value'].to_list()[0]):,}"
            trimmed_reads = f"{trimming_stats[trimming_stats['Metric'] == 'percent_passing']['Value'].to_list()[0] / 100:.1%}" # Convert percentage to fraction
            mapped_reads = f"{(mapping_stats[mapping_stats['Metric'] == 'mapped_reads']['Value'].to_list()[0]/trimming_stats[trimming_stats['Metric'] == 'total_reads']['Value'].to_list()[0]):.1%}"
            if self.met_df.empty:
                summary_stats = pd.DataFrame.from_dict({
                    "Metric": metric_names,
                    "Value": [total_reads, trimmed_reads, mapped_reads, np.nan, np.nan]
                }, dtype="object")
            else:
                summary_stats = pd.DataFrame.from_dict({
                    "Metric": metric_names,
                    "Value": [total_reads, trimmed_reads, mapped_reads,
                            f"{int(trimming_stats[trimming_stats['Metric'] == 'total_reads']['Value'].to_list()[0]/len(self.met_df)):,}",
                            f"{1 - (self.all_cells['unique_reads'].sum()/self.all_cells['total_reads'].sum()):.2f}"
                            ]
                }, dtype="object")
        summary_stats.iloc[0, summary_stats.columns.get_loc('Metric')] = "Total Reads"
        summary_stats.iloc[1, summary_stats.columns.get_loc('Metric')] = "Reads Passing Trimming"
        summary_stats.iloc[2, summary_stats.columns.get_loc('Metric')] = "Reads Passing Mapping"
        summary_stats.iloc[3, summary_stats.columns.get_loc('Metric')] = "Reads Per Passing Cell"
        summary_stats.iloc[4, summary_stats.columns.get_loc('Metric')] = "Saturation"
        return (DatapaneUtils.createTableIgnoreWarning(summary_stats[["Metric", "Value"]].style.pipe(
            DatapaneUtils.styleTable, title="Read Stats", hideColumnHeaders=True, boldColumn=['Metric'], numericCols=['Value'])))

    def construct_passing_cell_stats(self) -> dp.Table:
        metric_name = ['number_passing_cells',
                       'percent_reads_in_passing_cells',
                       'median_uniq_reads',
                       'CGcov_median',
                       'CHcov_median',
                       'CG_mC_Pct_median',
                       'CH_mC_Pct_median',
                       'sample_name'
                      ]
        metric_values = [len(self.met_df),
                         self.all_cells[self.all_cells['pass'] == 'pass']['total_reads'].sum()/self.all_cells['total_reads'].sum(),
                         int(self.met_df['unique_reads'].median()),
                         int(self.met_df['cg_cov'].median()),
                         int(self.met_df['ch_cov'].median()),
                         self.met_df['mcg_pct'].median()/100, # Convert percentage to fraction
                         self.met_df['mch_pct'].median()/100,
                         self.sample
                        ]
        if 'tss_enrich' in self.met_df.columns:
            metric_name.append('median_tss_enrich')
            metric_values.append(self.met_df['tss_enrich'].median())
        
        if self.met_df.empty:
            met_stats = pd.DataFrame.from_dict({
                'Metric': metric_name,
                'Value': [np.nan]* len(metric_name)
            }, dtype="object")
        else:
            met_stats = pd.DataFrame.from_dict({
                "Metric": metric_name,
                "Value": metric_values
            }, dtype="object")
        met_stats.to_csv(f"{self.out_dir}/csv/{self.sample}.passingCellSummaryStats.csv", index=False)
        met_stats = met_stats[met_stats['Metric'] != 'sample_name']
        met_stats = met_stats.reset_index(drop=True)

        # Reformat metric names and numeric values for HTML report (vs. .csv)
        for idx, col_name in enumerate(CELL_STAT_COL_NAMES):
            met_stats.iloc[idx, met_stats.columns.get_loc('Metric')] = col_name
            if col_name in ["Number of Passing Cells", "Median Unique Reads", "Median CG's Covered", "Median CH's Covered"]:
                met_stats.iloc[idx, 1] = f"{met_stats.iloc[idx, 1]:,}"
            if col_name in ["Reads in Passing Cells", "Median CG Methylation", "Median CH Methylation"]:
                met_stats.iloc[idx, 1] = f"{met_stats.iloc[idx, 1]:.1%}"
        if 'tss_enrich' in self.met_df.columns:
            idx = len(CELL_STAT_COL_NAMES)
            met_stats.iloc[idx, met_stats.columns.get_loc('Metric')] = "Median TSS Enrichment"
            met_stats.iloc[idx, 1] = f"{met_stats.iloc[idx, 1]:.2f}"

        return (DatapaneUtils.createTableIgnoreWarning(met_stats[["Metric", "Value"]].style.pipe(
            DatapaneUtils.styleTable, title="Cell Statistics", hideColumnHeaders=True, boldColumn=['Metric'], numericCols=['Value'])))

    def build_cell_covered_box(self) -> dp.Plot:
        tmp_df = self.met_df.rename(columns={"cg_cov": "CG", "ch_cov": "CH"})
        IN = tmp_df[['sample', 'CG', 'CH']].melt(id_vars=['sample'], var_name='variable', value_name='value')
        fig = px.box(IN, x='variable', y='value', color='sample', points='suspectedoutliers', color_discrete_map=self.color_map)
        update_fig_layout(fig, 'Coverage per Cell', '', 'Cytosines Covered', False, 800, 500, self.is_library_report, 'category', 'log')
        fig.write_image(f"{self.out_dir}/png/{self.sample}.cellCoveredBox.png")
        return dp.Plot(fig)

    def build_cg_cell_methyl_percent_box(self) -> dp.Plot:
        tmp_met_df = self.met_df.rename(columns={"mcg_pct": ""})
        IN = tmp_met_df[["sample", ""]].melt(id_vars=["sample"], var_name="", value_name="value")
        fig = px.box(IN, x='', y='value', color='sample', points='suspectedoutliers', color_discrete_map=self.color_map)
        update_fig_layout(fig, 'CG Methylation per Cell', '', 'Methylation %', False, 800, 500, self.is_library_report, 'category')
        fig.update_layout(yaxis=dict(range=[0, 100]))
        fig.write_image(f"{self.out_dir}/png/{self.sample}.cellMethylCgPercentBox.png")
        return dp.Plot(fig)

    def build_ch_cell_methyl_percent_box(self) -> dp.Plot:
        tmp_met_df = self.met_df.rename(columns={"mch_pct": ""})
        IN = tmp_met_df[["sample", ""]].melt(id_vars=["sample"], var_name="", value_name="value")
        fig = px.box(IN, x='', y='value', color='sample', points='suspectedoutliers', color_discrete_map=self.color_map)
        update_fig_layout(fig, 'CH Methylation per Cell', '', 'Methylation %', False, 800, 500, self.is_library_report, 'category')
        fig.update_layout(yaxis=dict(range=[0, IN["value"].max() + 1]))
        fig.write_image(f"{self.out_dir}/png/{self.sample}.cellMethylChPercentBox.png")
        return dp.Plot(fig)
    
    def build_tss_enrich_box(self) -> dp.Plot:
        tmp_met_df = self.met_df.rename(columns={"tss_enrich": ""})
        IN = tmp_met_df[['sample', '']].melt(id_vars=["sample"], var_name="variable", value_name="value")
        fig = px.box(IN, x='variable', y='value', color='sample', points='suspectedoutliers', color_discrete_map=self.color_map)
        update_fig_layout(fig, 'TSS Enrichment per Cell', '', 'TSS Fold Enrichment', False, 800, 500, self.is_library_report, 'category')
        fig.update_layout(yaxis=dict(range=[0, IN["value"].max() + 1]))
        fig.write_image(f"{self.out_dir}/png/{self.sample}.cellTssEnrichBox.png")
        return dp.Plot(fig)

    def build_mito_box(self) -> dp.Plot:
        tmp_met_df = self.all_cells.rename(columns={"mito_reads": ""})
        IN = tmp_met_df[['sample', '',]].melt(id_vars=["sample"], var_name="variable", value_name="value")
        fig = px.box(IN, x='variable', y='value', color='sample', points='suspectedoutliers')
        update_fig_layout(fig, 'Mitochondrial Reads', '', '% Mito Reads', False, 800, 500, self.is_library_report, 'category')
        fig.update_layout(yaxis=dict(range=[0, IN["value"].max() + 1]))
        fig.write_image(f"{self.out_dir}/png/{self.sample}.mitoReadsBox.png")
        return dp.Plot(fig)

    def build_mito_tss_ch_table(self) -> dp.Table:
        mito_tss_dist_df = pd.DataFrame({
            'Cells with % of Mito Reads > 1': [len(self.all_cells[self.all_cells["mito_reads"] > 1])],
            'Cells with TSS Enrichment > 3': [len(self.all_cells[self.all_cells["tss_enrich"] > 3]) if 'tss_enrich' in self.all_cells.columns else pd.NA],
            'Cells with CH Methylation % > 1': [len(self.met_df[self.met_df['mch_pct'] > 1])]
        })
        return (DatapaneUtils.createTableIgnoreWarning(mito_tss_dist_df.style.pipe(
            DatapaneUtils.styleTable, title="Mitochondrial Reads and TSS Enrichment Distribution", hideColumnHeaders=False, boldColumn=None)))

def update_fig_layout(fig,
                      title_text='',
                      xaxis_title='',
                      yaxis_title='',
                      autosize=False,
                      width=800,
                      height=500,
                      showlegend=True,
                      xaxis_type='linear',
                      yaxis_type='linear',
                      boxgroupgap=0.5):
    fig.update_layout(
        title={'text': title_text, 'x': 0.5, 'xanchor': 'center'},
        xaxis_title=xaxis_title,
        xaxis_type=xaxis_type,
        yaxis_type=yaxis_type,
        yaxis_title=yaxis_title,
        autosize=autosize,
        width=width,
        height=height,
        boxgroupgap=boxgroupgap,
        showlegend=showlegend
    )
