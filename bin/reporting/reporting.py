import pandas as pd
import plotly.express as px
import datapane as dp
import warnings
import re
import json
import glob
import os
import numpy as np
import plotly.graph_objects as go
import matplotlib.colors as colors
import seaborn as sns
from pathlib import Path
from matplotlib import pyplot as plt

def build_combined_plate_plot(cell_stats_complexity_df: pd.DataFrame, barcodes: Path, identifier: str, tgmt_barcodes: Path, libraryName: str, writeDir: Path) -> list:
    num_cells_dict = cell_stats_complexity_df[f"{identifier}_well"].value_counts().to_dict()
    datapane_list = []
    max_letter, max_number = get_max_number_and_letter(tgmt_barcodes)
    wellPlateTotalReadsDf = pd.DataFrame(0, columns=range(1, max_number+1), index=[chr(j) for j in range(65, ord(max_letter)+1)])
    wellPlateTotalCellsDf = wellPlateTotalReadsDf.copy()
    df_by_plate = cell_stats_complexity_df[["total", f"{identifier}_well"]]
    num_cells_dict = df_by_plate[f"{identifier}_well"].value_counts().to_dict()
    df_by_plate = df_by_plate.groupby(f"{identifier}_well", as_index=False)["total"].sum()
    for _, row in df_by_plate.iterrows():
        letter = row[f'{identifier}_well'][0]
        number = int(row[f'{identifier}_well'][1:])
        wellPlateTotalReadsDf.at[letter, number] = row["total"]
    
    for well in num_cells_dict:
        letter = well[0]
        number = int(well[1:])
        wellPlateTotalCellsDf.at[letter, number] = num_cells_dict[well]
    wellPlateTotalCellsDf.to_csv(writeDir / "csv" / f"{libraryName}.total_passing_cells_{identifier}.csv")
    wellPlateTotalReadsDf.to_csv(writeDir / "csv" / f"{libraryName}.total_reads_{identifier}.csv")
    datapane_list.append(plot_plate_plot(wellPlateTotalReadsDf, f"Total Reads({identifier})", 100))
    datapane_list.append(plot_plate_plot(wellPlateTotalCellsDf, f"Total Passing Cells({identifier})", False))

    return datapane_list

def get_max_number_and_letter(tgmt_barcodes: Path):
    barcodes_df = pd.read_csv(tgmt_barcodes, sep="\t", names=["barcode", "well"])
    well_list = barcodes_df["well"].to_list()
    max_number = 1
    max_letter = 'A'
    for well in well_list:
        well = well[1:]
        max_letter = well[0] if well[0] > max_letter else max_letter
        max_number = int(well[1:]) if int(well[1:]) > max_number else max_number
    return max_letter, max_number

def build_plate_plot(cell_stats_complexity_df: pd.DataFrame, barcodes: Path, writeDir: Path, sampleName: str) -> list:
    barcodes_df = pd.read_csv(barcodes, sep="\t", names=["barcode", "well"])
    well_list = barcodes_df["well"].to_list()
    max_number = 1
    max_letter = 'A'
    datapane_list = []
    num_cells_dict = cell_stats_complexity_df["tgmt_well"].value_counts().to_dict()
    for well in well_list:
        well = well[1:]
        max_letter = well[0] if well[0] > max_letter else max_letter
        max_number = int(well[1:]) if int(well[1:]) > max_number else max_number
    
    for i in range(1,4):
        wellPlateTotalReadsDf = pd.DataFrame(0, columns=range(1, max_number+1), index=[chr(j) for j in range(65, ord(max_letter)+1)])
        wellPlateTotalCellsDf = wellPlateTotalReadsDf.copy()
        df_by_plate = cell_stats_complexity_df[cell_stats_complexity_df["tgmt_well"].str[0] == str(i)]
        df_by_plate = df_by_plate[["total", "tgmt_well"]]
        df_by_plate['tgmt_well'] = df_by_plate['tgmt_well'].str[1:]
        num_cells_dict = df_by_plate["tgmt_well"].value_counts().to_dict()
        df_by_plate = df_by_plate.groupby("tgmt_well", as_index=False)["total"].sum()
        
        for _, row in df_by_plate.iterrows():
            letter = row['tgmt_well'][0]
            number = int(row['tgmt_well'][1:])
            wellPlateTotalReadsDf.at[letter, number] = row["total"]
        
        for well in num_cells_dict:
            letter = well[0]
            number = int(well[1:])
            wellPlateTotalCellsDf.at[letter, number] = num_cells_dict[well]
        wellPlateTotalCellsDf.to_csv(writeDir / "csv" / f"{sampleName}.total_passing_cells_tn5_plate_{i}.csv")
        wellPlateTotalReadsDf.to_csv(writeDir / "csv" / f"{sampleName}.total_reads_tn5_plate_{i}.csv")
        datapane_list.append(plot_plate_plot(wellPlateTotalReadsDf, f"Total Reads(tn5): Plate {i}", 100))
        datapane_list.append(plot_plate_plot(wellPlateTotalCellsDf, f"Total Passing Cells(tn5): Plate {i}", False))
    
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


class BuildReadsPage:
    """
    Contains functions for building the reads page of either a sample report or a combined sample report
    """
    def __init__(self, sampleName: str, writeDir: Path, cell_stats_complexity_df: pd.DataFrame,
                 is_library_report: bool, fragment_df: pd.DataFrame = pd.DataFrame()):
        self.cell_stats_complexity_df = cell_stats_complexity_df
        self.sampleName = sampleName
        self.outDir = writeDir
        self.fragment_df = fragment_df
        self.is_library_report = is_library_report

    def build_knee_plot(self) -> dp.Plot:
        figs = []
        for sampleName in self.cell_stats_complexity_df["sampleName"].unique():
            df = self.cell_stats_complexity_df[self.cell_stats_complexity_df["sampleName"]==sampleName]
            sorted_uniq = np.sort(df['uniq'])[::-1]
            figs.append(go.Scatter(x=np.arange(1, len(df)+1), y=sorted_uniq, mode='lines', name=sampleName))
        if not self.is_library_report:
            figs.append(go.Scatter(x=[min(np.arange(1, len(self.cell_stats_complexity_df)+1)),
                                   max(np.arange(1, len(self.cell_stats_complexity_df)+1))],
                                   y=[self.cell_stats_complexity_df['threshold'].iloc[0],
                                   self.cell_stats_complexity_df['threshold'].iloc[0]],
                                   mode='lines', line=dict(dash="dash", color="green")))


        layout = go.Layout(xaxis=dict(type='log'), yaxis=dict(type='log'))
        fig = go.Figure(data=figs, layout=layout)
        fig.update_layout(
            title={'text': 'Knee Plot', 'x': 0.5, 'xanchor': 'center'},
            xaxis_type='log',
            yaxis_type='log',
            xaxis_title='Barcode',
            yaxis_title='Unique Reads',
            width=800,
            height=500,
            showlegend=self.is_library_report
        )
        fig.write_image(f"{self.outDir}/png/{self.sampleName}.kneePlot.png")
        return dp.Plot(fig)

    def build_fragment_length_histogram(self) -> dp.Plot:
        fig = px.bar(self.fragment_df, x='bin', y='count')
        fig.update_layout(
            title={'text': 'Fragment Length Histogram', 'x': 0.5, 'xanchor': 'center'},
            xaxis_title='Fragment Length',
            yaxis_title='Fragment Counts',
            autosize=False,
            width=800,
            height=500
        )
        fig.write_image(f"{self.outDir}/png/{self.sampleName}.fragLenHist.png")
        return dp.Plot(fig)

    def create_total_complexity_plot(self):
        passing_df = self.cell_stats_complexity_df[self.cell_stats_complexity_df["pass_filter"] == "pass"]
        passing_df['sampleName'] = passing_df['sampleName'] + "(passing)"
        passing_fig = px.scatter(passing_df, x='percent', y='uniq', color='sampleName', opacity=0.8)

        failing_df = self.cell_stats_complexity_df[self.cell_stats_complexity_df["pass_filter"]=="fail"]
        failing_df['sampleName'] = failing_df['sampleName'] + "(failing)"
        failing_fig = px.scatter(failing_df, x='percent', y='uniq', color='sampleName', opacity=0.2)

        fig = go.Figure(data=passing_fig.data + failing_fig.data)
        fig.update_layout(
            xaxis_title='(Uniq/TotalReads) %',
            yaxis_title='Unique Reads',
            yaxis_type='log',
            autosize=False,
            width=800,
            height=500
        )
        fig.write_image(f"{self.outDir}/png/{self.sampleName}.complexityTotal.png")
        return dp.Plot(fig)

    def build_threshold_table(self) -> dp.Table:
        threshold_dict = {}
        for sampleName in self.cell_stats_complexity_df["sampleName"].unique():
            threshold_dict[sampleName] = [int(self.cell_stats_complexity_df[
                self.cell_stats_complexity_df["sampleName"] == sampleName
            ]["threshold"].to_list()[0])]
        threshold = pd.DataFrame(threshold_dict)
        threshold = threshold.T.reset_index()
        threshold.columns = ['Metric', 'Value']
        return (DatapaneUtils.createTableIgnoreWarning(threshold[["Metric", "Value"]].style.pipe(
            DatapaneUtils.styleTable, title="Unique reads threshold for cell calling", hideColumnHeaders=True, boldColumn=['Metric'], numericCols=['Value'])))

    def build_combined_passing_met_stats_table(self, combined_passing_met: pd.DataFrame) -> dp.Table:
        combined_passing_met.iloc[0, combined_passing_met.columns.get_loc('Metric')] = "Number of Passing Cells"
        combined_passing_met.iloc[1, combined_passing_met.columns.get_loc('Metric')] = "Median Total Reads"
        combined_passing_met.iloc[2, combined_passing_met.columns.get_loc('Metric')] = "Median Unique Reads"
        combined_passing_met.iloc[3, combined_passing_met.columns.get_loc('Metric')] = "Median of Total Covered C's"
        combined_passing_met.iloc[4, combined_passing_met.columns.get_loc('Metric')] = "Median of Percent of Unique over Total Reads(%)"
        combined_passing_met.iloc[5, combined_passing_met.columns.get_loc('Metric')] = "Median of Covered CG's"
        combined_passing_met.iloc[6, combined_passing_met.columns.get_loc('Metric')] = "Median of Covered CH's"
        combined_passing_met.iloc[7, combined_passing_met.columns.get_loc('Metric')] = "Median of CG Methylation Percent(%)"
        combined_passing_met.iloc[8, combined_passing_met.columns.get_loc('Metric')] = "Median of CH Methylation Percent(%)"
        combined_passing_met.to_csv(f"{self.outDir}/{self.sampleName}.combinedPassingCellStats.csv", index=False)

    def build_bc_parser_stats(self, bcParserMetrics: Path | bool) -> dp.Table:
        if bcParserMetrics:
            with open(bcParserMetrics) as f:
                bc_metrics_dict = json.load(f)
            read_info = bc_metrics_dict['reads']
            bc_metrics_df = pd.DataFrame({
                'Status': ['Pass', 'BarcodeError', 'LinkerError', 'SequenceError'],
                'Reads': [read_info['Pass'][0], read_info['BarcodeError'][0], read_info['LinkerError'][0], read_info['SequenceError'][0]],
                'Percent': [read_info['Pass'][1], read_info['BarcodeError'][1], read_info['LinkerError'][1], read_info['SequenceError'][1]]
            })
        else:
            bc_metrics_df = pd.DataFrame({
                'Status': ['Pass', 'BarcodeError', 'LinkerError', 'SequenceError'],
                'Reads': [np.nan, np.nan, np.nan, np.nan],
                'Percent': [np.nan, np.nan, np.nan, np.nan]
            })
        return (DatapaneUtils.createTableIgnoreWarning(bc_metrics_df.style.pipe(
            DatapaneUtils.styleTable, title="Barcode Metrics", hideColumnHeaders=False, boldColumn=['Status'], numericCols=['Reads', 'Percent'])))

    def build_summary_stats(self, mapping_stats: pd.DataFrame | bool, trimming_stats: pd.DataFrame | bool) -> dp.Table:
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
        summary_stats.to_csv(f"{self.outDir}/csv/{self.sampleName}.summaryStats.csv", index=False)
        summary_stats.iloc[0, summary_stats.columns.get_loc('Metric')] = "Total Reads"
        summary_stats.iloc[1, summary_stats.columns.get_loc('Metric')] = "Reads Passing Trimming"
        summary_stats.iloc[2, summary_stats.columns.get_loc('Metric')] = "Reads Passing Mapping"
        fig = px.bar(summary_stats, y='Metric', x='Value', orientation='h')
        fig.update_layout(
            title={'text': 'Summary Stats', 'x': 0.5, 'xanchor': 'center'},
            xaxis_title='Reads',
            yaxis_title='',
            legend_title='',
            autosize=False,
            width=800,
            height=500,
            showlegend=self.is_library_report
        )
        return (dp.Plot(fig))

    def build_summary_stats_table(self, mapping_stats: pd.DataFrame | bool, trimming_stats: pd.DataFrame | bool) -> dp.Table:
        if isinstance(mapping_stats, bool) and isinstance(trimming_stats, bool):
            summary_stats = pd.DataFrame.from_dict({
                "Metric": ["total_reads", "percent_passing_trimming", "percent_passing_mapping", "percent_reads_in_passing_cells"],
                "Value": [np.nan, np.nan, np.nan, np.nan]
            }, dtype="object")
        else:
            summary_stats = pd.DataFrame.from_dict({
                "Metric": ["total_reads", "percent_passing_trimming", "percent_passing_mapping", "percent_reads_in_passing_cells"],
                "Value": [int(trimming_stats[trimming_stats["Metric"] == "total_reads"]["Value"].to_list()[0]), trimming_stats[trimming_stats["Metric"] == "percent_passing"]["Value"].to_list()[0],
                          (mapping_stats[mapping_stats["Metric"] == "mapped_reads"]["Value"].to_list()[0]/trimming_stats[trimming_stats["Metric"] == "total_reads"]["Value"].to_list()[0])*100,
                          (self.cell_stats_complexity_df[self.cell_stats_complexity_df['pass_filter'] == 'pass']['total'].sum()/trimming_stats[trimming_stats["Metric"] == "total_reads"]["Value"].to_list()[0])*100]
            }, dtype="object")
        summary_stats.iloc[0, summary_stats.columns.get_loc('Metric')] = "Total Reads"
        summary_stats.iloc[1, summary_stats.columns.get_loc('Metric')] = "Percent Reads Passing Trimming"
        summary_stats.iloc[2, summary_stats.columns.get_loc('Metric')] = "Percent Reads Passing Mapping"
        summary_stats.iloc[3, summary_stats.columns.get_loc('Metric')] = "Percent Reads In Passing Cells"
        summary_stats = summary_stats.drop(3)
        return (DatapaneUtils.createTableIgnoreWarning(summary_stats[["Metric", "Value"]].style.pipe(
            DatapaneUtils.styleTable, title="Summary Stats", hideColumnHeaders=True, boldColumn=['Metric'], numericCols=['Value'])))

    def construct_passing_cell_stats(self, met_passing: pd.DataFrame) -> dp.Table:
        met_stats = pd.DataFrame.from_dict({
            "Metric": ['number_passing_cells', 'median_total_reads', 'median_uniq_reads', 'totalcov_median',
                       'median_uniq_over_total_percent', 'CGcov_median', 'CHcov_median', 'CG_mC_Pct_median', 'CH_mC_Pct_median', 'sample_name'],
            "Value": [len(met_passing), int(met_passing['total'].median()), int(met_passing['uniq'].median()), int(met_passing['Coverage'].median()), met_passing['percent'].median(),
                      met_passing['CG_Cov'].median(), met_passing['CH_Cov'].median(), met_passing['CG_mC_Pct'].median(), met_passing['CH_mC_Pct'].median(), self.sampleName]
        }, dtype="object")
        met_stats.to_csv(f"{self.outDir}/csv/{self.sampleName}.passingCellSummaryStats.csv", index=False)
        met_stats.iloc[0, met_stats.columns.get_loc('Metric')] = "Number of Passing Cells"
        met_stats.iloc[1, met_stats.columns.get_loc('Metric')] = "Median Total Reads"
        met_stats.iloc[2, met_stats.columns.get_loc('Metric')] = "Median Unique Reads"
        met_stats.iloc[3, met_stats.columns.get_loc('Metric')] = "Median of Total Covered C's"
        met_stats.iloc[4, met_stats.columns.get_loc('Metric')] = "Median of Percent of Unique over Total Reads(%)"
        met_stats.iloc[5, met_stats.columns.get_loc('Metric')] = "Median of Covered CG's"
        met_stats.iloc[6, met_stats.columns.get_loc('Metric')] = "Median of Covered CH's"
        met_stats.iloc[7, met_stats.columns.get_loc('Metric')] = "Median of CG Methylation Percent(%)"
        met_stats.iloc[8, met_stats.columns.get_loc('Metric')] = "Median of CH Methylation Percent(%)"
        met_stats = met_stats.drop(9)
        return (DatapaneUtils.createTableIgnoreWarning(met_stats[["Metric", "Value"]].style.pipe(
            DatapaneUtils.styleTable, title="Passing Cells Stats Per Cell", hideColumnHeaders=True, boldColumn=['Metric'], numericCols=['Value'])))


class BuildMethylPage:
    """
    Contains functions for building the methyl page of either a sample report or a combined sample report
    """
    def __init__(self, sampleName: str, writeDir: Path, met_df: pd.DataFrame, is_library_report: bool):
        self.sampleName = sampleName
        self.outDir = writeDir
        self.met_df = met_df
        self.is_library_report = is_library_report

    def build_cell_covered_box(self) -> dp.Plot:
        IN = self.met_df[['sampleName', 'CG_Cov', 'CH_Cov']].melt(id_vars=['sampleName'], var_name='variable', value_name='value')
        fig = px.box(IN, x='variable', y='value', color='sampleName', points='suspectedoutliers')
        fig.update_layout(
            title={'text': 'Cell Cytosines Covered', 'x': 0.5, 'xanchor': 'center'},
            xaxis_title='',
            yaxis_title='Cytosines Covered',
            yaxis_type='log',
            legend_title='',
            autosize=False,
            width=800,
            height=500,
            showlegend=self.is_library_report,
            boxgroupgap = 0.5
        )
        fig.write_image(f"{self.outDir}/png/{self.sampleName}.cellCoveredBox.png")
        return dp.Plot(fig)
    
    def build_total_and_unique_reads_box(self) -> dp.Plot:
        tmp_df = self.met_df.rename(columns={"total": "Total", "uniq": "Unique"})
        IN = tmp_df[["sampleName", "Total", "Unique"]].melt(id_vars=["sampleName"], var_name="variable", value_name="value")
        fig = px.box(IN, x='variable', y='value', color='sampleName', points='suspectedoutliers')
        fig.update_layout(
            title={'text': 'Total & Unique Reads per Cell', 'x': 0.5, 'xanchor': 'center'},
            xaxis_title='',
            yaxis_title='Reads',
            legend_title='',
            autosize=False,
            width=800,
            height=500,
            showlegend=self.is_library_report,
            boxgroupgap = 0.5
        )
        fig.write_image(f"{self.outDir}/png/{self.sampleName}.cellTotalReadsBox.png")
        return dp.Plot(fig)

    def build_uniq_over_total_percent_box(self) -> dp.Plot:
        IN = self.met_df[["sampleName", "percent"]].melt(id_vars=["sampleName"], var_name="variable", value_name="value")
        fig = px.box(IN, x='variable', y='value', color='sampleName', points='suspectedoutliers')
        fig.update_layout(
            title={'text': 'Unique over Total Reads(%)', 'x': 0.5, 'xanchor': 'center'},
            xaxis_title='',
            yaxis_title='Unique/Total %',
            legend_title='',
            yaxis=dict(range=[0, 100]),
            autosize=False,
            width=800,
            height=500,
            showlegend=self.is_library_report,
            boxgroupgap = 0.5
        )
        fig.write_image(f"{self.outDir}/png/{self.sampleName}.cellUniqueOverTotalPercentBox.png")
        return dp.Plot(fig)

    def build_cg_cell_methyl_percent_box(self) -> dp.Plot:
        IN = self.met_df[["sampleName", "CG_mC_Pct"]].melt(id_vars=["sampleName"], var_name="variable", value_name="value")
        fig = px.box(IN, x='variable', y='value', color='sampleName', points='suspectedoutliers')
        fig.update_layout(
            title={'text': 'CG Cell Methylation', 'x': 0.5, 'xanchor': 'center'},
            xaxis_title='',
            yaxis_title='Cell Methylation %',
            legend_title='',
            yaxis=dict(range=[0, 100]),
            autosize=False,
            width=800,
            height=500,
            showlegend=self.is_library_report,
            boxgroupgap=0.5
        )
        fig.write_image(f"{self.outDir}/png/{self.sampleName}.cellMethylCgPercentBox.png")
        return dp.Plot(fig)

    def build_ch_cell_methyl_percent_box(self) -> dp.Plot:
        IN = self.met_df[["sampleName", "CH_mC_Pct"]].melt(id_vars=["sampleName"], var_name="variable", value_name="value")
        fig = px.box(IN, x='variable', y='value', color='sampleName', points='suspectedoutliers')
        fig.update_layout(
            title={'text': 'CH Cell Methylation', 'x': 0.5, 'xanchor': 'center'},
            xaxis_title='',
            yaxis_title='Cell Methylation %',
            legend_title='',
            yaxis=dict(range=[0, IN["value"].max() + 1]),
            autosize=False,
            width=800,
            height=500,
            showlegend=self.is_library_report,
            boxgroupgap = 0.5
        )
        fig.write_image(f"{self.outDir}/png/{self.sampleName}.cellMethylChPercentBox.png")
        return dp.Plot(fig)

    def build_cell_cg_per_total(self) -> dp.Plot:
        IN = self.met_df[['sampleName', 'cg_total_ratio']].melt(id_vars=["sampleName"], var_name="variable", value_name="value")
        fig = px.box(IN, x='variable', y='value', color='sampleName', points='suspectedoutliers')
        fig.update_layout(
            title={'text': 'CG per Read', 'x': 0.5, 'xanchor': 'center'},
            xaxis_title='',
            yaxis_title='CGs per Raw Read per Cell',
            legend_title='',
            yaxis=dict(range=[0, 1]),
            autosize=False,
            width=800,
            height=500,
            showlegend=self.is_library_report,
            boxgroupgap = 0.5
        )
        fig.write_image(f"{self.outDir}/png/{self.sampleName}.cellCGperTotalReadBox.png")
        return dp.Plot(fig)
