# Outputs

Analysis workflow results for all samples in the run will appear in the output directory defined by `--outDir` (`ScaleMethyl.out` by default). 
For detailed information about the library and sample level QC reports see [qcReport.md](qcReport.md)

## All output files
Here is a description of all of the output files that are saved outside the work directory. Most users will only be interested in the reports in the `report` folder and the output files in the `samples` folder for downstream analysis.
| Directory | File | Description |
|-----------|------|-------------|
| `alignments` |  | Outputs of Bisulfite alignment, split per sample / Tn5-barcode |
| `alignments/dedup/<sample>.<tgmt_well>` | `<sample>.<tgmt_well>.dedup.bam` | Deduplicated alignment output, with single-cell barcode and UMI information in tags; Can be turned off with `--bamDedupOut false` |
|                                         | `<sample>.<tgmt_well>.cell_stats.tsv` | Tab delimited report of reads for each barcode (cell) in tgmt_well for this sample |
|                                         | `<sample>.<tgmt_well>.cellInfo.txt` | Tab delimited report of coverage (CG, CH, and overall) for each barcode (cell) in tgmt_well for this sample |
|                                         | `<sample>.<tgmt_well>.dedup_stats.tsv` | Report of overall dedup stats for tgmt_well for this sample |
|                                         | `<sample>.<tgmt_well>.frament_hist.tsv` | Data for fragment size histogram for tgmt_well for this sample |
| `barcodes` | | Outputs of barcode reports and fastq files |
| `barcodes/<sample>.<split>.demux` | `*.fastq.gz` | Demultiplexed sample fastq files; only included with `--fastqOut true` |
|                                   | `barcode_counts.i5.tsv` | Counts of detected i5 barcodes in this split |
|                                   | `barcode_counts.i7.tsv` | Counts of detected i7 barcodes in this split |
|                                   | `barcode_counts.tgmt.tsv` | Counts of detected combined i5/i7 counts (tgmt_well) barcodes in this split |
| `fastq` | | Fastq generation, QC and processing |
|         | `fastqc/*.html` | [fastqc](https://github.com/s-andrews/FastQC) report for each fastq file in the sequencing library |
|         | `trim/<sample>.<tgmt_well>/<sample>.<tgmt_well>.trim_stats.json` | JSON containing stats from trimming |
| `matrix_generation` | | Files needed for matrix generation step (Required for --startPostExtraction  option) |
| `metrics_for_reporting` |  | By-sample csv report information. (Required for --reportingOnly option) |
| `report` | | QC reports and statistics |
| `report/sample_reports/<sample>` | `<sample>.report.html` | An interactive standalone HTML report including key metrics/figures for each sample |
|                  | `csv/<sample>.*.csv` | Sample metrics in csv format |
|                  | `csv/<sample>.report_clusters.tsv` | cluster IDs for each passing barcode for umap in sample report
| `report/library_report/<library>` | `library.<libName>.report.html` | Barcode summary and demultiplexing statistics for the whole library (potentially multiple samples) |
|                  | `csv/<libName>.combinedPassingCellStats.csv` | Key metrics for passing cells for the whole library |
| `samples` | | Single-cell methylation outputs |
|           | `<sample>.allCells.csv` | Metrics per cell-barcode, including barcodes / well positions |
| `samples/genome_bin_matrix` | `<sample>.<context>.{score.}mtx.gz` | Context (CG/CH) binned genome-wide matrix files in (gzipped) Matrix Market format. CG context includes score. |
|                             | `<sample>.<context>.features.tsv` | Contains genome coordinates for score matrix for each context |
|                             | `<sample>.barcodes.tsv` | Contains full i5/i7/cell barcodes for passing cells |
| `samples/methylation_coverage` | `<format>/<sample>/{CG,CH}/<barcode>.*` | Per-cell methylation calls in bismark .cov, .allc or amethyst .h5 format |



### Genome-bin Methylation Matrix
Rows are the genomic bins, and columns are the cellular barcode, see `features.tsv` and `barcodes.tsv` respectively.

Custom binned bed files for the CG and CH matrix generation can be passed using the `genomeTiles` and `genomeTilesCh` options in the [genomes.json](genomes.md); only included with `--windowMatrixOut true`. CH context is also only calculated if --calculateCH is true.

`cellGlobalMet` = Global methylation rate of the cell (methylated / total covered )

| File | Description |
|------|-------------|
| `<sample>.CG.score.mtx.gz` | `diff = (methylation rate - cellGlobalMet)`; if `diff` > 0, `score = diff/(1 - cellGlobalMet)`; else `score = diff/cellGlobalMet`. The score gives a scale of `-1:1` where -1 reflects no methylation for a covered region while non-covered regions are assigned to zero, which in some cases can give cleaner clustering. |
| `<sample>.CH.mtx.gz` | methylation rate (methylated / total covered CH) scale of `0:1` | 

These matrices can be loaded into Seurat, ScanPy or similar tools for clustering, downstream analysis, and visualization.