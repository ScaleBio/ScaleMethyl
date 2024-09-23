# Outputs

Analysis workflow results for all samples in the run will appear in the output directory defined by `--outDir` (`ScaleMethyl.out` by default). 
For detailed information about the library and sample level QC reports see [qcReport.md](qcReport.md)

## Key output files
| Directory | File | Description |
|-----------|------|-------------|
| `report` | | QC reports and statistics |
| `report/sample_reports/<sample>` | `<sample>.report.html` | An interactive standalone HTML report including key metrics/figures for each sample |
|                  | `csv/<sample>.*.csv` | Sample metrics in csv format |
| `report/library_report/<library>` | `library.<libName>.report.html` | Barcode summary and demultiplexing statistics for the whole library (potentially multiple samples) |
|                  | `csv/<libName>.combinedPassingCellStats.csv` | Key metrics for passing cells for the whole library |
| `samples` | | Single-cell methylation outputs |
|           | `<sample>.allCells.csv` | Metrics per cell-barcode, including barcodes / well positions |
| `samples/genome_bin_matrix` | `<sample>.{CG,CH}.{score}.mtx.gz` | CG and CH binned genome-wide matrix files in Matrix Market format |
| `samples/methylation_coverage` | `<format>/<sample>/{CG,CH}/<barcode>.*` | Per-cell methylation calls in bismark .cov, .allc or amethyst .h5 format |
| `fastq` | | Fastq generation, QC and processing |
|         | `fastqc/*.html` | [fastqc](https://github.com/s-andrews/FastQC) report for each fastq file in the sequencing library |
| `barcodes` | `<sample>.<split>.demux/*.fastq.gz` | Demultiplexed sample fastq files; only included with `--fastqOut true` |
| `alignment` | | Outputs of Bisulfite alignment, split per sample / Tn5-barcode |
| `alignment/dedup` | `<sample>.<tn5>/*.bam` | Deduplicated alignment output, with single-cell barcode and UMI information in tags |

### Genome-bin Methylation Matrix
Rows are the genomic bins, and columns are the cellular barcode, see `features.tsv` and `barcodes.tsv` respectively.

Custom binned bed files for the CG and CH matrix generation can be passed using the `genomeTiles` and `genomeTilesCh` options in the [genomes.json](genomes.md); only included with `--matrixGenerationCG true` or `--matrixGenerationCH true` 

`cellGlobalMet` = Global methylation rate of the cell (methylated / total covered )

| File | Description |
|------|-------------|
| `<sample>.CG.score.mtx.gz` | `diff = (methylation rate - cellGlobalMet)`; if `diff` > 0, `score = diff/(1 - cellGlobalMet)`; else `score = diff/cellGlobalMet`. The score gives a scale of `-1:1` where -1 reflects no methylation for a covered region while non-covered regions are assigned to zero, which in some cases can give cleaner clustering. |
| `<sample>.CH.mtx.gz` | methylation rate (methylated / total covered CH) scale of `0:1` | 

These matrices can be loaded into Seurat, ScanPy or similar tools for clustering, downstream analysis, and visualization.
