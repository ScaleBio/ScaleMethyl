# Outputs

Analysis workflow results for all samples in the run will appear in the output directory defined by `--outDir` (`ScaleMethyl.out` by default). 
For detailed information about the library and sample level QC reports see [qcReport.md](qcReport.md)

## Key output files
| Directory | File | Description |
|-----------|------|-------------|
| `report` | `<sample>.report.html` | An interactive standalone HTML report including key metrics/figures for each sample |
| `report/csv` | `<sample>.allCells.csv` | Metrics per cell-barcode, including barcodes / well positions |
|              | `<sample>.total_{passing_cells,reads}_<tn5>_<plate_#>.csv` | Total passing cells and reads for each well coordinate on a sample level |
|              | `<sample>.mapping_stats.csv` | Sample-level mapping statistics derived from [BSBolt](https://github.com/NuttyLogic/BSBolt) alignment output  |
|              | `<sample>.trimming_stats.csv` | Sample-level read summary statistics after trimming  |
|              | `<sample>.passingCellSummaryStats.csv` | Sample-level summary statistics for all passing cells in sample   |
|              | `<sample>.fragment_hist.csv` | Read length histogram for a sample   |
| `library_report` | `library.<libName>.report.html` | Barcode summary and demultiplexing statistics for the whole library (potentially multiple samples) |
|                  | `<libName>.combinedPassingCellStats.csv` | Key metrics for passing cells for the whole library |
| `library_report/csv` | `<library>.total_{passing_cells,reads}_{i5,i7}.csv` | Total passing cells and reads for each well coordinate on a library level |
| `trim` | `<sample>.<tn5>/<sample>.<tn5>.trim_stats.json` | Detailed trimming reports for demultiplexed sample fastq files   |
| `fastqc` | `*.html` | [fastqc](https://github.com/s-andrews/FastQC) report for each fastq file in the sequencing library
| `barcodes/<libName>.<fastqNum>.demux` | `<sample>.*.fastq.gz` | Sample fastq files (Demultiplexed and barcode error-corrected); only included with `--fastqOut true` |
| `multiqc` | `<sample>.<tn5>*.multiqc_report.html` | [MultiQC](https://multiqc.info) report for fastq generation, fastQC and trimming
| `library_barcode_metrics` | `<libName>.metrics.json` | [bc_parser](bc_parser.md) metrics post demux and barcode correction 
| `alignment` | `<sample>.<tn5>/<sample>.<tn5>.bam` | [BSBolt](https://github.com/NuttyLogic/BSBolt) alignment output, including BAM file, with single-cell barcode and UMI information in tags for each sample.tn5; only included with `--bamOut true`
| `bamMerge` | `<sample>.<tn5>/<sample>.<tn5>.bam` | Merged BAM files with reads that have the same tn5(only if starting workflow from BAM). The @RG tags ID, SM, PL, and LB are assigned for each read; only included with `--bamMergeOut true`
| `bamDeDup` | `<sample>.<tn5>/<sample>.<tn5>.dedup.bam` | [sc_dedup](sc_dedup.md) Coordinate-sorted BAM files with duplicate reads removed.  Reads whose leftmost aligned fragment has the same leftmost position as a previously encountered read are considered duplicates; only included with `--bam-dedup-out = true`
| `matrix` | `<sample>.{CG,CH}.{score}.mtx.gz` | CG and CH binned genome-wide matrix files in Matrix Market format. Rows are the genomic bins, and columns are the cellular barcode, see `features.tsv` and `barcodes.tsv` respectively. Custom binned bed files for the CG and CH matrix generation can be passed using the `genomeTiles` and `genomeTilesCh` options in the [genomes.json](genomes.md); only included with `--matrixGenerationCG true` or `--matrixGenerationCH true`
| `cov` | `<barcode>.{CG,CH}.cov.gz` | Per-cell methylation calls in bismark .cov format
| `allc` | `{CG,CH}/<barcode>.allc.tsv.gz` | Per-cell methylation calls in bismark allCools .allc format; only included with `--allcOut` 

### Matrix Detailed Descriptions 

`cellGlobalMet` = Global methylation rate of the cell (methylated / total covered )

| File | Description |
|------|-------------|
| `<sample>.CG.score.mtx.gz` | `diff = (methylation rate - cellGlobalMet)`; if `diff` > 0, `score = diff/(1 - cellGlobalMet)`; else `score = diff/cellGlobalMet`. The score gives a scale of `-1:1` where -1 reflects missing values instead of also being assigned to zero, which in some cases can give cleaner clustering. |
| `<sample>.CH.mtx.gz` | methylation rate (methylated / total covered CH) scale of `0:1` | 

These matrices can be loaded into Seurat, ScanPy or similar tools for clustering, downstream analysis, and visualization.
