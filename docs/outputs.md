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
|              | `<sample>.passingCellMapMethylStats.csv` | Cell-level methylation site coverage statistics for all cells in a sample   |
|              | `<sample>.passingCellSummaryStats.csv` | Sample-level summary statistics for all passing cells in sample   |
|              | `<sample>.fragment_hist.csv` | Read length histogram for a sample   |
| `library_report` | `library.<libName>.report.html` | Barcode summary and demultiplexing statistics for the whole library (potentially multiple samples) |
|                  | `<libName>.combinedPassingCellStats.csv` | Key metrics for passing cells for the whole library |
| `library_report/csv` | `<library>.total_{passing_cells,reads}_{i5,i7}.csv` | Total passing cells and reads for each well coordinate on a library level |
| `trim` | `<sample>.<tn5>/<sample>.fq.gz_trimming_report_R{1,2}.txt` | Detailed trimming reports for demultiplexed sample fastq files   |
| `fastqc` | `*.html` | [fastqc](https://github.com/s-andrews/FastQC) report for each fastq file in the sequencing library
| `fastq/<libName>.<laneNum>.demux` | `<sample>.*.fastq.gz` | Sample fastq files (Demultiplexed and barcode error-corrected); only included with `--fastqOut true` |
| `multiqc` | `<sample>.<tn5>*.multiqc_report.html` | [MultiQC](https://multiqc.info) report for fastq generation, fastQC and trimming
| `library_barcode_metrics` | `<libName>.metrics.json` | [bc_parser](bc_parser.md) metrics post demux and barcode correction 
| `alignment` | `<sample>.<tn5>/<sample>.<tn5>.bam` | [BSBolt](https://github.com/NuttyLogic/BSBolt) alignment output, including BAM file, with single-cell barcode and UMI information in tags for each sample.tn5; only included with `--bamOut true`
| `bamMerge` | `<sample>.<tn5>/<sample>.<tn5>.bam` | Merged BAM files with reads that have the same tn5(only if starting workflow from BAM). The @RG tags ID, SM, PL, and LB are assigned for each read; only included with `--bamMergeOut true`
| `bamDeDup` | `<sample>.<tn5>/<sample>.<tn5>.dedup.nsrt.bam` | [sc_dedup](sc_dedup.md) Name-sorted BAM files with duplicate reads removed.  Reads whose leftmost aligned fragment has the same leftmost position as a previously encountered read are considered duplicates; only included with `--bam-dedup-out = true`
| `matrix` | `<sample>.{CG,CH}.{mtx,cov.mtx,score.mtx,ratio.mtx}` | [sciMETv2](https://github.com/adeylab/sciMETv2/blob/main/sciMET_meth2mtx.pl) CG & CH binned genome-wide matrix files. Rows are the genomic bin coordinate, and columns are the cellular barcode. Custom binned bed files for the matrix generation can be passed in the `genomeTiles` in the [genomes.json](genomes.md); only included with `--matrixGenerationCG true` or `--matrixGenerationCH true`
| `{cg,ch}_sort_cov` | `<sample>.<tn5>{CG,CH}.chroms.sort/chr*.bed.gz` | Bed sorted raw methylation calls for each barcode. Note these are not filtered for passing cells; only included with `--covOut true` 

### Matrix Detailed Descriptions 

`cellGlobalMet` = Global methylation rate of the cell (methylated / total covered )

| File | Description |
|------|-------------|
| `mtx` | methylation rate (methylated / total covered CG or CH) scale of `0:1`|
| `score.mtx` | `diff = (methylation rate - cellGlobalMet)`; if `diff` = 0, `score = diff/(1 - cellGlobalMet)`; else `score = diff/cellGlobalMet`. The score gives a scale of `-1:1` where -1 reflects missing values instead of also being assigned to zero, which in some cases can give cleaner clustering. 
| `cov.mtx` | total covered CG or CH in any methylation context | 
| `ratio.mtx` | methylation rate / `cellGlobalMet` | 

The binned methylation rate matrix (.mtx) can be loaded into Seurat, ScanPy or similar tools for visualization or downstream analysis.
