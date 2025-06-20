# Version 1.3
## 1.3.0
* Provides updated Amethyst output for v1.0.0, update your Amethyst to 1.0.0 to use final output.
* Adds support for BWA-meth/parabricks aligner and sets BWA-meth as new default, providing improved runtime performance and filtering spurious low-coverage methylation sites. 
* Numerous reporting improvements, including UMAP generation
# Version 1.2
## 1.2.3 
* Fix to Amethyst output (CG and CH context was being output switched) 
    - Introduced in v1.2.0 and affecting runs with --amethystOut True
    - Can use --startPostExtraction  True and --previousOutDir options to fix without re-running the entire workflow
## 1.2.2
* Fix to sample merging (especially when calculating CH contexts)
* Using local tempdir to avoid permission errors in reporting
## 1.2.1
* Fix remote URL in sample and library QC reports (HTML)
## 1.2.0
* Analysis workflow support for individual plates & subsequent plate merging (docs/mergeMultipleAlignments.md)
* Window-size input for matrix gen, rather than file-format upload (including blacklist)
* Post-alignment checkpoint - allows for quick re-run of matrix gen
* Optional CH context calculations
* allcools files are tabix indexed automatically
* Reporting updates

# Version 1.1.0
* Added Amethyst compatible .h5 cell-methylation outputs
* Changed [matrix format](docs/outputs.md#matrix-detailed-descriptions).
* Filter reads/fragments with artifically high CH methylation signal 
* Updated QC metrics and sample QC (HTML) [report](docs/qcReport.md)
* Updated [output directory structure](docs/outputs.md)
* Updated per-cell methylation signal extraction to handle indel errors in reads as well as improve runtime performance
* Combined all depenencies in a single conda environment / docker container
* Use cutadapt directly for adapter trimming

# Version 1.0.4
* Added optional [single-cell output formats](docs/analysisParameters.md#optional-intermediate-file-outputs): ALLC and Bismark coverage files

# Version 1.0.3
* Upgraded datapane version and downgraded urllib version in containers

# Version 1.0.2
* Updated the testing dataset including examples for runParams.yml, samples.csv and samplesheet.csv's
* Option for TSS enrichment `--runTssEnrich`
* Add option to pass a different bin size to CG and CH matrix generation
    - See `tssWin` and `backgroundWin` [genomes.md ](docs/genomes.md)

# Version 1.0.1
* Add TSS enrichment to sample QC reports
    - See `genomeTiles` and `genomeTilesCh` [genomes.md ](docs/genomes.md)
* Add options for filtering passing cells by a percent uniquely mapped reads
    - See `--minUniqTotal` and `--maxUniqTotal`
* Add options to turn off CG or CH matrix generation 
    - See `--matrixGenerationCG` and `--matrixGenerationCH`
* Add library structure for splitting by i7 rather than i5 during bcl-convert. See [lib_split_on_i7.json](references/lib_split_on_i7.json)
* Add multiqc for trimmed fastqs
* README updates

# Version 1.0
Initial Release
