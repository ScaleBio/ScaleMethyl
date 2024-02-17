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
