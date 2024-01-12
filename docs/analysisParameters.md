# Analysis parameters

All analysis parameters can be set in a `runParams.yml` file, which is then passed to nextflow with `nextflow run -params-file runParams.yml`. 
Alternatively each option is this file can also be set on the nextflow command-line directly, overwriting the value in the parameter file. E.g.
`nextflow run --samples=samples.foo.csv`

*Note* that `nextflow` options are given with a single `-`, while workflow parameters (e.g. `samples`) are given with a double dash `--`.


## Inputs
The workflow can start from an Illumina sequencer runFolder (bcl files), a directory with fastq files or directories with bam files. Specify either

* `runFolder : "path/to/runFolder"` <br>
OR
* `fastqDir : "path/to/fastqs"` <br>
OR
* `bam1Dir : "path/to/bams"`  
* `bam2Dir : "path/to/bams"`  

where fastqDir is a directory containing all bcl-converted fastq files (no `Undetermined*fastq.gz`). See [Fastq Generation](fastqGeneration.md) for details on file names, etc.

When starting from a sequencer run folder the workflow uses Illumina [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) for automatic fastq generation.

When starting from bam files refer to [mergeBam](mergeBam.md) for details on additional parameters, convention on filenames, etc. 

### Sample Information
* `samples : "samples.csv"`

A [file](examples/samples.csv) listing all samples in the analysis with their names, barcode sequences and optional sample settings

### Reference genome
* `genome : "/genomes/grch38/genome.json"`

Path to a [genome.json](docs/genomes.md) file that contains the location of all sequence and index files as well as other parameters for the reference genome to use. 

## Optional and Advanced parameters
Run `nextflow run path/to/ScaleMethyl --help` for a description of available options and see the example [runParams.yml](examples/runParams.yml) file.

System options (compute resource requirements, etc.) as well as all parameter defaults, are in the workflow [nextflow.config](../nextflow.config).

### Optional and Advanced Parameters

#### Fastq handling 
* `bclConvertParams` controls whether we split by lane or index2, see [Fastq generation](docs/fastqGeneration.md).

#### Optional intermediate file outputs 
* `fastqOut` controls publishing fastq files from bcl-convert and bcParser barcode corrected sample demultiplexed fastqs to the output directory. If fastq files are not needed for custom analysis, disabling fastq output will save compute time and storage space. 
* `trimOut` controls publishing post barcode demultiplexed trimmed fastq files to the output directory in the `trim` folder. 
* Setting `bamOut` to false will suppress alignment (.bam file) output from bsbolt. Also setting `bamDedupOut` to false will suppress deduplicated bam file output from scDedup. If the user does not specifically need bam files for custom downstream analysis, disabling BAM output will save compute time and storage space.
* `bamMergeOut` controls publishing of merged bam files to the output directory in the `bamMerge`, see [merge]
* `covOut` controls publishing methylation extraction bed files to the output direcotry `{cg,ch}_sort_cov` folders. These are split by sample and can be parsed for passing cells and reformatted for downstream analysis tools. 

#### Cell calling 
The default cell thresholding algorithm used by the ScaleMethyl workflow imposes a one-dimensional threshold on the coordinate collapsed uniquely mapped reads per cell barcode.

Parameters:

* All cell-barcodes with over `minUniqCount` [1000] unique reads are considered as possible cells 
* The read-count of top cells, topCount, is estimated as the `topCellPercentile` [99] of read-counts of cell-barcodes above `minUniqCount` 
* The cell threshold is set a fixed `minCellRatio` [20] of the `topCount`; `topCount` / `minCellRatio` 

You can adjust `--minUniqCount`, `--topCellPercentile`, `--minCellRatio` in the runParams.yml. These can be overridden by providing the threshold column in the [samples.csv](docs/samplesCsv.md).

#### Filter cells by percent unique reads per barcode
In the sample reports, you will find our plot of unique reads vs percent unique reads compared (compared to total reads per barcode) `report/<sample>/png/<sample>.complexityTotal.png`. To filter passing cells for a range in the percent unique reads (x-axis) in addition to the unique read threshold, provide these options in addition to cell calling options above:
* `--minUniqTotal` [1] minimum percent unique reads per barcode for cell calling. 
* `--maxUniqTotal` [100] maximum percent unique reads per barcode for cell calling.

#### Rerunning report generation with a different threshold

This assumes that you have ran the workflow and would like to apply a different threshold than being called by the cell thresholding algorithm.

If you resume the run with a samples.csv that you have added the threshold to, it will rerun the workflow from the beginning. However, these matrixes can be easily filtered for the passing cells if needed.

If you do not wish to rerun only the report with a more stringent threshold and can filter your matrices for downstream analysis, there is a reporting only entry point to the workflow. 

* `reporting : true` option to only run reporting from a `resultDir` for a completed run. 
* `resultDir` would be the outDir of your original run. 
* `outDir` is your new outDir for the new reports. It is recommended that you do not use the same folder as the original run, as you will need to reference your originals for the run summary stats. 
* `samples.csv` must be provided with the additional threshold column, see [samples.csv](docs/samplesCsv.md) for details.


### Parallel execution
The workflow supports extra parallelism in a few different ways. One parameter is `splitFastq`, which controls whether samples are split by i5 barcode in the `bcl-convert` step rather than lane. This results in n fastq files (per library/lane), where n is the number of i5 barcodes used in the experiment. This allows the workflow to split up the demultiplexing step more than 1-4 lanes. To use `splitFastq`, you must also turn off the lane splitting by default with `bclConvertParams : "--no-lane-splitting true"`. `splitFastq` and no-lane-splitting options are both on by default. 

To run the workflow splitting by lane and rather than i5 set the following
* `splitFastq : false`  
and 
* `bclConvertParams : "--no-lane-splitting false"`  

There is another parameter, `splitBarcodeParsing`, which controls whether `bcParser` produces a set of demuxed and barcode corrected fastq files for each TN5 barcode. This results in 288 sets of fastq files which then go through alignment, deduplication and extraction in parallel.

Since our methylation workflow is quite compute intensive, our recommendation is to keep the defaults for `splitFastq` and `splitBarcodeParsing`(both of which are set to true by default). The caveat is that this configuration will launch around 288 jobs in parallel for each sample, so your infrastructure must be able to support that.

### Library Structure Definition
* `libStructure : "lib.json"`

The library structure JSON file defines 
* Barcode locations in paired-end reads Where in the reads cell-barcodes are found
* List of allowed barcode sequences for [tagmentation](`references/tgmt.txt`) barcodes, [i7](`references/i7.txt`) and [i5](`references/i5.txt`) MET PCR barcodes.
* Which parts of the reads represent genomic DNA and which should be masked (e.g. PCR barcodes, mosaic end sequences, 10 base randomer)

The default file, for our standard product configuration, is included in `references/lib.json`.

Previous alpha kit PCR indexes can be accessed from `references/prev_pcr/lib.json` with `libStructure : "prev_pcr/lib.json"` in the runParams.yml.

