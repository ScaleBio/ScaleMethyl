# Analysis parameters

All analysis parameters can be set in a `runParams.yml` file, which is then passed to nextflow with `nextflow run -params-file runParams.yml`. 
Alternatively each option is this file can also be set on the nextflow command-line directly, overwriting the value in the parameter file. E.g.
`nextflow run --samples=samples.foo.csv`

*Note* that `nextflow` options are given with a single `-`, while workflow parameters (e.g. `samples`) are given with a double dash `--`.


## Inputs
The workflow can either start from an Illumina sequencer runFolder (bcl files), a directory with fastq files or directories with bam files. Specify either

* fastqDir : "path/to/fastqs" <br>
OR
* runFolder : "path/to/runFolder" <br>
OR
* bam1Dir : "path/to/bams"
* bam2Dir : "path/to/bams"

where fastqDir is a directory containing all input fastq files (no Undetermined*fastq.gz). See [Fastq Generation](fastqGeneration.md) for details on file names, etc.

When starting from a sequencer run folder the workflow uses Illumina [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) for automatic fastq generation.

When starting from bam files refer to [mergeBam](mergeBam.md) for details on additional parameters, convention on filenames, etc. 

### Sample Information
* samples : "samples.csv"

A [file](examples/samples.csv) listing all samples in the analysis with their names, barcode sequences and optional sample settings

### Reference genome
* genome : "/genomes/grch38/genome.json"

Path to a [genome.json](docs/genomes.md) file that contains the location of all sequence and index files as well as other parameters for the reference genome to use. 

## Optional and Advanced parameters
Run `nextflow run path/to/ScaleMethyl --help` for a description of available options and see the example [runParams.yml](examples/runParams.yml) file.

System options (compute resource requirements, etc.) as well as all parameter defaults, are in the workflow [nextflow.config](../nextflow.config).

### Selected optional parameters
* Setting `bamOut` to false will suppress alignment (.bam file) output from bsbolt. Also setting `bamDedupOut` to false will suppress deduplicated bam file output from scDedup. If the user does not specifically need bam files for custom downstream analysis, disabling BAM output will save compute time and storage space
* `fastqOut` controls publishing fastq files from bcl-convert and bcParser to the output directory. As in the case of publishing bam files, if fastq files are not needed for custom analysis, disabling fastq output will save compute time and storage space
* `bclConvertParams` controls whether we split by lane or not

### Parallel execution
The workflow supports extra parallelism in a few different ways. One parameter is `splitFastq`, which controls whether samples are split by i5 barcode in the `bcl-convert` step. This results in n fastq files (per library/lane), where n is the number of i5 barcodes used in the experiment. There is another parameter, `splitBarcodeParsing`, which controls whether `bcParser` produces a set of demuxed and barcode corrected fastq files for each TN5 barcode. This results in 288 sets of fastq files which then go through alignment, deduplication and extraction in parallel.

Since our methylation workflow is quite compute intensive, our recommendation is to keep the defaults for `splitFastq` and `splitBarcodeParsing`(both of which are set to true by default). The caveat is that this configuration will launch around 288 jobs in parallel for each sample, so your infrastructure must be able to support that.

### Library Structure Definition
* libStructure : "lib.json"

The library structure JSON file defines 
* Where in the reads cell-barcodes are found
* What the list of allowed barcode sequences is
* Which parts of the reads represent genomic DNA and which should be masked (e.g. RT and ligation barcodes)

The default file, for our standard product configuration, is included in `references/lib.json`.
