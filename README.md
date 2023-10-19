# ScaleBio scMethylation Workflow 

This is a Nextflow workflow to run analysis of ScaleBio single-cell methylation sequencing libraries. It processes data from sequencing reads to alignments, single-cell outputs (binned methylation cov matrix, etc.), and QC reports.

## Getting started
* First install [Nextflow](http://www.nextflow.io) (22.04 or later)
* Download this workflow to your machine
* Install [dependencies](docs/dependencies.md)
* Download / configure a reference [genome](docs/genomes.md) for your samples
* Create a [samples.csv](docs/samplesCsv.md) table for your samples
* Create [runParams.yml](docs/analysisParameters.md), specifying inputs and analysis options for your run
* Launch the workflow for your run

## Inputs
* Sequencing reads
    * Path to the Illumina Sequencer RunFolder (bcl files)
    * Path to the fastq files, generated outside (before) this workflow, see [Fastq generation](docs/fastqGeneration.md).
* Sample Table
    * A .csv file listing all samples in the analysis. See [samples.csv](docs/samplesCsv.md).
* Reference Genome
    * The workflow requires a reference genome, including a [BSBolt](https://github.com/NuttyLogic/BSBolt) index for alignment, and gene annotation. See [Reference Genomes](docs/genomes.md)

## Outputs
The workflow produces per-sample and per-library QC reports (`html`), alignments (`bam`), "insert" (`mtx`) and more; See [Outputs](docs/outputs.md) for a full list.


## Workflow Execution
### Workflow test
Add info on small test run

### Nextflow Command-line
**Note** that `nextflow` options are given with a single `-` (e.g. `-profile`), while workflow parameters (e.g. `--outDir`) are given with a double dash `--`.

See the [Nextflow command-line documentation](https://www.nextflow.io/docs/latest/cli.html) for the options to run `nextflow` on different systems (including HPC clusters and cloud compute).

## Configuration
### Specifying Analysis Parameters
Analysis parameters (inputs, options, etc.) can be defined either in a [runParams.yml](docs/examples/runParams.yml) file or directly on the nextflow command-line. See [analysisParameters](docs/analysisParameters.md) for details on the options.

### Config File
In addition to the analysis parameters, a user-specific nextflow configuration file can be used for system settings (compute and storage resources, resource limits, storage paths, etc.):

`-c path/to/user.config`

See [Nextflow configuration](https://www.nextflow.io/docs/latest/config.html) for the way different configuration files, parameter files and the command-line interact.

## Dependency Management
Different options to provide all required dependencies are described [here](docs/dependencies.md). Follow one approach there and then run nextflow with the corresponding `-profile`.

## Running in the cloud
Nextflow itself supports running using [AWS](https://www.nextflow.io/docs/latest/aws.html), [Azure](https://www.nextflow.io/docs/latest/azure.html) and [Google Cloud](https://www.nextflow.io/docs/latest/google.html). 

# Versions and Updates
See the [Change log](changelog.md)
