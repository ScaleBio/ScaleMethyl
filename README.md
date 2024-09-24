# ScaleBio Seq Suite: Methylation Workflow 

This is a Nextflow workflow to run analysis of ScaleBio Single Cell DNA methylation sequencing libraries. It processes data from sequencing reads to alignments, single-cell outputs (methylation coverage, binned genome matrix, etc.), and QC reports.

## Getting started
* First install [Nextflow](http://www.nextflow.io) (23.10 or later)
* Download this workflow to your machine
* Setup [dependencies](docs/dependencies.md)
* Launch the small pipeline [test run](#workflow-test)
* Download / configure a reference [genome](docs/genomes.md) for your samples
* Create a [Sample Barcode Table](docs/samplesCsv.md)
* Create [runParams.yml](docs/analysisParameters.md), specifying inputs and analysis options for your run
* Launch the workflow for your data

## Requirements
* Linux system with GLIBC >= 2.17 (such as CentOS 7 or later)
* Java 11 or later
* 64GB of RAM and 12 CPU cores
    * For large datasets a distributed compute system (HPC or cloud) is strongly recommended
* Working storage space 5x the input data size
    * E.g. Approximately 4TB temporary storage for deep sequencing of a _small kit_

## Inputs
* Sequencing reads, either
    * `--runFolder`: Path to the Illumina Sequencer RunFolder (`bcl` files)
    * `--fastqDir`: Path to the fastq files, generated outside (before) this workflow; see [Fastq generation](docs/fastqGeneration.md).
* Sample Barcode Table
    * `--samples`: A .csv file listing all samples in the analysis; See [samples.csv](docs/samplesCsv.md).
* Reference Genome
    * `--genome`: A reference genome, including a [BSBolt](https://github.com/NuttyLogic/BSBolt) index for alignment, and gene annotation; See [Reference Genomes](docs/genomes.md)

## Outputs
The workflow produces per-sample and per-library QC reports (`html`), alignments (`bam`), per-cell methylation coverage calls (bismark-like `cov`), genomic-bin methylation score matrix files (`mtx`) and more; See [Outputs](docs/outputs.md) for a full list.

## Workflow Execution
### Workflow test
A small test run, with all input data stored online can be run with the following command:

`nextflow run /PATH/TO/ScaleMethyl -profile PROFILE -params-file /PATH/TO/ScaleMethyl/docs/examples/runParams.yml --outDir ScaleMethyl.out`

`-profile docker` is the preferred option if the system supports _Docker_ containers;  See [Dependency Management](#dependency-management) for alternatives.

With this command, nextflow will automatically download the example data from the internet (AWS S3), so please ensure that the compute nodes have internet access and storage space. Alternatively you can manually download the data first (using [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html))
```
aws s3 sync s3://scale.pub/testData/methylation/reference/ reference --no-sign-request
aws s3 sync s3://scale.pub/testData/methylation/downsampled_pbmcs_v1.1/ fastqs --no-sign-request
```
and then run with
`nextflow run /PATH/TO/ScaleMethyl/ -profile PROFILE -params-file /PATH/TO/ScaleMethyl/docs/examples/runParams.yml --genome reference/genome.json --fastqDir fastqs --outDir ScaleMethyl.out

Note that this test run is merely a quick and easy way to verify that the pipeline executes properly and does not represent a complete dataset.


### Nextflow
See the [Nextflow command-line documentation](https://www.nextflow.io/docs/) for the options to run `nextflow` on different systems; including HPC clusters and cloud compute.

## Configuration
### Specifying Analysis Parameters
Analysis parameters (inputs, options, etc.) can be defined either in a [runParams.yml](docs/examples/runParams.yml) file or directly on the nextflow command-line. See [analysisParameters](docs/analysisParameters.md) for details on the options.

**Note** that `nextflow` options are given with a single `-` (e.g. `-profile`), while workflow parameters (e.g. `--outDir`) are given with a double dash `--`.

### Config File
In addition to the analysis parameters, a user-specific nextflow configuration file can be used for system settings (compute and storage resources, resource limits, storage paths, etc.):

`-c path/to/user.config`

See [Nextflow configuration](https://www.nextflow.io/docs/latest/config.html) for the way different configuration files, parameter files and the command-line interact.


## Dependency Management
The Nextflow workflow can automatically use pre-built docker containers with all dependencies included. Activating the included `-profile docker` enables the required Nextflow settings. For details and alternatives see [Dependencies](docs/dependencies.md).

## Running in the cloud
Nextflow itself supports running using [AWS](https://www.nextflow.io/docs/latest/aws.html), [Azure](https://www.nextflow.io/docs/latest/azure.html) and [Google Cloud](https://www.nextflow.io/docs/latest/google.html). 

In addition [Nextflow tower](https://tower.nf) offers another simple way to manage and execute nextflow workflows in Amazon AWS.

# Versions and Updates
See the [change log](changelog.md)

# License
By purchasing product(s) and downloading the software product(s) of ScaleBio, You accept all of the terms of the [License Agreement](LICENSE.md). If You do not agree to these terms and conditions, You may not use or download any of the software product(s) of ScaleBio.

