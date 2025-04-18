# ScaleBio tools
In addition to third-party and open-source software the workflow also uses executable tools developed by ScaleBio:

- [bc_parser](bc_parser.md) 
	- Extracts and error corrects cell-barcodes and UMIs from the original (input) fastq files
	- Splits (demultiplexes) the input fastq files into sample fastq files based on cell-barcodes (TN5)
	- Barcode and read-level metrics
- [sc_dedup](sc_dedup.md)
	- BAM deduplication aware of the multi-level combinatorial cell-barcodes.
 	- Barcode and read-level metrics.

## Installation
These tools are included in the `scaleMethyl` docker container image; when running the workflow with `-profile docker` (or another container engine, such as _singularity_, _podman_ etc.), they will be automatically available.

 These tool are however currently not available through Conda, so if running without a container system, they need to be installed first. A download script to get static pre-compiled binaries for linux (*x86_64*) is [included](../envs/download-scale-tools.sh); Simply run:

```/PATH/TO/ScaleMethyl/envs/download-scale-tools.sh```
 
  This will install the binaries in `ScaleMethyl/bin` (inside the Nextflow workflow directory), from where they will be available during workflow execution.
