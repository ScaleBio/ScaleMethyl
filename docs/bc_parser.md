# bc_parser

`bc_parser` is a ScaleBio developed tool for combinatorial barcode handling, error-correction and demultiplexing. It supports a wide variety of library designs and sequencing configurations.

# Executing bc_parser
For a complete list of options and flags that bc_parser uses, run `bc_parser --help`

## Inputs
* Sequencing reads (Fastq files) `--reads` or `--read1 --read2 ...`
    * Including index read fastqs if part of the library design
    * If fastq file names follow the standard `bcl-convert` naming convention (including `_R1` for read 1, etc.), the full set of files can be given with a single option e.g. `--reads S1*.fastq.gz`.
* Library structure definition (json) `--library`
    * This defines the location and sequence-lists for all barcodes

## Outputs
When run with the `splitFastq` option in `nextflow.config` set to true, bcParser will produce output files split by TN5 barcode. Since there are 288 TN5 barcodes, bcParser will produce 288 sets of output files, each organized according to the well coordinate for the corresponding TN5 barcode
That makes the output easily compatible with downstream analyses that require knowledge of the barcode when processing groups of reads. 

## Reads
`bc_parser` can output the processed reads and barcodes in multiple formats. The most common case is output of transformed fastq files per sample, with error-corrected barcodes in the read-name. This is activated with `--write-fastq`

The barcodes in the read-name will be a combination of all barcoding levels (e.g. tagmentation barcode and droplet barcode) that together define a unique cell. Specifically, for the scaleMethyl workflow, the read-names in fastq files will reflect the barcodes defined in the library.json file provided to the workflow. These follow the structure `i7 (10bp)+i5 (10bp)+tgmt (8bp)`

### Metrics
The output directory also contains an overall metrics file `metrics.json` as well as `.tsv` files with the counts of all observed barcode sequences

  
