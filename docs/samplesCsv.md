# Sample Barcode Table

A sample barcode table (e.g. [samples.csv](examples/samples.csv)) file is used to list the samples included in an analysis run, their sample barcode (TN5) sequences and optional sample-specific analysis parameters.

It is a comma separated file (csv), with a header line (column names), followed by one sample per line. 
The first column is required to be `sample` and contains the name of each sample. All other columns are optional. _Column names are case sensitive!_

 Column | Description | Example | Default
:---- | ---- | ---- | :----:
sample | Sample name | Foobar-2 |
barcodes | TN5-plate wells used for this sample (optional) | 1A01-1H12;2A01-2D12 | 1A01-3H12
libName | Name for the sequencing library / fastq files (optional) | ScaleMethyl | ScaleMethyl
libIndex | Subset of i7 barcodes used for this experiment (optional) | TTGATATGAA;CCTAAGCGGT | All barcodes in references/i7.txt
libIndex2 | Subset of i5 barcodes used for this experiment (optional) | TATCATGATC;GAGCATATGG | All barcodes in references/i5.txt
threshold | Unique reads number beyond which a cell will be considered as a passing cell (optional) | 1000 | Will be calculated based on a heuristic
split | Split sample fastq files by TN5 barcode during analysis for parallelization. Also controlled by the `splitFastq` workflow argument (see `nextflow.config`). (optional) | true | true (value of `splitFastq` parameter)

* `sample` and `libName` should consist only of letters, numbers and dash (`-`) and start with a letter.
* A single `libName` should be used for all samples in a Scale Methylation sequencing library; not a different one per sample.
    * When running from fastq file input, `libName` should match the first part of the fastq file name for this sample, e.g.: `ScaleMethyl` for `ScaleMethyl_*.fastq.gz`.

## Demultiplexing samples
During analysis the sequencing data is first converted into library fastq files (`libName` column). If multiple samples were included in one sequencing library, these are then demultiplexed based on the sample (TN5) barcodes.

The TN5 wells used for each sample are given in `barcodes` as either
* An individual value (`1A01`); meaning Plate1, Row A, Column 1
    * Plate number ranges from 1-3, rows from A-H, and columns from 1-12 (288 total)
* A range of wells (`1A01-3H12`)
    * Ranges are interpreted first by plate then in **row-wise** order
    * `1A01-1B12` refers to rows A and B of the first Tn5 barcode plate
    * `1A01-2D12` refers to the entire first plate and the first 4 rows of the second plate
* A list of values or ranges, separated by semicolon (`1A01;1A02-1H02`)
* If samples were loaded by columns, list all wells individually
     * e.g. Column 1 of plate 2 = (`2A01;2B01;2C01;2D01;2E01;2F01;2G01;2H01`)

### Example
sample | barcodes
-- | --
One | 1A01-2D12
Two | 2E01-2H12
Three | 3A01-3A06;3B01-3B06;3C01-3C06;3D01-3D06;3E01-3E06;3F01-3F06;3G01-3G06;3H01-3H06
Four | 3A07-3A12;3B07-3B12;3C07-3C12;3D07-3D12;3E07-3E12;3F07-3F12;3G07-3G12;3H07-3H12

* *Sample One:* The first TN5 barcode plate and the first half of the 2nd TN5 plate split row-wise
* *Sample Two:* The last half of the 2nd TN5 barcode plate, after row D
* *Sample Three:* The first half of the 3rd TN5 barcode plate split column-wise (columns 1-6)
* *Sample Four:* The last half of the 3rd TN5 barcode plate split column-wise (columns 7-12)


