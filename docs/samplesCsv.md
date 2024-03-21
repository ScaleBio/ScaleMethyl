# Sample table

A sample table (e.g. [samples.csv](examples/samples.csv)) file is used to list the samples included in an analysis run, their sample barcode (TN5) sequences and optional sample-specific analysis parameters.

It is a comma separated file (csv), with a header line (column names), followed by one sample per line. 
The first column is required to be `sample` and contains the name of each sample. All other columns are optional. _Column names are case sensitive!_

 Column | Description | Example | Default
:---- | ---- | ---- | :----:
sample | Sample name | Foobar-2 |
barcodes | TN5-plate wells used for this sample (optional) | 1A01-3H12 | 1A01-3H12
libName | Name for the sequencing library / fastq files (optional) | ScaleMethyl | ScaleMethyl
libIndex | Subset of i7 barcodes used for this experiment (optional) | TTGATATGAA;CCTAAGCGGT | All barcodes in references/i7.txt
libIndex2 | Subset of i5 barcodes used for this experiment (optional) | TATCATGATC;GAGCATATGG | All barcodes in references/i5.txt
threshold | Unique reads number beyond which a cell will be considered as a passing cell (optional) | 1000 | Will be calculated based on a heuristic
split | Indicates whether bc_parser will produce output files split by TN5 barcode. Also controlled by the splitBarcodeParsing argument in nextflow config (optional) | true | true, because splitBarcodeParsing is set to true by default

* `sample` and `libName` should consist only of letters, numbers and dash (-)
* When running from pre-existing fastq file input, `libName` should match the first part of the fastq file name for this sample, e.g.: `Foo1` for `Foo1_*.fastq.gz`.
* If you have indexed on the PCR barcodes (i7 and or i5), you can pass the entire TN5 range `1A01-3H12` to `barcodes`
* When providing sequences in the `libIndex` and `libIndex2` columns, separate them using `;`
    * If `libIndex` and `libIndex2` are not provided, the assumption is that all the i5 and i7 barcodes have been used 

## Demultiplexing samples
During analysis the sequencing data is first converted into library fastq files (`libName` column). If multiple samples were included in one sequencing library, these are then demultiplexed based on the sample (TN5) barcodes. E.g.

sample | barcodes
-- | --
Foo | 1A01-3D12
Bar | 1E01-3H12

* The above is an equal split of tagmentation barcodes across 3 plates. The split is in the middle between row D and row E.

The TN5 wells used for each sample are given in `barcodes` as either
* An individual value (`1A01`)
* A range of wells (`1A01-3H12`)
    * Wells are sorted first by plate number, then by row and finally column, i.e. `1A01-3A12`.
        * Plate number ranges from 1-3, rows from A-H, and columns from 1-12 which indicate three plates for the TN5 barcodes (288 total)
    * Note that all ranges are read in **row-wise** order; e.g. 1A01-1B02, refers to A01-H01 (all of column 1 in plate 1) plus A02-B02 (2 wells; rows A and B column 2).
* A list of values or ranges, separated by semicolon (`;`) (`1A01;1A02-1H02`)
* If you load a sample by columns, seperate the tagment barcode locations with a (`;`)
     * e.g. Column 1 of plate 2 = (`2A01;2B01;2C01;2D01;2E01;2F01;2G01;2H01`)
