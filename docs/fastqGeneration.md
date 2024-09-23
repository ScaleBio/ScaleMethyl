# Fastq generation
The workflow can be started from a sequencer runFolder (bcl) [`--runFolder`]. In that case fastq files will be generated in the workflow, using Illumina `bcl-convert`.

Alternatively, it is also possible to generate fastq files upstream, for example when the ScaleMethyl library is multiplexed together with other libraries during sequencing at a core facility. In that case we recommend using [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) version 3.9 or later.

## samplesheet.csv
An example samplesheet for bcl-convert [ScaleMethyl_small-kit_samplesheet.csv](examples/samplesheets/ScaleMethyl_small-kit_samplesheet.csv) with typical options is included. This generates one set of fastq files (R1, R2) for each well of the PCR (final) plates. 
**Note** that this is not a split per sample loaded into the Tagmentation (first) plate.

The [reverse complement](examples/samplesheets/ScaleMethyl_small-kit_samplesheet_revComp.csv) version of the samplesheet.csv can be used with [instruments](https://knowledge.illumina.com/software/general/software-general-reference_material-list/000001800) that require index2 to be specified in the opposite orientation.

The equivalent for the large kit (8 plates): [ScaleMethyl_large-kit_samplesheet.csv](examples/samplesheets/ScaleMethyl_large-kit_samplesheet.csv)

## Index reads
For ScaleBio Methylation libraries the tagmentation barcodes are included in read 2, while the MET i5 and i7 PCR barcodes are in index reads 1 and 2. Since the index reads are needed for demultiplexing and barcode correction, we need to tell `bcl-convert` to generate index read fastqs using the `samplesheet.csv` setting: \
`CreateFastqForIndexReads,1`

# Using pre-generated fastq files as workflow input
Set `--fastqDir` to the directory containing the fastq files. The input fastq file names should follow the pattern `<LibraryName>_..._<Read>_...fastq.gz`, where
* `LibraryName` is the name of the sequencing library, which needs to match the `libName` column in `samples.csv`
* `Read` is one of `R1`, `R2`, `I1`, `I2`

This can be achieved by matching the `sample_ID` in the Illumina _bcl-convert_ _samplesheet.csv_ to the `libName` in the ScaleBio `samples.csv`.