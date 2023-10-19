# Fastq generation
The workflow can be started from a sequencer runFolder (bcl) [`--runFolder`]. In that case fastq files will be generated internally, using Illumina `bcl-convert`.

Alternatively it is also possible to generate fastq files upstream, for example when the ScaleMethyl library is multiplexed together with other libraries during sequencing at a core facility. In that case we recommend using [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html), but older versions of the Illumina software are also possible.

An example [samplesheet.csv](examples/samplesheet.csv) with typical options is included. Here all 288 TN5 barcode sequences are merged into one set of fastq files.

## Index reads
For ScaleBio Methylation libraries the i7 and TN5 barcodes are included in read 1, while the i5 barcode is in read 2. Index reads are needed for demultiplexing and barcode correction. Hence we need to tell `bcl-convert` to generate index read fastqs using the `samplesheet.csv` setting: \
`CreateFastqForIndexReads,1`

# Using pre-generated fastq files as workflow input
Set `--fastqDir` to the directory containing the fastq files. 
There are two ways that you can split up the barcode demultiplexing and correction step.

The first option is split by lane, where file names should follow the pattern `<LibraryName>_<S#>_<Lane>_<Read>_...fastq.gz`, where
* `Name` is the library name (`ScaleMethyl` by default, can be set in the `libName` column in `samples.csv`)
* `Read` is one of `R1`, `R2`, `I1`, `I2`
* `Lane` is the lane number from using `--no-lane-splitting false` all lanes are combined for each sample following the demultiplexing step
An example [samplesheet.csv](examples/samplesheet.csv) for this default configuration included.

The second option is split by PCR indexes (i7, i5, or both), where file names should follow the pattern `<LibraryName>_<IndexString>_<S#>_<Read>_...fastq.gz`, where
* `Name` is the library name (`ScaleMethyl` by default, can be set in the `libName` column in `samples.csv`)
* `Read` is one of `R1`, `R2`, `I1`, `I2`
* `IndexString` is a letter string (ex. i5 sequence) that is used to split the sample fastq into smaller files for demultiplexing and each sample is later merged together. 
* No `Lane` splitting in this configuration `--no-lane-splitting true`
An example [samplesheet_pcr_split.csv](examples/samplesheet_pcr_split.csv) for this configuration included.

# Using the ScaleMethyl workflow to generate fastq files
The ScaleMethyl workflow also has the capability to generate fastq files when given a runFolder as input. The `makeBclConvertSamplesheet` process in main.nf calls a python script which generates a bcl-convert samplesheet from the samples csv provided as input and the RunInfo.xml file. The samples csv contains information about whether a subset or a full plate of TN5, i5 and i7 barcodes were used. Based on whether `splitFastq` is set to true or not, the bcl-convert samplesheet generated contains samples split by i5. The `makeBclConvertSamplesheet` process will be invoked by default if a runFolder is provided and a bcl-convert samplesheet is not.
