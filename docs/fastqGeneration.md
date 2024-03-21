# Fastq generation
The workflow can be started from a sequencer runFolder (bcl) [`--runFolder`]. In that case fastq files will be generated internally, using Illumina `bcl-convert`.

Alternatively, it is also possible to generate fastq files upstream, for example when the ScaleMethyl library is multiplexed together with other libraries during sequencing at a core facility. In that case we recommend using [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html), but older versions of the Illumina software are also possible.

Here is an example samplesheet for bcl-convert [samplesheet_pcr_split_smallkit.csv](examples/samplesheet_pcr_split_smallkit.csv) with typical options included. All 288 TN5 barcode sequences are merged into one set of fastq files.

## Index reads
For ScaleBio Methylation libraries the tagmentation barcodes are included in read 2, while the MET i5 and i7 PCR barcodes are in index reads 1 and 2. Since the index reads are needed for demultiplexing and barcode correction, we need to tell `bcl-convert` to generate index read fastqs using the `samplesheet.csv` setting: \
`CreateFastqForIndexReads,1`

# Using pre-generated fastq files as workflow input
Set `--fastqDir` to the directory containing the fastq files. 
There are two ways that you can split up the barcode demultiplexing and correction step.

The first default option is split by PCR indexes (i7, i5, or both), where file names should follow the pattern `<LibraryName>_<IndexString>_<S#>_<Read>_...fastq.gz`, where
* `LibraryName` is the name of the library, which needs to match the `libName` column in `samples.csv`
* `IndexString` is a letter string (ex. i5 sequence) that is used to split the sample fastq into smaller files for demultiplexing and each
  sample is later merged together
* `S#`is the character S followed by a number representing the total number of unique i5/i7 index combinations (small kit = 192; large kit = 768)
    * e.g. `S1` = well `1A`; `S96` = well `12H`
* `Read` is one of `R1`, `R2`, `I1`(i7), `I2`(i5)
    * `Read` is required to be in the 4th field (after the 3rd `_`) 
* No `Lane` splitting in this configuration `--no-lane-splitting true`  
    - See [samplesheet_pcr_split_largekit.csv](examples/samplesheet_pcr_split_largekit.csv) and [samplesheet_pcr_split_smallkit.csv](examples/samplesheet_pcr_split_smallkit.csv) for example bcl-convert samplesheets split by PCR.  

The second option is split by lane, where file names should follow the pattern `<LibraryName>_<S#>_<Lane>_<Read>_...fastq.gz`, where
* `LibraryName` is the name of the library, which needs to match the `libName` column in `samples.csv`
* `Read` is one of `R1`, `R2`, `I1`(i7), `I2`(i5)
    * `Read` is required to be in the 4th field (after the 3rd `_`) 
* `Lane` is the flow cell lane number
    * When using `--no-lane-splitting false` all lanes are combined for each sample following the demultiplexing step  
    - See [samplesheet_largekit.csv](examples/samplesheet_largekit.csv) and [samplesheet_smallkit.csv](examples/samplesheet_smallkit.csv) for example bcl-convert samplesheets split by lane.  

# Using the ScaleMethyl workflow to generate fastq files
The ScaleMethyl workflow also has the capability to generate fastq files when given a runFolder as input. The `makeBclConvertSamplesheet` process in main.nf calls a python script which generates a bcl-convert samplesheet from the samples csv provided as input and the RunInfo.xml file. The samples csv contains information about whether a subset or a full plate of TN5, i5 and i7 barcodes were used. Based on whether `splitFastq` is set to true or not, the bcl-convert samplesheet generated contains samples split by i5. The `makeBclConvertSamplesheet` process will be invoked by default if a runFolder is provided and a bcl-convert samplesheet is not.  

To generate the samplesheet.csv outside of the workflow, this can be done with `bclConvertSheet.py` using the [scaleMethyl environment](../scaleMethyl.conda.yml).  

```
bin/bclConvertSheet.py -h 
usage: bclConvertSheet.py [-h] [--splitFastq] SAMPLES.csv LIBRARY.json RUNINFO.xml

Create bcl_convert samplesheet.csv from workflow samples.csv

positional arguments:
  SAMPLES.csv   CSV with samples and index sequences for ScaleMethyl workflow run
  LIBRARY.json  Library structure definition
  RUNINFO.xml   Sequencer runinfo (in runfolder)

optional arguments:
  -h, --help    show this help message and exit
  --splitFastq  Split sample by PCR index barcode sequence
```
