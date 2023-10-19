# Reference genomes
The workflow requires a genome reference and annotation to run. All files and settings for a genome are defined in a [genome.json](examples/genomes.json) file. When launching the workflow, the reference genome is selected by passing the path to a specific `genome.json` in the `--genome` parameter.

The `genome.json` file includes

Field |  Description | Required? | Example
:-- | -- | -- | --
name | The name of the species / genome-version | Required | human 
bsbolt_index | Path to the BSBolt index directory | Required | `/PATH/TO/bsbolt.ref` 
genomeTiles | Path to binned genome sorted bed file | Required | `/PATH/TO/mainChr.50kbp.bed` 

* All files (`bsbolt_index`, `bins`, ...) can be specified either as
    - an absolute path (`/path/to/genome`)
    - a relative path starting from the location of the `genome.json` file (`genes/bins.bed`)
    - a AWS S3 url (s3://path/to/genome)

## BSBolt index
The provided BSBolt index needs to be built with BSBolt version `>= 1.5.0`. See the BSbolt [documentation](https://bsbolt.readthedocs.io/en/latest/bsb_index/) for additional options. An example command would be
```
bsbolt Index -G {fasta reference} -DB {database output}
```

## Annotation
All annotated transcripts in the provided _GTF_ file are included in the analysis (i.e. in the output methylation matrices). 
Genomes are binned into 50kb bins for that are used to create binned methylation matrices. To use bed files with different binning, pass to _genomeTiles_ in the `genome.json`

## Pre-built genomes
Pre-build reference genome for human is available for download:
* Human: http://scale.pub.s3.amazonaws.com/genomes/grch38.tgz
* Mouse: http://scale.pub.s3.amazonaws.com/genomes/mm39.tgz
* Human/Mouse Barnyard: http://scale.pub.s3.amazonaws.com/genomes/grch38_mm39.tgz

Download these to your analysis server, unpack them and then use e.g.
`--genome /PATH/TO/grch38/grch38.json`

