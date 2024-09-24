# Reference genomes
The workflow requires a genome reference and annotation to run. All files and settings for a genome are defined in a [genome.json](examples/genomes.json) file. When launching the workflow, the reference genome is selected by passing the path to a specific `genome.json` in the `--genome` parameter.

The `genome.json` file includes

Field |  Description | Required? | Example
:-- | -- | -- | --
name | The name of the species / genome-version | Required | human 
bsbolt_index | Path to the BSBolt index directory | Required | `/PATH/TO/bsbolt.ref` 
genomeTiles | Path to binned genome sorted bed file for CG matrix | Required | `/PATH/TO/50kbp.bed` 
genomeTilesCh | Path to binned genome sorted bed file for CH matrix | Required | `/PATH/TO/100kbp.bed` 
bsbolt_chrs | Path to tsv chromosome labels (mito, filter) to filter from deduplicated BAM | Required | `/PATH/TO/bsbolt_chrs.tsv` 
tssWin | Path to bed sorted 200nt windows centered at TSS | Required | `/PATH/TO/tss.bed` 
backgroundWin | Path to bed sorted 200nt windows centered at TSS -1kb upstream | Required | `/PATH/TO/background.bed` 

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
All genomic non-overlapping bins in the provided _genomeTiles_ and _genomeTilesCh_ bed files are used as features for CG and CH methylation matrix generation. 
* To use bed files with different sized bins (default 50kb) for CG matrix, pass to _genomeTiles_ in the `genome.json`. 
* To pass a different size bins (default 250kb) for CH matrix, pass to _genomeTilesCh_ in the `genome.json`.

## Pre-built genomes
Pre-build reference genome for human is available for download:
* Human: http://scale.pub.s3.amazonaws.com/genomes/methyl/grch38.tgz
* Mouse: http://scale.pub.s3.amazonaws.com/genomes/methyl/mm39.tgz
* Human/Mouse Barnyard: http://scale.pub.s3.amazonaws.com/genomes/methyl/grch38_mm39.tgz

Download these to your analysis server, unpack them and then use e.g.
`--genome /PATH/TO/grch38/grch38.json`

