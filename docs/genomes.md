# Reference genomes
The workflow requires a genome reference and annotation to run. All files and settings for a genome are defined in a [genome.json](examples/genomes.json) file. When launching the workflow, the reference genome is selected by passing the path to a specific `genome.json` in the `--genome` parameter.

The `genome.json` file includes

Field |  Description | Required? | Example
:-- | -- | -- | --
name | The name of the species / genome-version | Required | human 
bsbolt_index | Path to the BSBolt index directory | Required (one index required) | `/PATH/TO/bsbolt.ref` 
bwa_index | Path to the bwa-meth (bwa-mem2) index directory | Required (one index required) | `/PATH/TO/bwa-meth/bwa-mem2.ref` 
parabricks_index | Path to the parabricks (bwa-meth/bwa-mem) index directory | Required (one index required) | `/PATH/TO/bwa-meth/bwa-mem.ref`
ref_fasta | Name of fasta file within bwa-meth/parabricks index | Required (with bwa_index or parabricks_index) | `hg38.fa`
genomeTiles | Path to binned genome sorted bed file for CG matrix | Optional | `/PATH/TO/50kbp.bed` 
genomeTilesCh | Path to binned genome sorted bed file for CH matrix | Optional | `/PATH/TO/100kbp.bed` 
filter_chrs | Path to tsv chromosome labels (mito, filter) to filter from deduplicated BAM | Required | `/PATH/TO/filter_chrs.tsv` 
tssWin | Path to bed sorted 200nt windows centered at TSS | Required unless runTssEnrich: false | `/PATH/TO/tss.bed` 
backgroundWin | Path to bed sorted 200nt windows centered at TSS -1kb upstream | Required unless runTssEnrich: false | `/PATH/TO/background.bed` 

* All files (`bsbolt_index`, `bins`, ...) can be specified either as
    - an absolute path (`/path/to/genome`)
    - a relative path starting from the location of the `genome.json` file (`genes/bins.bed`)
    - a AWS S3 url (s3://path/to/genome)

## bwa-meth index (default)
The bwa-meth index needs to be built with bwa-meth version `>= 0.2.7` and bwa-mem2 version `>= 2.2.1`. See the bwa-meth [documentation](https://github.com/brentp/bwa-meth) for additional options. An example command would be
```
bwameth.py index-mem2 {fasta reference}
```

## parabricks index (GPU)
The parabricks index needs to be built with bwa-meth version `>= 0.2.7` and bwa-mem version `>= 0.7.19`. See the bwa-meth [documentation](https://github.com/brentp/bwa-meth) for additional options. An example command would be
```
bwameth.py index {fasta reference}
```

## BSBolt index (legacy)
The BSBolt index needs to be built with BSBolt version `>= 1.5.0`. See the BSbolt [documentation](https://bsbolt.readthedocs.io/en/latest/bsb_index/) for additional options. An example command would be
```
bsbolt Index -G {fasta reference} -DB {database output}
```

## Annotation
All genomic non-overlapping bins in the provided _genomeTiles_ and _genomeTilesCh_ bed files are used as features for CG and CH methylation matrix generation. 
* To use bed files with different sized bins (default 50kb) for CG matrix, pass to _genomeTiles_ in the `genome.json`. 
* To pass a different size bins (default 250kb) for CH matrix, pass to _genomeTilesCh_ in the `genome.json`.

## Pre-built genomes
Pre-built reference genomes are available for download:
* Human: 
    - bwa-meth (default): https://s3.us-east-2.amazonaws.com/scale.pub/genomes/methyl/grch38-bwa.tgz
    - parabricks (GPU): https://s3.us-east-2.amazonaws.com/scale.pub/genomes/methyl/grch38-parabricks.tgz
    - bsbolt (legacy): https://s3.us-east-2.amazonaws.com/scale.pub/genomes/methyl/grch38-bsbolt.tgz
* Mouse: 
    - bwa-meth (default): https://s3.us-east-2.amazonaws.com/scale.pub/genomes/methyl/mm39-bwa.tgz
    - parabricks (GPU): https://s3.us-east-2.amazonaws.com/scale.pub/genomes/methyl/mm39-parabricks.tgz
    - bsbolt (legacy): https://s3.us-east-2.amazonaws.com/scale.pub/genomes/methyl/mm39-bsbolt.tgz
* Human/Mouse Barnyard: 
    - bwa-meth (default): https://s3.us-east-2.amazonaws.com/scale.pub/genomes/methyl/grch38_mm39-bwa.tgz
    - parabricks (GPU): https://s3.us-east-2.amazonaws.com/scale.pub/genomes/methyl/grch38_mm39-parabricks.tgz
    - bsbolt (legacy): https://s3.us-east-2.amazonaws.com/scale.pub/genomes/methyl/grch38_mm39-bsbolt.tgz

Download these to your analysis server, unpack them and then use e.g.
`--genome /PATH/TO/unpacked/folder/genome.json`

