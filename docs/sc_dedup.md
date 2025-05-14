# sc_dedup
`sc_dedup` is a ScaleBio developed tool for removing duplicate reads from an aligned BAM file. The tool is barcode-aware. Reads from different barcodes aren't considered duplicates of one another.

As used in scaleMethyl, sc_dedup uses the leftmost position of the leftmost aligned fragment in a mate-pair and compares it to the corresponding position of previously
encountered reads in the <IN.BAM> file. For singletons and mate-pairs with one unmapped fragment, the leftmost alignment position of the mapped fragment will be
used.  The rightmost fragment in a mate-pair is not considered. After the first read at a position is encountered, future reads whose alignments 
start at that position will be discarded as duplicates. See [nextflow.config](../nextflow.config) for parameter settings used by the dedup process in the workflow. 


## Executing sc_dedup
For a complete list of flags and options sc_dedup recognizes, please run `sc_dedup --help`

## Inputs
* <IN.BAM> /path/to/bamFile Where bamFile contains read alignments (coordinate sorted) 
    
## Outputs
When run with the `bamDedupOut` option in nextflow.config set to true, sc_Dedup will publish 288 folders under the `bamDeDup` folder in the output directory.  Each folder will contain one bam file (one for each TN5 barcode).
Output filenames will begin using the <IN.BAM> filename (without the `.bam`).
* <filename>.dedup.bam               The file will contain the input reads that weren't discared as duplicates.

### Metrics
The output directory will contain the following metrics files
  *     <Filename>.cell_stats.tsv    The read metrics for each cell barcode in the <IN.BAM> file.
  *     <filename>.dedup_stats.tsv   The summary read statistics for all reads in the
            <filename.dedup.bam> file
  *     <filename>.fragment_hist.tsv Readlength histogram of all reads in the 
            <filename.dedup.bam> file
    
### Options
The `--genome <genome.tsv>` option specifies the reference genome chromosome information for metrics. In order to get accurate counts of total and unique reads, it is important to include the line `chrM    mito` in the <genome.tsv> specified by this option.  An example genome.tsv can be found at the location referred to by the filter_chrs setting in [genomes.json](examples/genome.json)

