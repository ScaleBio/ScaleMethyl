#!/usr/bin/env Rscript

args <- commandArgs()

help <- function(){
    cat("met_tssEnrich.R : calculate per cell tss enrichment ratio\n")
    cat("Usage: \n")
    cat("--tssBed        : TSS window bed file (TSS+/-100bp) tss and background must have the same #bin and bin size [required]\n")
    cat("--backgroundBed : Background window bed where tss has been shifted 1kb (ex TSS -1.2kb to -1kb)              [required]\n")    
    cat("--outName       : Output prefix (default is complexity file 1 prefix)                                       [required]\n")
    cat("--annoFile      : barcode annotation with pass_filter column (pass/fail)                                    [required]\n")
    cat("--bamFile       : deduplicated bam file                                                                     [required]\n")
    cat("\n")
    q()
}

if( !is.na(charmatch("-h",args)) || !is.na(charmatch("--help",args))){
    help()
} else {
    tssBed        <- sub('--tssBed=',       '', args[grep('--tssBed=', args)] )
    backgroundBed <- sub('--backgroundBed=','', args[grep('--backgroundBed=', args)] )    
    outName       <- sub('--outName=',      '', args[grep('--outName=', args)] )
    annoFile      <- sub('--annoFile=',     '', args[grep('--annoFile=', args)] )
    bamFile       <- sub('--bamFile=',      '', args[grep('--bamFile=', args)] )
}

library(GenomicAlignments)
library(GenomicRanges)
library(GenomicFeatures)
library(plyr)
library(rtracklayer)

# Load genomic regions
tss        <- import(tssBed)
background <- import(backgroundBed)
# Toss chromosomes x,y,mt
seqlevels(tss,pruning.mode="coarse") <- seqlevels(tss)[grep("chrY|chrM|hs_Y|hs_M|mm_Y|mm_M",seqlevels(tss), invert=TRUE)]
seqlevels(background,pruning.mode="coarse") <- seqlevels(background)[grep("chrY|chrM|hs_Y|hs_M|mm_Y|mm_M",seqlevels(background), invert=TRUE)]

# Get barcodes for each bam
anno <- read.table(annoFile, header=TRUE, sep=",")
anno <- anno[anno$pass_filter == "pass",]

# Load bam file
param  <- ScanBamParam(what='mapq')
bd     <- readGAlignments(bamFile, param=param, use.names=TRUE)

# Add the barcode from read_id
mcols(bd)$barcode <- sub(":.*$", "", names(bd))
# Get only reads associated with passing cells
bd                <- bd[mcols(bd)$barcode %in% anno$BC]
anno              <- anno[anno$BC %in% mcols(bd)$barcode,]

# Loop over each barcode
enrich_df <- do.call(rbind, lapply(unique(mcols(bd)$barcode), tssEnrich <- function(bc){
    # Filter the current barcode reads
    bd_sub             <- bd[mcols(bd)$barcode %in% bc]
    # Count overlapping reads with tss regions
    tss_counts         <- sum( countOverlaps(bd_sub, tss, ignore.strand = TRUE) >0 )
    #Count overlapping reads with the background regions
    background_counts  <- sum( countOverlaps(bd_sub, background, ignore.strand = TRUE) >0 )
    # Calculate tss enrichemnt for the current barcode
    tss_enrich         <- round( tss_counts/background_counts, 4 )
    # Export as a dataframe to row bind into the final dataframe
    data.frame(BC = bc, tss_enrich = tss_enrich, tss_counts = tss_counts, background_counts = background_counts)
}))

write.table(enrich_df, file=sub("$", ".tss_enrich.csv", outName), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

print("done")
