/**
* Module that merges bam files in two bsbolt output directories
*
* Process:
*     MergeBam
*/

// Merge two bam files together using samtools
process MergeBam {
input: 
	tuple val(sample), val(tgmt), path("bam1.${tgmt}.bam"), path("bam2.${tgmt}.bam")
output: 
	tuple val("${sample}.${tgmt}"), path("${sample}.${tgmt}.bam"), emit: bam
publishDir "${params.outDir}/bamMerge/${sample}.${tgmt}", pattern: "${sample}.${tgmt}.bam", enabled: params.bamMergeOut
tag "$sample"
script:
	sthreads = task.cpus
"""
	samtools addreplacerg -r ID:bam1 -r LB:bam1 bam1.${tgmt}.bam | samtools sort -@ $sthreads -o bam1.bam -
	samtools addreplacerg -r ID:bam2 -r LB:bam2 bam2.${tgmt}.bam | samtools sort -@ $sthreads -o bam2.bam -
	samtools merge --threads $sthreads -o ${sample}.${tgmt}.bam bam1.bam bam2.bam
"""
}

workflow INPUT_BAM_READS{
take:
    bam1Dir // Path to the first bsbolt output directory
    bam2Dir // Path to the second bsbolt output directory

main:
    // The bam input directory structure is expected to be exactly the same as the
    // output directory structure created by bsbolt in the align process
    bam1 = Channel.fromPath("$bam1Dir/*/*bam", checkIfExists: true)
    bam2 = Channel.fromPath("$bam2Dir/*/*bam", checkIfExists: true)
    b1 = bam1.map { file ->
        def ns = file.getBaseName().toString().tokenize('.')
        return tuple(ns[0], ns[1], file)
    }.groupTuple(by:[0,1])
    b2 = bam2.map { file ->
        def ns = file.getBaseName().toString().tokenize('.')
        return tuple(ns[0], ns[1], file)
    }.groupTuple(by:[0,1])
    // b1 -> [sample name, well coordinate, bam file for that sample name and well coordinate]
    b1.dump(tag: 'b1')
    // b2 -> [sample name, well coordinate, bam file for that sample name and well coordinate]
    b2.dump(tag: 'b2')
    bams = b1.join(b2, by: [0,1], failOnDuplicate: true)
    bams.dump(tag: 'bams')
    MergeBam(bams)
    // MergeBam.out.bam -> [sample name.well coordinate, sample name.well coordinate.bam]
    MergeBam.out.bam.dump(tag: 'MergeBam.out.bam')

emit:
    mergedBam = MergeBam.out.bam // Per TN5 merged BAM file
}
