/**
* Module that deduplicates aligned bam files and extracts methylation stats
*
* Process:
*     Dedup
*     Extract
*/

// Deduplicate reads in aligned bam files using sc_dedup
process Dedup {
input: 
	tuple val(sample), path(bam)
	path(chroms)
output: 
	tuple val(sample), path("${sample}.dedup.bam"), emit: bam
	tuple val(sample), path("${sample}.cell_stats.tsv"), emit: complexity
	tuple val(sample), path("${sample}.dedup_stats.tsv"), emit: stats
	tuple val(sample), path("${sample}.fragment_hist.tsv"), emit: fragHist
	tuple val(sample), env(TOTAL_READS), emit: totalReads
publishDir { outDir }, pattern: "${sample}.dedup.bam", enabled: params.bamDedupOut
tag "$sample"
label 'process_medium'
script:
    outDir = file(params.outDir) / "alignments" / "dedup" / sample
	sthreads = task.cpus
"""
	samtools sort --threads $sthreads $bam -o ${sample}.coodsort.bam
	sc_dedup ${sample}.coodsort.bam --duplicate-key ${params.dedupKey} --min-mapq ${params.minMapq} --barcode-input Qname --out-prefix ${sample} --write-threads $sthreads --genome $chroms
    # count number of passing reads
    TOTAL_READS=\$(cat ${sample}.cell_stats.tsv | awk 'BEGIN { total_reads = 0 } NR > 1 { total_reads += \$3 } END { print total_reads }')
"""
}

// Extract methylation stats from deduplicated bam files
process Extract {
input: 
	tuple val(sample), path(bam), val(totalReads)
output: 
	tuple val(sample), path("${sample}.met_CG.parquet"), emit: metCG
	tuple val(sample), path("${sample}.met_CH.parquet"), emit: metCH
	tuple val(sample), path("${sample}.cellInfo.txt"), emit: report
tag "$sample"
script:
    nprocs = Math.max(task.cpus - 1, 1)
"""
	met_extract.py $bam --sample $sample --threshold ${params.chReadsThreshold / 100} --subprocesses $nprocs
"""
}

workflow DEDUP_AND_EXTRACT {
take:
    dedupBamInput // Aligned bam files from bsbolt
    genome // Reference genome

main:
    // Run scDedup for removing duplicate reads from bam files produced by bsbolt
    Dedup(dedupBamInput, genome.bsbolt_chrs)
    
    // Extract methylation stats from deduplicated bam files
    Extract(
        Dedup.out.bam
            .join(Dedup.out.totalReads)
            .filter { sample, bam, totalReads -> totalReads as int > 0 }
    )
    
    // Collect all cell_stats.tsv files for a sample
    complexTabs = Dedup.out.complexity.map({it[1]}).map { file ->
        def ns = file.getName().toString().tokenize('.').get(0)
        return tuple(ns, file)
    }.groupTuple()
    // Collect all fragment_hist.tsv files for a sample
    fragTabs = Dedup.out.fragHist.map({it[1]}).map { file ->
        def ns = file.getName().toString().tokenize('.').get(0)
        return tuple(ns, file)
    }.groupTuple()
    readsFrags = complexTabs.join(fragTabs)
    // readFrags -> [sample name, [list of cell_stats.tsv files], [list of fragment_hist.tsv files]]
    readsFrags.dump(tag: 'readFrags')
    sampleExtractStats = Extract.out.report.map({it[1]}).map { file ->
        def ns = file.getName().toString().tokenize('.')
        return tuple(ns.get(0), file)
    }.groupTuple()
    // sampleExtractStats -> [sample name, [list of cellInfo.txt files]]
    sampleExtractStats.dump(tag: 'sampleExtractStats')
    dedupStats = Dedup.out.stats.map({it[1]}).map { file ->
            def ns = file.getName().toString().tokenize('.')
            return tuple(ns.get(0), file)
    }.groupTuple()
    dedupStats.dump(tag: 'dedupStats')

emit:
    readsFrags = readsFrags // Cell stats and fragment histogram files from scDedup
    sampleExtractStats = sampleExtractStats // cellInfo files from Extract
    dedupStats = dedupStats // scDedup stats
    dedupBam = Dedup.out.bam // Deduplicated bam files
    metCG = Extract.out.metCG // CG methylation stats
    metCH = Extract.out.metCH // CH methylation stats
}
