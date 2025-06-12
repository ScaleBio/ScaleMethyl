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
publishDir { outDir }, pattern: "${sample}.dedup_stats.tsv"
publishDir { outDir }, pattern: "${sample}.cell_stats.tsv"
publishDir { outDir }, pattern: "${sample}.fragment_hist.tsv"
tag "$sample"
label 'process_medium'
script:
    outDir = file(params.outDir) / "alignments" / "dedup" / sample
	sthreads = task.cpus
"""
	samtools sort --threads $sthreads $bam -O bam -o ${sample}.coodsort.bam
	sc_dedup ${sample}.coodsort.bam --duplicate-key ${params.dedupKey} --min-mapq ${params.minMapq} --barcode-input Qname --out-prefix ${sample} --write-threads $sthreads --genome $chroms
    # count number of passing reads
    TOTAL_READS=\$(cat ${sample}.cell_stats.tsv | awk 'BEGIN { total_reads = 0 } NR > 1 { total_reads += \$3 } END { print total_reads }')
"""
}

// Extract methylation stats from deduplicated bam files
process Extract {
    
input: 
	tuple val(sample), path(bam), val(totalReads)
    val(aligner)
    path(index)
    val(fastaFile)
output: 
	tuple val(sample), path("${sample}.met_CG.parquet"), emit: metCG
	tuple val(sample), path("${sample}.met_CH.parquet"), emit: metCH, optional: true
	tuple val(sample), path("${sample}.cellInfo.txt"), emit: report
tag "$sample"
publishDir file(params.outDir) / "matrix_generation", pattern: "${sample}.met_*.parquet", mode: 'copy', enabled: params.parquetOut
publishDir {outDir}, pattern: "${sample}.cellInfo.txt", mode: 'copy'
script:
    outDir = file(params.outDir) / "alignments" / "dedup" / sample
    nprocs = Math.max(task.cpus - 1, 1)
    contexts = (params.calculateCH) ? "CG,CH" : "CG"
    options = (aligner == 'bwa-meth' || aligner == 'parabricks') ? "--ref $index/$fastaFile" :  ""
"""
	met_extract.py $bam --sample $sample --threshold ${params.chReadsThreshold / 100} --subprocesses $nprocs --contexts $contexts --aligner $aligner $options
"""
}

workflow DEDUP_AND_EXTRACT {
take:
    dedupBamInput // Aligned bam files from aligner
    genome // Reference genome

main:
    // Run scDedup for removing duplicate reads from bam files produced by aligner
    Dedup(dedupBamInput, genome.filter_chrs)
    
    // Extract methylation stats from deduplicated bam files
    if(params.aligner == 'bwa-meth') {
        index = genome.bwa_index
        fastaFile = genome.ref_fasta
    } else if(params.aligner == 'parabricks') {
        index = genome.parabricks_index
        fastaFile = genome.ref_fasta
    } else {
        index = []
        fastaFile = ""
    }
    Extract(
        Dedup.out.bam
            .join(Dedup.out.totalReads)
            .filter { sample, bam, totalReads -> totalReads as int > 0 },
        params.aligner,
        index,fastaFile
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
