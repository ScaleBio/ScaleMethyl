/**
* Module that performs alignment of sequencing reads to a reference genome.
* Process:
*     AlignBsBolt
*     AlignBWA
*     AlignParabricks
*/

process AlignBsBolt {
input: 
	path(index) // index directory
	tuple val(sample), path(pairs)
output: 
	tuple val(sample), path("${sample}.bam"), emit: bam
	tuple val(sample), path("${sample}.*.log"), emit: log
publishDir { outDir }, pattern: "${sample}.bam", enabled: params.bamOut
tag "$sample"
script:
    outDir = file(params.outDir) / "alignments" / "bam" / sample
	athreads = task.cpus - 2 
	othreads = 2
"""
    bsbolt Align -F1 ${pairs[1]} -F2 ${pairs[0]} -t $athreads -OT $othreads -O $sample -DB $index >> ${sample}.bsbolt.log 2>> ${sample}.bsbolt.log
"""
}

process AlignBWAMeth {
input: 
	path(index) // index directory
	tuple val(sample), path(pairs)
    val(fastaName) // Name of the fasta file to be used for bwa-meth
output: 
	tuple val(sample), path("${sample}.bam"), emit: bam
	tuple val(sample), path("${sample}.*.log"), emit: log
publishDir { outDir }, pattern: "${sample}.bam", enabled: params.bamOut
tag "$sample"
script:
    outDir = file(params.outDir) / "alignments" / "bam" / sample
	athreads = task.cpus - 2 
	othreads = 2
"""
    # Touch the index files to update the timestamp for bwa-meth
    touch -d \"\$(date -R -r $index/${fastaName}.bwameth.c2t) - 2 hours\" -m $index/${fastaName}.bwameth.c2t
    bwameth.py --threads $athreads --reference $index/$fastaName ${pairs[1]} ${pairs[0]} 2>> ${sample}.bwa_meth.runtime | samtools view -@ $othreads -bS - -o ${sample}.bam
    samtools flagstat -O tsv -@ $task.cpus ${sample}.bam > ${sample}.bwa_meth.log
"""
}

process AlignParabricks {
input: 
	path(index) // index directory
	path(fastqR1) // List of fastq files
    path(fastqR2)
    val(sampleNames)
    val(fastaName) // Name of the fasta file to be used for parabricks
output: 
	path("outs/*.bam"), emit: bam
	path("outs/*.log"), emit: log
publishDir { outDir }, pattern: "outs/*.bam", enabled: params.bamOut, saveAs: { bamFile -> 
    def bamFileName = bamFile.tokenize('/').get(1)
    def sample=bamFileName.tokenize('.').get(0)+"."+bamFileName.tokenize('.').get(1)
    return "${sample}/${sample}.bam"
}

script:
    outDir = file(params.outDir) / "alignments" / "bam" 
    def fixedSampleNames = sampleNames.join(',')
    def sampleNumber = sampleNames.size()
"""
    mkdir -p outs
    # Build fastq list
    for i in `seq 1 $sampleNumber`; do
        sample=`echo $fixedSampleNames | cut -d',' -f\${i}`
        echo "\${sample}_R2.fq.gz \${sample}_R1.fq.gz '@RG\\tID:\${sample}\\tPU:\${sample}\\tSM:sample\\tLB:lib'" >> fastq_list.txt
        mkdir -p outs/\${sample}
    done
    # Touch the index files to update the timestamp for parabricks
    touch -d \"\$(date -R -r $index/${fastaName}.bwameth.c2t) - 2 hours\" -m $index/${fastaName}.bwameth.c2t
    pbrun fq2bam_meth --ref $index/$fastaName --in-fq-list fastq_list.txt  --out-bam full_out.bam --bwa-nstreams 1 --bwa-cpu-thread-pool $task.cpus --align-only --num-gpus 1
    
    # split output bam into sample bams
    samtools split -f "outs/%!.bam" -@ $task.cpus full_out.bam 

    # build stats file for each sample bam and move sample bams into proper output folders
    for i in `ls outs/*.bam`; do
        sample=`basename \$i`
        samtools flagstat -O tsv -@ $task.cpus \$i > outs/\${sample%%.bam}.bwa_meth.log
    done
"""
}

workflow ALIGNMENT {
take:
    trimFastq // from input_reads/Trim
    genome // reference genome
main: // Run bsbolt(aligner) on trimmed fastq files
    if(params.aligner == "bsbolt") {
        AlignOut=AlignBsBolt(genome.bsbolt_index, trimFastq)
        AlignOutBams = AlignOut.bam
        AlignOutLogs = AlignOut.log
    } else if (params.aligner == "bwa-meth") {
        AlignOut=AlignBWAMeth(genome.bwa_index, trimFastq, genome.ref_fasta)
        AlignOutBams = AlignOut.bam
        AlignOutLogs = AlignOut.log
    } else if (params.aligner == "parabricks") {
        // collate the samples into groups
        sampleChunks = trimFastq.collect(flat:false)
                                    .map{ items -> 
                                    def chunkSize = Math.ceil(items.size() / params.parabricksNumGpu) as int
                                    items.collate(chunkSize) }
                                    .flatMap().dump(tag:'sampleChunks')
        sampleNames = sampleChunks.map{ it.collect{it[0]} }
        fastqNamesR1 = sampleChunks.map{ it.collect{it[1][0]} }
        fastqNamesR2 = sampleChunks.map{ it.collect{it[1][1]} }
        AlignParabricks(genome.parabricks_index, fastqNamesR1, fastqNamesR2, sampleNames, genome.ref_fasta)
        AlignOutBams = AlignParabricks.out.bam.flatten().map{file -> tuple(file.getName().toString().tokenize('.').get(0)+"."+file.getName().toString().tokenize('.').get(1), file)}
        AlignOutLogs = AlignParabricks.out.log.flatten().map{file -> tuple(file.getName().toString().tokenize('.').get(0)+"."+file.getName().toString().tokenize('.').get(1), file)}
        
    } else {
        ParamLogger.throwError("Invalid aligner specified. Please choose from bsbolt, bwa-meth or parabricks")
    }

emit:
    alignedBam = AlignOutBams // Post alinger aligned bams
    alignLog = AlignOutLogs // aligner log file
}