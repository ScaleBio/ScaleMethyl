/**
* Module that performs runfolder to fastq conversion, sample demux, barcode error correction, trimming,
* alignment and qc of demuxed fastq files
* 
* Process:
*     MakeBclConvertSamplesheet
*     BclConvert
*     Fastqc
*     BarcodeDemux
*     MergeDemux
*     Trim
*/

def constructMatchingKey(fname) {
	def identifyingCharRemoved = fname.replaceAll(/_[RI][12][_\.]/, "")
	return identifyingCharRemoved
}

// Extract sampleName (and optionally subsample, split by TN5 well) for one bcParser output
// @bcParserOut -> library name, fastq file
def constructSampleName (bcParserOut, splitFastq) {
	def libName = bcParserOut[0]
	def file = bcParserOut[1]
    def sampleName = ""
    def subsampleName = ""
    // If splitFastq is true, sample name is the first part of the filename before the first '.'
    // subsampleName is the first part of the filename before the first '_' and after the first '.'
        // subsampleName indicates unique identifier which corresponds to the fastq-num we provide as input to bcParser
    if (splitFastq) {
        def tok = file.getName().toString().tokenize('_')
        sampleName = tok[0]
        subsampleName = "${sampleName}.${(tok[1].tokenize('_'))[0]}"
    }
    // If splitFastq is false, sample name is the first part of the filename before the first '_'
    // There is not a well coordinate in the filename when splitFastq is false so we tokenize on '_' to get the sample name
    else {
        sampleName = file.getName().toString().tokenize('_')[0]
        subsampleName = sampleName
    }
	return tuple(libName, sampleName, subsampleName, file)
}

// Create a bcl-convert samplesheet for libraries in samples.json
process MakeBclConvertSamplesheet {
input: 
	path(samplesCsv)
	path(libStructDir) // Directory containing the library structure definition (barcode sequence lists are loaded from here)
	val(libStructName) // Filename of the library structure definition .json
	path(runinfo) // RunInfo.xml from Sequencer RunFolder (Read-lengths etc.)
output: 
	path("samplesheet.csv")
publishDir file(params.outDir) / "fastq", mode: 'copy'
label 'process_single'
script:
	opts = ""
	libJson = "$libStructDir/$libStructName"
	if (params.splitFastq) {
		opts = "--splitFastq"
	}
"""
	bcl_convert_sheet.py $samplesCsv $libJson $runinfo $opts > samplesheet.csv
"""
}

// Run bcl-convert, used when starting from a sequencing run-folder
// Requires a separate bcl-convert samplesheet
process BclConvert {
input: 
	path(run)
	path(samplesheet)
	val(bclConvertParams)
output: 
	path("fastq/*fastq.gz"), emit: fastq
	path("fastq/Reports/*"), emit: stats
publishDir params.outDir, pattern: 'fastq/Reports/*', mode: 'copy'
publishDir params.outDir, pattern: 'fastq/*.fastq.gz', enabled: params.fastqOut
script:
"""
	bcl-convert --sample-sheet $samplesheet --bcl-input-directory $run \
    --bcl-num-conversion-threads 4 --bcl-num-compression-threads 4 --bcl-num-decompression-threads 4 \
    $bclConvertParams --output-directory fastq
"""
}

// Run fastqc on either input fastq files or bcl-convert generated fastq files
process Fastqc {
input:
	path(fqFiles)
	path(adapters)
output:
	path("fastqc/*.html")
publishDir file(params.outDir) / "fastq", mode: 'copy'
label 'process_single'
script:
"""
	mkdir -p fastqc
	fastqc --adapters $adapters -o fastqc $fqFiles --threads ${task.cpus} -f fastq
"""
}

// Demux fastq files and do barcode error correction
process BarcodeDemux { 
input: 
	path(sheet) // samples.csv 
	path(libStructDir) // Directory containing the library type definition file (barcode sequence lists are loaded from here) 
	val(libStructJson) // Filename of the library definition .json 
	tuple val(libName), val(matchingKey), path(fqFiles), val(fastqNum) // Input fastq file 
output: 
	tuple val(libName), path("$outDir/*_S[1-9]*_R1*.fastq.gz"), emit: read1Fastq
	tuple val(libName), path("$outDir/*_S[1-9]*_R2*.fastq.gz"), emit: read2Fastq 
	path("$outDir/*_S0_*.fastq.gz"), emit: unknown optional true 
	path("$outDir/*.tsv") 
	tuple val(libName), path("$outDir/metrics.json"), emit: metrics
publishDir file(params.outDir) / "barcodes", pattern: "$outDir/*{txt,tsv,json}", mode: 'copy'
publishDir file(params.outDir) / "barcodes", pattern: "$outDir/*gz", enabled: params.fastqOut
/* Optionally publish barcode demuxed fastq files */
tag "$libName" 
label 'process_medium'
script: 
	if (params.splitFastq) {
		outDir = "${libName}.${fastqNum}.demux" 
	} else {
		outDir = "${libName}.demux"
	} 
	libStruct = "$libStructDir/$libStructJson" 
""" 
	bc_parser -v --write-fastq --demux $sheet --lib-name $libName --reads ${fqFiles.join(" ")} --lib-struct $libStruct --out $outDir --fastq-num $fastqNum --write_barcode_sequence
""" 
}

// Merge bc_parser stats from multiple jobs for a single library
process MergeDemux {
input:
    tuple val(libName), path("demux_metrics*.json") // bc_parser metrics from all bc_parser runs for a library
    path(libStructDir) // Directory containing the library type definition file (barcode sequence lists are loaded from here) 
    val(libStructJson) // Filename of the library definition .json 
output:
    tuple val(libName), path("${libName}.metrics.json"), emit: barcodeMetrics
publishDir file(params.outDir) / "barcodes", mode:'copy'
tag "$libName"
label 'process_single'
script:
	libStruct = "$libStructDir/$libStructJson" 
    """
    merge_bc_parser_output.py --bc_jsons demux_metrics* --lib_json $libStruct --libName $libName
    """	
}

// Run cutadapt on fastq files demultiplexed by bcParser
// One task gets all fastq files for a given sample / TN5 well, which are combined here
process Trim {
input:
	tuple val(sample), path(pairs1), path(pairs2)
output:
	tuple val(sample), path("${sample}_{R1,R2}.fq.gz"), emit: fastq
	tuple val(sample), path("*.trim_stats.json"), emit: stats
	path("*.trim_log"), emit: log
publishDir { outDir }, pattern: "*.fq.gz", enabled: params.trimOut
publishDir { outDir }, pattern: "*.json"
tag "$sample"
script:
    tthreads = Math.max(task.cpus - 2, 1)
    outDir = file(params.outDir) / "fastq" / "trim" / sample
"""
	cutadapt -j $tthreads -a CTATCTCTTATA -A AGATCGGAAGAGC -U -10 -m20 -o ${sample}_R1.fq.gz -p ${sample}_R2.fq.gz <(zcat cat ${pairs1}) <(zcat ${pairs2}) --json ${sample}.trim_stats.json | tee ${sample}.trim_log
"""
}

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

workflow INPUT_READS {
take:
    samples // sample and libName combination
    samplesCsv // samples.csv file
    libJson // Library structure definition file
    fastqDir // Directory containing fastq files
    runFolder // Directory containing sequencing run data

main:
    // Combination of libName and samples in a channel
	libNameAndSamples = samples.map{ [it.libName, it.sample] }
    // If run folder is provided, bcl-convert samplesheet will be generated(if not provided) and
    // bcl-convert will be executed to produce fastq files
    if (runFolder != null) {
        if (params.fastqSamplesheet == null) {
            // Generate bcl-convert samplesheet
            MakeBclConvertSamplesheet(samplesCsv, libJson.getParent(), libJson.getName(), 
                                      file("$runFolder/RunInfo.xml", checkIfExists:true))
            fqSheet = MakeBclConvertSamplesheet.out
        }
        else {
            fqSheet = file(params.fastqSamplesheet)
        }
        // Execute bcl-convert to generate fastq files
        BclConvert(file(runFolder), fqSheet, params.bclConvertParams)
        fqs = BclConvert.out.fastq.flatten().filter({ !it.name.contains("Undetermined") })
    } else if (params.fastqDir != null) {
        fqs = Channel.fromPath("$params.fastqDir/*fastq.gz", checkIfExists: true).filter({ !it.name.contains("Undetermined") })
    } else {
        ParamLogger.throwError("Must specify either 'runFolder' or 'fastqDir'")
    }
    // fqs -> [single fastq file]
    // n channels for n fastq files
    fqs.dump(tag:'fqs')
    if (params.fastqc) {
        Fastqc(fqs, params.adapters)
    }
    // Organize fastq files by sample
    // fqFiles -> (libName, matchingKey, file)
    fqFiles = fqs.map { file ->
        def fname = file.getName().toString()
        def libName = fname.tokenize('_')[0]
        def matchingKey = ""
        // Construct matching key to join on later on in the workflow
        if (params.splitFastq) {
            matchingKey = constructMatchingKey(fname)
        } else { // Matching key is just the library name if params.splitFastq is false
            matchingKey = libName
        }
        tuple(libName, matchingKey, file)
    }
    fqFiles.dump(tag:'fqFiles')

    pairedFqFiles = fqFiles.groupTuple(by:[0,1])
    pairedFqFiles.dump(tag:'pairedFqFiles')

    // This keeps only those sets of fastq files that correspond to a sample from our samples.csv
    // it.libName takes only the library name from samples csv
    // unique because samples csv has multiple entries for multiple samples, but only one library name
    // cross ensures that the right library names are paired with the right fastq files
    // The last map is to remove the extra library name element(it[0]) in the tuple that we get after the cross,
        // and add a unique counter to the tuple which will be used by bcParser for output file naming
    def nObjects = 0
    fqSamples = samples.map({it.libName}).unique().cross(pairedFqFiles)
    fqSamples = fqSamples.map{it[1] << nObjects++}
    fqSamples.dump(tag:'fqSamples')

    // Check that for each library (unique libName in samples.csv), we have a complete set of fastq files
    // checkFastq -> (libName, sampleName, Fastqs, counter)
    // This slightly strange channel construction (including it.id) works around the lack of 'left-join'
        // in nextflow (this way we can distinguish 'left-only' and 'right-only' join results)
    checkFastq = libNameAndSamples.unique{ it[0] }.join(fqSamples, remainder: true)
    checkFastq.dump(tag:'checkFastq')
    checkFastq.map {
        // If samples.csv has a library name which does not have corresponding fastq files
        // then 3rd element of checkFastq will be null
        if (it[3] == null) {
            ParamLogger.throwError("Library ${it[0]} does not have matching fastq files. None of the provided fastq filenames start with \"${it[0]}_\"")
        }
        // Check that each library has index1, read1 and read2 fastq files
        if (it[3].any {fq -> fq.getName().contains("_R1_") or fq.getName().contains("_R1.")} == false) {
            ParamLogger.throwError("Library ${it[0]} does not have read1 fastq file")
        }

        if (it[3].any {fq -> fq.getName().contains("_R2_") or fq.getName().contains("_R2.")} == false) {
            ParamLogger.throwError("Library ${it[0]} does not have read2 fastq file")
        }

        if (it[3].any {fq -> fq.getName().contains("_I1_") or fq.getName().contains("_I1.")} == false) {
            ParamLogger.throwError("Library ${it[0]} does not have index fastq file")
        }
    }
    
    BarcodeDemux(samplesCsv, libJson.getParent(), libJson.getName(), fqSamples)
    barcodeDemuxFastqs = (BarcodeDemux.out.read1Fastq.groupTuple()).join(BarcodeDemux.out.read2Fastq.groupTuple(), remainder: true)
    // Check that there exists fastq files(R1, R2) post bcParser, for every sample that corresponds to a libName in samples.csv
    // bcParserOutCheck -> (libName, [sample names corresponding to that libName], [[R1 fastq files]], [[R2 fastq files]])
    bcParserOutCheck = libNameAndSamples.groupTuple().join((barcodeDemuxFastqs), remainder: true)
    bcParserOutCheck.dump(tag: 'bcParserOutCheck')
    bcParserOutCheck.map {
        if (it[2] == null) {
            ParamLogger.throwError("Library ${it[0]} does not have any passing R1 and R2 fastq files after bcParser")
        }
        // iterate over every sample name corresponding to a libName
        for (sample in it[1]) {
            // extract the sample name from the id
            sampleName = sample.tokenize('.')[0]
            exists = false
            combinedR1R2Fastqs = it[2].flatten() + it[3].flatten()
            // iterate over every fastq file for that libName and check that the sample name is present in the fastq file name
            for (file in combinedR1R2Fastqs) {
                if (file.getName().toString().contains(sampleName)) {
                    exists = true
                }
            }
            if (!exists) {
                log.warn("Sample ${sampleName} does not have passing R1 and R2 files from bcParser.")
            }
        }
    }
    
    // read1Fqs -> (libName, sample name, subsample name, read1 fastq files)
    read1Fqs = BarcodeDemux.out.read1Fastq.transpose().map{
        constructSampleName(it, params.splitFastq)
    }.groupTuple(by:[0,1,2])
 
    // read2Fqs -> (libName, sample name, subsample name, read2 fastq files)
    read2Fqs = BarcodeDemux.out.read2Fastq.transpose().map{
        constructSampleName(it, params.splitFastq)
    }.groupTuple(by:[0,1,2])

    // readFqs -> (libName, sample name, subsample name, [read1 fastq files], [read2 fastq files])
    readFqs = read1Fqs.join(read2Fqs, by:[0,1,2])
    readFqs.dump(tag: 'readFqs')

    // if splitFastq = true: demuxFqs -> [sample name.well coordinate, [demuxed R1 fastq files], [demuxed R2 fastq files]]
    // if splitFastq = false: demuxFqs -> [sample name, [demuxed fastq files]]
    // Drop libName and matchingKey from the tuple
    demuxFqs = readFqs.map{ tuple(it[2], it[3], it[4]) }
    demuxFqs.dump(tag: 'demuxFqs')

    organizedDemux = BarcodeDemux.out.metrics.groupTuple()
    // merge all bcParser metrics together into one json file for each library
    MergeDemux(organizedDemux, libJson.getParent(), libJson.getName())
    
    // Run trim_galore on fastq files post demux
    Trim(demuxFqs)
    
    // if splitFastq = true: Trim.out.fastq -> [sample name.well coordinate, [trimmed fastq files from that well coordinate for that sample]]
    // if splitFastq = false: Trim.out.fastq -> [sample name, [trimmed fastq files for that sample]]
    Trim.out.fastq.dump(tag: 'Trim.out.fastq')

emit:
    trimLog = Trim.out.stats // cutadapt statistics
    trimFastq = Trim.out.fastq // Trimmed fastq files
    mergedBcParserMetrics = MergeDemux.out.barcodeMetrics // bcParser metrics for all libName
}
