nextflow.enable.dsl=2
import groovy.json.JsonSlurper
import groovy.json.JsonOutput

// Load .json file into object
def loadJson(json) {
	def jsonSlurper  = new JsonSlurper()
	return jsonSlurper.parse( json )
}

def expandPath(path, baseDir) {
	if (path == null) { return null}
	if (path =~ /^s3/) { return file(path)}
	return baseDir.resolve(path)
}

def loadGenome(json) {
	def baseDir = json.getParent()
	genome = loadJson(json)
	genome.bsbolt_index = expandPath(genome.bsbolt_index, baseDir)
	genome.bsbolt_chrs = expandPath(genome.bsbolt_chrs, baseDir)
	genome.genomeTiles = expandPath(genome.genomeTiles, baseDir)
	return genome
}

// Return 0 for null or non-integer strings
def toIntOr0(str) {
	if ((str != null) && str.isInteger())
		return str as int;
	else
		return 0
}

fqDir = params.fastqDir
runDir = params.runFolder
fastqSamplesheet = params.fastqSamplesheet
genome = loadGenome(file(params.genome))
libJson = expandPath(params.libStructure, file("${projectDir}/references/"))
bam1Dir = params.bam1Dir
bam2Dir = params.bam2Dir

// Prepare samples.csv with defaults, rename legacy columns, etc.
process regularizeSamplesCsv {
input: 
	path("samples.in.csv")
output: 
	path("samples.csv")
publishDir "${params.outDir}", mode: 'copy'
label 'py_process'
cache 'deep'
script:
	opts=""
	if (params.splitBarcodeParsing) {
		opts += "--splitBarcodeParsing"
	}
"""
	regularizeSamplesCsv.py samples.in.csv $opts > samples.csv
"""
}

// Create a bcl-convert samplesheet for libraries in samples.json
process makeBclConvertSamplesheet {
input: 
	path(samplesCsv)
	path(libStructDir) // Directory containing the library structure definition (barcode sequence lists are loaded from here)
	val(libStructName) // Filename of the library structure definition .json
	path(runinfo) // RunInfo.xml from Sequencer RunFolder (Read-lengths etc.)
output: 
	path("samplesheet.csv")
publishDir "${params.outDir}/fastq", mode: 'copy'
label 'py_process'
script:
	opts = ""
	libJson = "$libStructDir/$libStructName"
	if (params.splitFastq) {
		opts = "--splitFastq"
	}
"""
	bclConvertSheet.py $samplesCsv $libJson $runinfo $opts > samplesheet.csv
"""
}


/*Run bcl-convert, used when starting from a sequencing run-folder
  Requires a separate bcl-convert samplesheet*/
process bclconvert {
input: 
	path(run)
	path(samplesheet)
	val(bclConvertParams)
output: 
	path("fastq/*fastq.gz"), emit: fastq
	path("fastq/Reports/*"), emit: stats
	publishDir "${params.outDir}/", pattern: 'fastq/Reports/*', mode: 'copy'
	publishDir "${params.outDir}/", pattern: 'fastq/*.fastq.gz', enabled: params.fastqOut
script:
"""
	bcl-convert --sample-sheet $samplesheet --bcl-input-directory $run \
    --bcl-num-conversion-threads 4 --bcl-num-compression-threads 4 --bcl-num-decompression-threads 4 \
    $bclConvertParams --output-directory fastq
"""
}

process fastqc {
input:
	path(fqFiles)
	path(adapters)
output:
	path("fastqc/*.html")
	publishDir "${params.outDir}", mode: 'copy'
script:
"""
	mkdir -p fastqc
	fastqc --adapters ${adapters} -o fastqc ${fqFiles.join(" ")} --threads ${task.cpus} -f fastq --extract
"""
}

process barcodeDemux { 
input: 
	path(sheet) // samples.csv 
	path(libStructDir) // Directory containing the library type definition file (barcode sequence lists are loaded from here) 
	val(libStructJson) // Filename of the library definition .json 
	tuple(val(libName), val(laneNum), path(fqFiles)) // Input fastq file 
output: 
	tuple(val(libName), path("$outDir/*_S[1-9]*.fastq.gz"), emit: fastq) 
	path("$outDir/*_S0_*.fastq.gz"), emit: unknown optional true 
	path("$outDir/*.tsv") 
	tuple(val(libName), path("$outDir/metrics.json"), emit: metrics)
	publishDir "${params.outDir}/fastq/", pattern: "$outDir/*{txt,tsv,json}", mode: 'copy'
	publishDir "${params.outDir}/fastq/", pattern: "$outDir/*gz", enabled: params.fastqOut
	tag "$libName" 
script: 
	outDir = "${libName}.${laneNum}.demux" 
	libStruct = "$libStructDir/$libStructJson" 
""" 
	bc_parser --discard-readnames -v --write-fastq --demux $sheet --lib-name $libName --fastq-num $laneNum --reads ${fqFiles.join(" ")} --lib-struct $libStruct --out $outDir
""" 
}

process trim {
input:
	tuple(val(sample), path(pairs1), path(pairs2)) // Input read pairs (fastq)
output:
	tuple(val(sample), path("${sample}_val_{1,2}.fq.gz"), emit: fastq)
	tuple(val(sample), path("${sample}.fq.gz_trimming_report_R2.txt"), emit: report)
	publishDir "${params.outDir}/trim/${sample}", pattern: "${sample}_val_{1,2}.fq.gz", enabled: params.trimOut
script:
	tthreads = task.cpus
"""
	cat ${pairs1} > ${sample}_R1.fq.gz
	cat ${pairs2} > ${sample}_R2.fq.gz
	trim_galore -a CTATCTCTTATA -a2 AGATCGGAAGAGC --three_prime_clip_R2 10 --basename ${sample} -j $tthreads --paired ${sample}_R1.fq.gz ${sample}_R2.fq.gz
	rm ${sample}_R{1,2}.fq.gz
	mv ${sample}_R2.fq.gz_trimming_report.txt ${sample}.fq.gz_trimming_report_R2.txt
	mv ${sample}_R1.fq.gz_trimming_report.txt ${sample}.fq.gz_trimming_report_R1.txt
"""
}

process align {
input: 
	path(index) // Bismark index directory
	tuple(val(sample), path(pairs))
output: 
	tuple(val(sample), path("${sample}.bam"), emit: bam)
	tuple(val(sample), path("${sample}.bsbolt.log"), emit: log)
	publishDir "${params.outDir}/bsbolt/${sample}", pattern: "${sample}.bam", enabled: params.bamOut
script:
	athreads = task.cpus - 2 
	othreads = 2
"""
	bsbolt Align -F1 ${pairs[1]} -F2 ${pairs[0]} -t $athreads -OT $othreads -O $sample -DB $index >> ${sample}.bsbolt.log 2>> ${sample}.bsbolt.log
"""
}

process mergeBam {
input: 
	tuple(val(sample), val(tgmt), path("bam1.${tgmt}.bam"), path("bam2.${tgmt}.bam")) // Input bam files
	path(rgs)
output: 
	tuple(val("${sample}.${tgmt}"), path("${sample}.${tgmt}.bam"), emit: bam)
	publishDir "${params.outDir}/bamMerge/${sample}.${tgmt}", pattern: "${sample}.${tgmt}.bam", enabled: params.bamMergeOut
script:
	sthreads = task.cpus
"""
	samtools addreplacerg -r ID:bam1 -r LB:bam1 bam1.${tgmt}.bam | samtools sort -@ $sthreads -o bam1.bam -
	samtools addreplacerg -r ID:bam2 -r LB:bam2 bam2.${tgmt}.bam | samtools sort -@ $sthreads -o bam2.bam -
	samtools merge --threads $sthreads -o ${sample}.${tgmt}.bam bam1.bam bam2.bam
"""
}

process dedup {
input: 
	tuple(val(sample), path(bam))
	val(minUniq)
	val(minUniqTotal)
	path(chroms)
output: 
	tuple(val(sample), path("${sample}.dedup.nsrt.bam"), emit: bam)
	tuple(val(sample), path("${sample}.cell_stats.tsv"), emit: complexity)
	tuple(val(sample), path("${sample}.dedup_stats.tsv"), emit: stats)
	tuple(val(sample), path("${sample}.fragment_hist.tsv"), emit: fragHist)
	tuple(val(sample), path("${sample}.uniqcount.annot"), emit: annot)
	publishDir "${params.outDir}/bamDeDup/${sample}", pattern: "${sample}.dedup.nsrt.bam", enabled: params.bamDedupOut
script:
	sthreads = task.cpus
"""
	samtools sort --threads $sthreads $bam -o ${sample}.coodsort.bam
	sc_dedup ${sample}.coodsort.bam --duplicate-key ${params.dedupKey} --min-mapq ${params.minMapq} --barcode-input QnameStart --out-prefix ${sample} --write-threads $sthreads --genome $chroms
	samtools sort --threads $sthreads -n ${sample}.dedup.bam -o ${sample}.dedup.nsrt.bam
	grep -v Barcode ${sample}.cell_stats.tsv | awk -vUniq=${minUniq} -vPct=${minUniqTotal} '{if(\$4 >= Uniq) {print \$1}}' > ${sample}.uniqcount.annot
"""
}

// will produce a large amount of files in the sample dir 3 per files per BC in the bam that should remove after bed files per chrom are created
process extract {
input: 
	tuple(val(sample), path(annot), path(bam))
output: 
	tuple(val(sample), path("${sample}*CG.chroms"), emit: metGC )
	tuple(val(sample), path("${sample}*CH.chroms"), emit: metCH)
	tuple(val(sample), path("${sample}*cellInfo.txt"), emit: report)
script:
	ethreads = task.cpus
"""
	sciMET_BSBolt_extract.pl -O $sample -t $ethreads -C $annot $bam
	rm -rf $sample
"""
}

process sortExtractCG {
input:
	tuple(val(sample), path(metGC))
output: 
	tuple(val(sample), path("${sample}*CG.chroms.sort"), emit: metGC )
	publishDir "${params.outDir}/cg_sort_cov/${sample}", pattern: "${sample}*CG.chroms.sort", enabled: params.covOut
script:
	ethreads = task.cpus
"""
	mkdir ${sample}.CG.chroms.sort
	for f in ${metGC}/*bed; do outFile=`basename \$f`; sort -k2,2n \$f | gzip > ${sample}.CG.chroms.sort/\${outFile}.gz; done
"""
}

process sortExtractCH {
input:
	tuple(val(sample), path(metCH))
output: 
	tuple(val(sample), path("${sample}*CH.chroms.sort"), emit: metCH )
	publishDir "${params.outDir}/ch_sort_cov/${sample}", pattern: "${sample}*CH.chroms.sort", enabled: params.covOut
script:
	ethreads = task.cpus
"""
	mkdir ${sample}.CH.chroms.sort
	for f in ${metCH}/*bed; do outFile=`basename \$f`; sort -k2,2n \$f | gzip > ${sample}.CH.chroms.sort/\${outFile}.gz; done
"""
}


process merge_demux {
input:
    tuple(val(libName), path("demux_metrics*.json"))
    path(libStructDir) // Directory containing the library type definition file (barcode sequence lists are loaded from here) 
    val(libStructJson) // Filename of the library definition .json 
output:
    tuple(val(libName), path("${libName}.metrics.json"), emit: barcode_metrics)
publishDir "${params.outDir}/library_barcode_metrics", mode:'copy'
label 'py_process'
script:
	libStruct = "$libStructDir/$libStructJson" 
    """
    mergeBCparserOutput.py --bc_jsons demux_metrics* --lib_json $libStruct --libName $libName
    """	
}


process generateReport {
input:
	tuple(val(sample), path(cell_stats), path(frags), path(cell_info_stats), path(mapping_stats), path(trimming_stats), val(threshold))
	path(references)
output:
	tuple(val(sample), path("${sample}/csv/${sample}.allCells.csv") , emit: complex )
	tuple(val(sample), path("${sample}/${sample}.report.html"), emit: report)
	tuple(val(sample), path("filter_cells/${sample}*.cells"), emit: annot)
	tuple(val(sample), path("${sample}/csv/${sample}.passingCellsMapMethylStats.csv"), emit: cell_map_methyl_stats)
	tuple(val(sample), path("${sample}/csv/${sample}.passingCellSummaryStats.csv"), emit: passing_cell_methyl_stats)
	tuple(val(sample), path("${sample}/csv/${sample}.trimming_stats.csv"), emit: trimming_stats, optional: true)
	tuple(val(sample), path("${sample}/csv/${sample}.mapping_stats.csv"), emit: mapping_stats, optional: true)
	tuple(val(sample), path("${sample}/csv/${sample}.summaryStats.csv"), emit: summary_stats)
	tuple(val(sample), path("${sample}/png"))
	tuple(val(sample), path("${sample}/csv"))
publishDir "${params.outDir}/report", pattern: "${sample}/*.html", mode:'copy'
publishDir "${params.outDir}/report", pattern: "${sample}/png", mode:'copy'
publishDir "${params.outDir}/report", pattern: "${sample}/csv", mode:'copy'
label 'py_process'
script:
	"""
	library_complexity_and_stats.py --sampleName ${sample} --outDir ${sample} --tgmt_barcodes ${references}/tgmt.txt --topCellPercentage ${params.topCellPercentage} --threshold $threshold\
	--minCellRatio ${params.minCellRatio} --minReads ${params.minReads} --minUniqTotal ${params.minUniqTotal} --maxUniqTotal ${params.maxUniqTotal} --i5_barcodes ${references}/i5.txt --i7_barcodes ${references}/i7.txt
	"""
}

process combinedSampleReport {
input:
	path(complexity_stats)
	path(methyl_stats)
	path(passing_met_stats)
	tuple(val(libName), path("demuxMetrics.json"))
	path(references)
output:
	tuple(val(libName), path("library_report/library.${libName}.report.html"), emit: combined_report)
	tuple(val(libName), path("library_report/library.${libName}.combinedPassingCellStats.csv"), emit: combined_passing_stats)
	tuple(val(libName), path("library_report/png"))
	tuple(val(libName), path("library_report/csv"))
publishDir "${params.outDir}", pattern: "library_report/png", mode:'copy'
publishDir "${params.outDir}", pattern: "library_report/csv", mode:'copy'
publishDir "${params.outDir}", pattern: "library_report/*.html", mode:'copy'
publishDir "${params.outDir}", pattern: "library_report/*.csv", mode:'copy'
label 'py_process'
script:
	"""
	combined_sample_report.py --libraryName ${libName} --outDir library_report --complexity_stats ${complexity_stats} --methyl_stats ${methyl_stats} \
	--passing_met_stats ${passing_met_stats} --i5_barcodes ${references}/i5.txt --i7_barcodes ${references}/i7.txt --tgmt_barcodes ${references}/tgmt.txt
	"""
}

// Binned mtx generation
process mtxCG {
input:
	path(tiles) // Non overlapping genome bins
	tuple(val(sample), path(cells), path(metGC))
output: 
	tuple(val(sample), path("${sample}*CG.mtx"), emit: mtx, optional: true )
	tuple(val(sample), path("${sample}*CG.cov.mtx"), emit: cov, optional: true )
	tuple(val(sample), path("${sample}*CG.ratio.mtx"), emit: ratio, optional: true )
	tuple(val(sample), path("${sample}*CG.score.mtx"), emit: score, optional: true )
script:
"""
	sciMET_meth2mtx.pl -F $metGC -B $tiles -O ${sample} -C $cells
"""
}

process mtxCH {
input:
	path(tiles) // Non overlapping genome bins
	tuple(val(sample), path(cells), path(metCH))
output: 
	tuple(val(sample), path("${sample}*CH.mtx"), emit: mtx, optional: true )
	tuple(val(sample), path("${sample}*CH.cov.mtx"), emit: cov, optional: true )
	tuple(val(sample), path("${sample}*CH.ratio.mtx"), emit: ratio, optional: true )
	tuple(val(sample), path("${sample}*CH.score.mtx"), emit: score, optional: true )
script:
"""
	sciMET_meth2mtx.pl -F $metCH -B $tiles -O ${sample} -C $cells
"""
}

// merge mtx files
process mergeMtxCG {
input:
	tuple(val(sample), path(metGC))
	tuple(val(sample), path(metGC_score))
	tuple(val(sample), path(metGC_cov))
	tuple(val(sample), path(metGC_ratio))
output: 
	tuple(val(sample), path("${sample}.CG.mtx"), emit: mtx )
	tuple(val(sample), path("${sample}.CG.score.mtx"), emit: score_mtx )
	tuple(val(sample), path("${sample}.CG.cov.mtx"), emit: cov_mtx )
	tuple(val(sample), path("${sample}.CG.ratio.mtx"), emit: ratio_mtx )
publishDir "${params.outDir}/matrix", pattern: '*{mtx}', mode: 'copy'
script:
"""
	sciMET_mergeMtx.pl ${sample}.CG.mtx $metGC
	sciMET_mergeMtx.pl ${sample}.CG.score.mtx $metGC_score
	sciMET_mergeMtx.pl ${sample}.CG.cov.mtx $metGC_cov
	sciMET_mergeMtx.pl ${sample}.CG.ratio.mtx $metGC_ratio
"""
}

// merge CH mtx files
process mergeMtxCH {
input:
	tuple(val(sample), path(metCH))
	tuple(val(sample), path(metCH_score))
	tuple(val(sample), path(metCH_cov))
	tuple(val(sample), path(metCH_ratio))
output: 
	tuple(val(sample), path("${sample}.CH.mtx"), emit: mtx )
	tuple(val(sample), path("${sample}.CH.score.mtx"), emit: score_mtx )
	tuple(val(sample), path("${sample}.CH.cov.mtx"), emit: cov_mtx )
	tuple(val(sample), path("${sample}.CH.ratio.mtx"), emit: ratio_mtx )
publishDir "${params.outDir}/matrix", pattern: '*.mtx', mode: 'copy'
script:
"""
	sciMET_mergeMtx.pl ${sample}.CH.mtx $metCH
	sciMET_mergeMtx.pl ${sample}.CH.score.mtx $metCH_score
	sciMET_mergeMtx.pl ${sample}.CH.cov.mtx $metCH_cov
	sciMET_mergeMtx.pl ${sample}.CH.ratio.mtx $metCH_ratio
"""
}

workflow {
	fqs = null
	regularizeSamplesCsv(file(params.samples))
	samples = regularizeSamplesCsv.out.splitCsv(header:true, strip:true)
	// Indicates starting from bam files
	if (params.bam1Dir != null && params.bam2Dir != null) {
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
		b1.dump(tag: 'b1')
		b2.dump(tag: 'b2')
		bams = b1.join(b2, by: [0,1], failOnDuplicate: true)
		bams.dump(tag: 'bams')
		mergeBam(bams, params.bamRgHeader)
		mergeBam.out.bam.dump(tag: 'mergeBam.out.bam')
		dedup(mergeBam.out.bam, params.minUniqCount, params.minUniqTotal, genome.bsbolt_chrs)
		trimming_and_mapping = false
	}
	else {
		if (runDir != null) {
			if (fastqSamplesheet == null) {
				makeBclConvertSamplesheet(regularizeSamplesCsv.out, libJson.getParent(), libJson.getName(), 
							  file("$runDir/RunInfo.xml", checkIfExists:true))
				fqSheet = makeBclConvertSamplesheet.out
			}
			else {
				fqSheet = file(fastqSamplesheet)
			}
			bclconvert(file(runDir), fqSheet, params.bclConvertParams)
			fqs = bclconvert.out.fastq.flatten().filter({ !it.name.contains("Undetermined") })
		} else if (fqDir != null) {
			fqs = Channel.fromPath("$fqDir/*fastq.gz", checkIfExists: true)
		} else {
			throwError("Must specify either 'runFolder' or 'fastqDir'")
		}
		// fqs -> [single fastq file]
		// n channels for n fastq files
		fqs.dump(tag:'fqs')
		if (params.fastqc) {
			fastqc(fqs, params.adapters)
		}
		// sampleFqs -> [sample name, identifier[either sample name or lane number], [fastq files corresponding to that sample name and identifier]]
		sampleFqs = fqs.map { file ->	
			def ns = file.getName().toString().tokenize('_')
			return tuple(ns[0], ns[2].substring(1), file) //(ns[0], ns[1], ns[3].substring(1) // if it were split by i5 and lane
		}.groupTuple(by:[0,1])
		sampleFqs.dump(tag: 'sampleFqs')
		libName = sampleFqs.map {
			it[0]
		}
		barcodeDemux(regularizeSamplesCsv.out, libJson.getParent(), libJson.getName(), sampleFqs)
		// barcodeDemux.out.fastq -> [sample name, [all demuxed fastq files corresponding to one sample name and identifier]]
		barcodeDemux.out.fastq.dump(tag: 'barcodeDemux.out.fastq')
		organized_demux = barcodeDemux.out.metrics.groupTuple()
		merge_demux(organized_demux, libJson.getParent(), libJson.getName())
		demuxFqs = barcodeDemux.out.fastq.flatMap({it[1]}).map { file ->
			def ns = file.getName().toString().tokenize('_')
			return tuple(ns[0], ns[3].substring(0,3), file)
		}.groupTuple(by:[0,1])
		// pairedFqs -> [sample name.well coordinate, [fastq files from that well coordinate for that sample]]
		pairedFqs = demuxFqs.map{tuple(it[0], it[2][0], it[2][1])}.groupTuple()
		pairedFqs.dump(tag: 'pairedFqs')
		trim(pairedFqs)
		trim.out.fastq.dump(tag: 'trim.out.fastq')
		align(genome.bsbolt_index, trim.out.fastq)
		dedup(align.out.bam, params.minUniqCount, params.minUniqTotal, genome.bsbolt_chrs)
		trimming_and_mapping = true
	}
	annot = dedup.out.annot.filter ({ it[1].size() > 0 }).join(dedup.out.bam)
	annot.dump(tag: 'annot')
 	extract(annot)
 	sortExtractCG(extract.out.metGC)
	sortExtractCG.out.metGC.dump(tag: 'sortExtractCG.out.metGC')
	sortExtractCH(extract.out.metCH)
	sortExtractCH.out.metCH.dump(tag: 'sortExtractCH.out.metCH')
 	
	complexTabs = dedup.out.complexity.map({it[1]}).map { file ->
 		def ns = file.getName().toString().tokenize('.').get(0)
 		return tuple(ns, file)
 	}.groupTuple()
	fragTabs = dedup.out.fragHist.map({it[1]}).map { file ->
 		def ns = file.getName().toString().tokenize('.').get(0)
 		return tuple(ns, file)
 	}.groupTuple()
	readsFrags = complexTabs.join(fragTabs)
	readsFrags.dump(tag: 'readFrags')
 	sampleExtractStats = extract.out.report.map({it[1]}).map { file ->
		def ns = file.getName().toString().tokenize('.')
		return tuple(ns.get(0), file)
	}.groupTuple()
	sampleExtractStats.dump(tag: 'sampleExtractStats')
	
	if (trimming_and_mapping) {
		// mapping stats
		sampleMapStats = align.out.log.map({it[1]}).map { file ->
			def ns = file.getName().toString().tokenize('.')
			return tuple(ns.get(0), file)
		}.groupTuple()
		sampleMapStats.dump(tag: 'sampleMapStats')
		trim.out.report.dump(tag: 'trim.out.report')
		sampleTrimStats = trim.out.report.map({it[1]}).map { file ->
			def ns = file.getName().toString().tokenize('.')
			return tuple(ns.get(0), file)
		}.groupTuple()
		sampleTrimStats.dump(tag: 'sampleTrimStats')
		report_stats = readsFrags.join(sampleExtractStats).join(sampleMapStats).join(sampleTrimStats)
		report_stats.dump(tag: 'report_stats')
		demuxMetrics = merge_demux.out.barcode_metrics
	}
	else {
		report_stats = readsFrags.join(sampleExtractStats)
		report_stats = report_stats.map {
			tuple(it[0], it[1], it[2], it[3], [], [])
		}
		report_stats.dump(tag: 'report_stats')
		demuxMetrics = tuple("combined", [])
	}
	
	threshold = samples.map{[it.sample, toIntOr0(it.threshold)]}
	report_stats = report_stats.join(threshold)
	
	generateReport(report_stats, libJson.getParent())
	
	allMetStats = generateReport.out.cell_map_methyl_stats.map({it[1]}).collect()
	allMetStats.dump(tag: 'allMetStats')
	allComplexStats = generateReport.out.complex.map({it[1]}).collect()
	allComplexStats.dump(tag: 'allComplexStats')
	passing_cell_met_stats = generateReport.out.passing_cell_methyl_stats.map({it[1]}).collect()
	
	combinedSampleReport(allComplexStats, allMetStats, passing_cell_met_stats, demuxMetrics, libJson.getParent())
	//Different logic required for tokenization when we split in bcParser vs when we don't
	if (params.splitBarcodeParsing) {
		cgPass = generateReport.out.annot.filter ({ it[1].size() > 0 }).flatMap({it[1]}).map {file ->
			def ns = file.getName().toString().tokenize('.')
			return tuple(ns[0] + '.' + ns[1], file)
		}.groupTuple().join(sortExtractCG.out.metGC)
		chPass = generateReport.out.annot.filter ({ it[1].size() > 0 }).flatMap({it[1]}).map {file ->
	 		def ns = file.getName().toString().tokenize('.')
	 		return tuple(ns[0] + '.' + ns[1], file)
		}.groupTuple().join(sortExtractCH.out.metCH)
	}
	else {
		cgPassTmp = generateReport.out.annot.filter ({ it[1].size() > 0 }).flatMap({it[1]}).map {file ->
	 		    def ns = file.getName().toString().tokenize('.')
	 		    return tuple(ns[0], file)
		}
		cgPass = sortExtractCG.out.metGC.cross(cgPassTmp)
		cgPass = cgPass.map {
			def ns = it[1][1].getName().toString().tokenize('.')
			return tuple(ns[0] + '.' + ns[1], it[1][1], it[0][1])
		}
		chPassTmp = generateReport.out.annot.filter ({ it[1].size() > 0 }).flatMap({it[1]}).map {file ->
	 		    def ns = file.getName().toString().tokenize('.')
	 		    return tuple(ns[0], file)
		}
		chPass = sortExtractCH.out.metCH.cross(chPassTmp)
		chPass = chPass.map {
			def ns = it[1][1].getName().toString().tokenize('.')
			return tuple(ns[0] + '.' + ns[1], it[1][1], it[0][1])
		}
	}

	// cgPass -> [sample name.well coordinate, [barcodes for passing cells for that well coordinate], sorted bed file] 
	cgPass.dump(tag: 'cgPass')
 	mtxCG(genome.genomeTiles,cgPass)
	
	// chPass -> [sample name.well coordinate, [barcodes for passing cells for that well coordinate], sorted bed file] 
	chPass.dump(tag: 'chPass')
	mtxCH(genome.genomeTiles,chPass)
	
	// merge CG mtx files
	sampleCGmtx = mtxCG.out.mtx.map({it[1]}).map { file ->
 		def ns = file.getName().toString().tokenize('.').get(0)
 		return tuple(ns, file)
 	}.groupTuple()
	sampleCGscore = mtxCG.out.score.map({it[1]}).map { file ->
		def ns = file.getName().toString().tokenize('.').get(0)
		return tuple(ns, file)
	}.groupTuple()
	sampleCGcov = mtxCG.out.cov.map({it[1]}).map { file ->
		def ns = file.getName().toString().tokenize('.').get(0)
		return tuple(ns, file)
	}.groupTuple()
	sampleCGratio = mtxCG.out.ratio.map({it[1]}).map { file ->
		def ns = file.getName().toString().tokenize('.').get(0)
		return tuple(ns, file)
	}.groupTuple()

	mergeMtxCG(sampleCGmtx, sampleCGscore, sampleCGcov, sampleCGratio)
	
	// merge CH mtx files
	sampleCHmtx = mtxCH.out.mtx.map({it[1]}).map { file ->
		def ns = file.getName().toString().tokenize('.').get(0)
		return tuple(ns, file)
	}.groupTuple()
	sampleCHscore = mtxCH.out.score.map({it[1]}).map { file ->
		def ns = file.getName().toString().tokenize('.').get(0)
		return tuple(ns, file)
	}.groupTuple()
	sampleCHcov = mtxCH.out.cov.map({it[1]}).map { file ->
		def ns = file.getName().toString().tokenize('.').get(0)
		return tuple(ns, file)
	}.groupTuple()
	sampleCHratio = mtxCH.out.ratio.map({it[1]}).map { file ->
		def ns = file.getName().toString().tokenize('.').get(0)
		return tuple(ns, file)
	}.groupTuple()
	
	mergeMtxCH(sampleCHmtx, sampleCHscore, sampleCHcov, sampleCHratio)
}

// Function that gets invoked when workflow completes
// Publishes a file named workflow_info.json that contains information about the pipeline execution
workflow.onComplete {
	testing = false
	def data = ["Workflow Information":
		    ["Execution status": "${ workflow.success ? 'OK' : 'failed' }",
		     "Pipeline completion timestamp": "$workflow.complete",
		     "Git repo URL": "$workflow.repository",
		     "Configuration files": "$workflow.configFiles",
		     "Container": "$workflow.container",
		     "Command line executed": "$workflow.commandLine",
		     "Configuration profile": "$workflow.profile",
		     "Start timestamp": "$workflow.start",
		     "Stop timestamp": "$workflow.complete",
		     "Run name": "$workflow.runName",
		     "Revision": "$workflow.revision",
		     "Commit ID": "$workflow.commitId",
		     "Duration": "$workflow.duration"
                    ]
		   ]
	def params_data = ["Parameters": [:]]
	for (p in params) {
		if (!p.key.contains('-')) {
			params_data."Parameters".put("$p.key", "$p.value")
		}
		if (p.key.equals('testing_run')) {
			testing = true
			data."Workflow Information".put("Exit status", "$workflow.exitStatus")
			data."Workflow Information".put("Error message", "$workflow.errorMessage")
		}
    	}
	def reference_data = ["Reference Genome": [:]]
	for (p in genome) {
		reference_data."Reference Genome".put("$p.key", "$p.value")
	}
	def manifest_data = ["Manifest": [:]]
	for (p in workflow.manifest.getProperties()) {
		p = p.toString()
		def split_str = p.split("=")
		if (split_str[0].equals("name") || split_str[0].equals("version")) {
			manifest_data."Manifest".put(split_str[0], split_str[1])
		}
	}
	def json_str = JsonOutput.toJson(data+params_data+reference_data+manifest_data)
	def json_beauty = JsonOutput.prettyPrint(json_str)
	workflow_info = file("$params.outDir/workflow_info.json")
	workflow_info.write(json_beauty)
}

