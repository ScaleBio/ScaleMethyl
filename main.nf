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
	genome.genomeTilesCh = expandPath(genome.genomeTilesCh, baseDir)
	genome.tssWin = expandPath(genome.tssWin, baseDir)
	genome.backgroundWin = expandPath(genome.backgroundWin, baseDir)
	return genome
}

// Return 0 for null or non-integer strings
def toIntOr0(str) {
	if ((str != null) && str.isInteger())
		return str as int;
	else
		return 0
}

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

// Run bcl-convert, used when starting from a sequencing run-folder
// Requires a separate bcl-convert samplesheet
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

// Run fastqc on either input fastq files or bcl-convert generated fastq files
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

// Demux fastq files and do barcode error correction
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

// Merge bc_parser stats from multiple jobs for a single library
process merge_demux {
input:
    tuple(val(libName), path("demux_metrics*.json")) // bc_parser metrics from all bc_parser runs for a library
    path(libStructDir) // Directory containing the library type definition file (barcode sequence lists are loaded from here) 
    val(libStructJson) // Filename of the library definition .json 
output:
    tuple(val(libName), path("${libName}.metrics.json"), emit: barcode_metrics)
publishDir "${params.outDir}/library_barcode_metrics", mode:'copy'
tag "$libName"
label 'py_process'
script:
	libStruct = "$libStructDir/$libStructJson" 
    """
    mergeBCparserOutput.py --bc_jsons demux_metrics* --lib_json $libStruct --libName $libName
    """	
}

// Trim demuxed fastq files produced by bc_parser and run multiqc
process trim {
input:
	tuple(val(sample), path(pairs1), path(pairs2))
output:
	tuple(val(sample), path("${sample}_val_{1,2}.fq.gz"), emit: fastq)
	tuple(val(sample), path("${sample}.fq.gz_trimming_report_R2.txt"), emit: report)
	tuple(val(sample), path("${sample}.fq.gz_trimming_report_R1.txt"))
	path("multiqc_report.html")
publishDir "${params.outDir}/trim/${sample}", pattern: "${sample}_val_{1,2}.fq.gz", enabled: params.trimOut
publishDir "${params.outDir}/trim/${sample}", pattern: "${sample}.fq.gz_trimming_report_R{1,2}.txt"
publishDir "${params.outDir}/multiqc/${sample.tokenize('.')[0]}", pattern: "multiqc_report.html", saveAs: {filename -> "${sample}_$filename"}, mode: 'copy'
script:
	tthreads = task.cpus
"""
	cat ${pairs1} > ${sample}_R1.fq.gz
	cat ${pairs2} > ${sample}_R2.fq.gz
	trim_galore --fastqc -a CTATCTCTTATA -a2 AGATCGGAAGAGC --three_prime_clip_R2 10 --basename ${sample} -j $tthreads --paired ${sample}_R1.fq.gz ${sample}_R2.fq.gz
	rm ${sample}_R{1,2}.fq.gz
	mv ${sample}_R2.fq.gz_trimming_report.txt ${sample}.fq.gz_trimming_report_R2.txt
	mv ${sample}_R1.fq.gz_trimming_report.txt ${sample}.fq.gz_trimming_report_R1.txt
	multiqc .
"""
}

// Align trimmed fastq files to a reference genome using bsbolt
process align {
input: 
	path(index) // Bismark index directory
	tuple(val(sample), path(pairs))
output: 
	tuple(val(sample), path("${sample}.bam"), emit: bam)
	tuple(val(sample), path("${sample}.bsbolt.log"), emit: log)
publishDir "${params.outDir}/alignment/${sample}", pattern: "${sample}.bam", enabled: params.bamOut
tag "$sample"
script:
	athreads = task.cpus - 2 
	othreads = 2
"""
	bsbolt Align -F1 ${pairs[1]} -F2 ${pairs[0]} -t $athreads -OT $othreads -O $sample -DB $index >> ${sample}.bsbolt.log 2>> ${sample}.bsbolt.log
"""
}

// Merge two bam files together using samtools
process mergeBam {
input: 
	tuple(val(sample), val(tgmt), path("bam1.${tgmt}.bam"), path("bam2.${tgmt}.bam"))
	path(rgs)
output: 
	tuple(val("${sample}.${tgmt}"), path("${sample}.${tgmt}.bam"), emit: bam)
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

// Deduplicate reads in aligned bam files using sc_dedup
process dedup {
input: 
	tuple(val(sample), path(bam))
	val(minUniq)
	path(chroms)
output: 
	tuple(val(sample), path("${sample}.dedup.nsrt.bam"), emit: bam)
	tuple(val(sample), path("${sample}.cell_stats.tsv"), emit: complexity)
	tuple(val(sample), path("${sample}.dedup_stats.tsv"), emit: stats)
	tuple(val(sample), path("${sample}.fragment_hist.tsv"), emit: fragHist)
	tuple(val(sample), path("${sample}.uniqcount.annot"), emit: annot)
publishDir "${params.outDir}/bamDeDup/${sample}", pattern: "${sample}.dedup.nsrt.bam", enabled: params.bamDedupOut
tag "$sample"
script:
	sthreads = task.cpus
"""
	samtools sort --threads $sthreads $bam -o ${sample}.coodsort.bam
	sc_dedup ${sample}.coodsort.bam --duplicate-key ${params.dedupKey} --min-mapq ${params.minMapq} --barcode-input QnameStart --out-prefix ${sample} --write-threads $sthreads --genome $chroms
	samtools sort --threads $sthreads -n ${sample}.dedup.bam -o ${sample}.dedup.nsrt.bam
	grep -v Barcode ${sample}.cell_stats.tsv | awk -vUniq=${minUniq} '{if(\$4 >= Uniq) {print \$1}}' > ${sample}.uniqcount.annot
"""
}

// Extract methylation stats from deduplicated bam files
process extract {
input: 
	tuple(val(sample), path(annot), path(bam))
output: 
	tuple(val(sample), path("${sample}*CG.chroms"), emit: metGC )
	tuple(val(sample), path("${sample}*CH.chroms"), emit: metCH)
	tuple(val(sample), path("${sample}*cellInfo.txt"), emit: report)
tag "$sample"
script:
	ethreads = task.cpus
"""
	sciMET_BSBolt_extract.pl -O $sample -t $ethreads -C $annot $bam
	rm -rf $sample
"""
}

// Sort bed files containing CG methylation stats
process sortExtractCG {
input:
	tuple(val(sample), path(metGC)) // CG methylation stats
output: 
	tuple(val(sample), path("${sample}*CG.chroms.sort"), emit: metGC )
publishDir "${params.outDir}/cg_sort_cov/${sample}", pattern: "${sample}*CG.chroms.sort", enabled: params.covOut
tag "$sample"
script:
	ethreads = task.cpus
"""
	mkdir ${sample}.CG.chroms.sort
	for f in ${metGC}/*bed; do outFile=`basename \$f`; sort -k2,2n \$f | gzip > ${sample}.CG.chroms.sort/\${outFile}.gz; done
"""
}

// Sort bed files containing CH methylation stats
process sortExtractCH {
input:
	tuple(val(sample), path(metCH))
output: 
	tuple(val(sample), path("${sample}*CH.chroms.sort"), emit: metCH )
publishDir "${params.outDir}/ch_sort_cov/${sample}", pattern: "${sample}*CH.chroms.sort", enabled: params.covOut
tag "$sample"
script:
	ethreads = task.cpus
"""
	mkdir ${sample}.CH.chroms.sort
	for f in ${metCH}/*bed; do outFile=`basename \$f`; sort -k2,2n \$f | gzip > ${sample}.CH.chroms.sort/\${outFile}.gz; done
"""
}

// Concatenate files produced on a per well coordinate basis and format them for the purpose of publishing/reporting
// Only do the latter if splitBarcodeParsing = false
process merge_stats {
input:
	tuple(val(sample), path(cell_stats), path(frags), path(cell_info_stats), path(dedup_stats), path(mapping_stats), path(trimming_stats))
output:
	tuple(val(sample), path("concat/${sample}.cell_stats.csv"), path("concat/${sample}.cellInfo.csv"), emit: merged_stats)
	tuple(val(sample), path("concat/${sample}.fragment_hist.csv"), emit: fragment_hist)
	tuple(val(sample), path("concat/${sample}.dedup_stats.csv"), emit: dedup_stats)
	tuple(val(sample), path("concat/${sample}.mapping_stats.csv"), path("concat/${sample}.trimming_stats.csv"), emit: trimming_and_mapping_stats, optional: true)
tag "$sample"
label 'py_process'
publishDir "${params.outDir}/", mode:'copy'
script:
"""
	merge_stats.py --sampleName ${sample} --outDir concat
"""
}

// Generate per cell metrics that will be used for reporting and also perform thresholding
process generateMetrics {
input:
	tuple(val(sample), path(cell_stats), path(cell_info_stats), val(threshold))
	path(references)
output:
	tuple(val(sample), path("${sample}.allCells.csv") , emit: complex )
	tuple(val(sample), path("filter_cells/${sample}*.cells"), emit: annot, optional: true)
	tuple(val(sample), path("${sample}.passingCellsMapMethylStats.csv"), emit: cell_map_methyl_stats)
publishDir "${params.outDir}/report", pattern: "${sample}/csv", mode:'copy'
tag "$sample"
label 'py_process'
script:
	"""
	generateMetrics.py --sampleName ${sample} --tgmt_barcodes ${references}/tgmt.txt --topCellPercentage ${params.topCellPercentile} --threshold $threshold --cellInfo ${cell_info_stats} --cellStats ${cell_stats}\
	--minCellRatio ${params.minCellRatio} --minReads ${params.minUniqCount} --minUniqTotal ${params.minUniqTotal} --maxUniqTotal ${params.maxUniqTotal} --i5_barcodes ${references}/i5.txt --i7_barcodes ${references}/i7.txt\
	"""
}

process tssEnrich {
input:
	tuple(val(sample), path(cell_map_methyl_stats), path(bam), val(sample_well_coordinate))
	path(tssBed)
	path(backgroundBed)
output:
	tuple(val(sample), path("${sample_well_coordinate}.tss_enrich.csv"), emit: enrich)
tag "$sample"
script:
	"""
	met_tssEnrich.R --tssBed=${tssBed} --backgroundBed=${backgroundBed} --outName=${sample_well_coordinate} --annoFile=${cell_map_methyl_stats} --bamFile=${bam}
	"""
}

// Generate datapane html report for each sample
process generateSampleReport {
input:
	tuple(val(sample), path(allCells), path(methyl_stats), path(fragment_hist), path(dedup_stats), path(mapping_stats), path(trimming_stats), path(tss_enrich))
	path(references)
	val(trimming_and_mapping)
	val(reporting_only)
output:
	tuple(val(sample), path("${sample}/${sample}.report.html"), emit: report)
	tuple(val(sample), path("${sample}/csv/${sample}.passingCellSummaryStats.csv"), emit: passing_cell_methyl_stats)
	tuple(val(sample), path("${sample}/csv/${sample}.passingCellsMapMethylStats.csv"), emit: cell_map_methyl_stats)
	tuple(val(sample), path("${sample}.mergedTssEnrich.csv"), emit: merged_tss_enrich, optional: true)
	tuple(val(sample), path("${sample}/png"))
	tuple(val(sample), path("${sample}/csv"))
publishDir "${params.outDir}/report", pattern: "${sample}/*.html", mode:'copy'
publishDir "${params.outDir}/report", pattern: "${sample}/png", mode:'copy'
publishDir "${params.outDir}/report", pattern: "${sample}/csv", mode:'copy'
publishDir "${params.outDir}/concat", pattern: "${sample}.mergedTssEnrich.csv", mode:'copy'
tag "$sample"
label 'py_process'
script:
	opts = ""
	if (trimming_and_mapping) {
		opts = "--trimmingStats ${trimming_stats} --mappingStats ${mapping_stats} "
	}
	if (reporting_only) {
		opts = opts + "--reporting_only"
	}
	"""
	generateSampleReport.py --sampleName ${sample} --outDir ${sample} --tgmt_barcodes ${references}/tgmt.txt --allCells ${allCells} --passingCellsMapMethylStats ${methyl_stats} \
	--fragmentHist ${fragment_hist} --outDir ${sample} --tssEnrich ${tss_enrich} --dedupStats ${dedup_stats} $opts
	"""
}

// Generate combined report for the entire run
process combinedSampleReport {
input:
	tuple(val(libName), path(methyl_stats), path(complexity_stats), path(passing_met_stats), path("demuxMetrics.json"))
	path(references)
output:
	tuple(val(libName), path("${libName}/library.${libName}.report.html"), emit: combined_report)
	tuple(val(libName), path("${libName}/library.${libName}.combinedPassingCellStats.csv"), emit: combined_passing_stats)
	tuple(val(libName), path("${libName}/png"))
	tuple(val(libName), path("${libName}/csv"))
publishDir "${params.outDir}/library_report", pattern: "${libName}/png", mode:'copy'
publishDir "${params.outDir}/library_report", pattern: "${libName}/csv", mode:'copy'
publishDir "${params.outDir}/library_report", pattern: "${libName}/*.html", mode:'copy'
publishDir "${params.outDir}/library_report", pattern: "${libName}/*.csv", mode:'copy'
tag "$libName"
label 'py_process'
script:
	"""
	export TMPDIR=\$(mktemp -p `pwd` -d)
	combined_sample_report.py --libraryName ${libName} --outDir ${libName} --complexity_stats ${complexity_stats} --methyl_stats ${methyl_stats} \
	--passing_met_stats ${passing_met_stats} --i5_barcodes ${references}/i5.txt --i7_barcodes ${references}/i7.txt --tgmt_barcodes ${references}/tgmt.txt
	"""
}

// Generate binned matrix for each barcode
process mtxCG {
input:
	path(tiles) // Non overlapping genome bins
	tuple(val(sample), path(cells), path(metGC))
output: 
	tuple(val(sample), path("${sample}*CG.mtx"), emit: mtx, optional: true )
	tuple(val(sample), path("${sample}*CG.cov.mtx"), emit: cov, optional: true )
	tuple(val(sample), path("${sample}*CG.ratio.mtx"), emit: ratio, optional: true )
	tuple(val(sample), path("${sample}*CG.score.mtx"), emit: score, optional: true )
tag "$sample"
script:
"""
	sciMET_meth2mtx.pl -F $metGC -B $tiles -O ${sample} -C $cells
"""
}

// Generate binned matrix for each barcode
process mtxCH {
input:
	path(tiles) // Non overlapping genome bins
	tuple(val(sample), path(cells), path(metCH))
output: 
	tuple(val(sample), path("${sample}*CH.mtx"), emit: mtx, optional: true )
	tuple(val(sample), path("${sample}*CH.cov.mtx"), emit: cov, optional: true )
	tuple(val(sample), path("${sample}*CH.ratio.mtx"), emit: ratio, optional: true )
	tuple(val(sample), path("${sample}*CH.score.mtx"), emit: score, optional: true )
tag "$sample"
script:
"""
	sciMET_meth2mtx.pl -F $metCH -B $tiles -O ${sample} -C $cells
"""
}

// Generate merged CG matrix with all barcodes
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
tag "$sample" 
script:
"""
	sciMET_mergeMtx.pl ${sample}.CG.mtx $metGC
	sciMET_mergeMtx.pl ${sample}.CG.score.mtx $metGC_score
	sciMET_mergeMtx.pl ${sample}.CG.cov.mtx $metGC_cov
	sciMET_mergeMtx.pl ${sample}.CG.ratio.mtx $metGC_ratio
"""
}

// Generate merged CH matrix with all barcodes
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
tag "$sample"
script:
"""
	sciMET_mergeMtx.pl ${sample}.CH.mtx $metCH
	sciMET_mergeMtx.pl ${sample}.CH.score.mtx $metCH_score
	sciMET_mergeMtx.pl ${sample}.CH.cov.mtx $metCH_cov
	sciMET_mergeMtx.pl ${sample}.CH.ratio.mtx $metCH_ratio
"""
}

workflow {
	// Initialize ParamLogger class in lib/ParamLogger.groovy
	// Pretty prints select parameters for a run and initializes a logger
	ParamLogger.initialise(workflow, params, log)
	libJson = expandPath(params.libStructure, file("${projectDir}/references/"))
	// If it is not a reporting only run
	if (!params.reporting) {
		genome = loadGenome(file(params.genome))
		regularizeSamplesCsv(file(params.samples))
		samples = regularizeSamplesCsv.out.splitCsv(header:true, strip:true)
		// Indicates starting from bam files
		if (params.bam1Dir != null && params.bam2Dir != null) {
			// The bam input directory structure is expected to be exactly the same as the
			// output directory structure created by bsbolt in the align process
			bam1 = Channel.fromPath("$params.bam1Dir/*/*bam", checkIfExists: true)
			bam2 = Channel.fromPath("$params.bam2Dir/*/*bam", checkIfExists: true)
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
			mergeBam(bams, params.bamRgHeader)
			// mergeBam.out.bam -> [sample name.well coordinate, sample name.well coordinate.bam]
			mergeBam.out.bam.dump(tag: 'mergeBam.out.bam')
			dedup(mergeBam.out.bam, params.minUniqCount, genome.bsbolt_chrs)
			// If starting from bam files, trimming and mapping stats do not exist since those steps are not run
			trimming_and_mapping = false
		}
		else {
			// If run folder is provided, bcl-convert samplesheet will be generated(if not provided) and
			// bcl-convert will be executed to produce fastq files
			if (params.runFolder != null) {
				if (params.fastqSamplesheet == null) {
					// Generate bcl-convert samplesheet
					makeBclConvertSamplesheet(regularizeSamplesCsv.out, libJson.getParent(), libJson.getName(), 
								              file("$params.runFolder/RunInfo.xml", checkIfExists:true))
					fqSheet = makeBclConvertSamplesheet.out
				}
				else {
					fqSheet = file(params.fastqSamplesheet)
				}
				// Execute bcl-convert to generate fastq files
				bclconvert(file(params.runFolder), fqSheet, params.bclConvertParams)
				fqs = bclconvert.out.fastq.flatten().filter({ !it.name.contains("Undetermined") })
			} else if (params.fastqDir != null) {
				fqs = Channel.fromPath("$params.fastqDir/*fastq.gz", checkIfExists: true).filter({ !it.name.contains("Undetermined") })
			} else {
				throwError("Must specify either 'runFolder' or 'fastqDir'")
			}
			// fqs -> [single fastq file]
			// n channels for n fastq files
			fqs.dump(tag:'fqs')
			if (params.fastqc) {
				fastqc(fqs, params.adapters)
			}
			// sampleFqs -> [library name, identifier[either sample number or lane number], [fastq files corresponding to that library name and identifier]]
			sampleFqs = fqs.map { file ->	
				def ns = file.getName().toString().tokenize('_')
				return tuple(ns[0], ns[2].substring(1), file) //(ns[0], ns[1], ns[3].substring(1) // if it were split by i5 and lane
			}.groupTuple(by:[0,1])
			sampleFqs.dump(tag: 'sampleFqs')
			barcodeDemux(regularizeSamplesCsv.out, libJson.getParent(), libJson.getName(), sampleFqs)
			
			// barcodeDemux.out.fastq -> [sample name, [all demuxed fastq files corresponding to one sample name and identifier]]
			barcodeDemux.out.fastq.dump(tag: 'barcodeDemux.out.fastq')
			organized_demux = barcodeDemux.out.metrics.groupTuple()
			// merge all bcParser metrics together into one json file for each library
			merge_demux(organized_demux, libJson.getParent(), libJson.getName())
			
			demuxFqs = barcodeDemux.out.fastq.flatMap({it[1]}).map { file ->
				def ns = file.getName().toString().tokenize('_')
				return tuple(ns[0], ns[3].substring(0,3), file)
			}.groupTuple(by:[0,1])
			// if splitBarcodeParsing = true: pairedFqs -> [sample name.well coordinate, [demuxed fastq files from that well coordinate for that sample]]
			// if splitBarcodeParsing = false: pairedFqs -> [sample name, [demuxed fastq files for that sample]]
			pairedFqs = demuxFqs.map{tuple(it[0], it[2][0], it[2][1])}.groupTuple()
			pairedFqs.dump(tag: 'pairedFqs')
			// Run trim_galore on fastq files post demux
			trim(pairedFqs)
			
			// if splitBarcodeParsing = true: trim.out.fastq -> [sample name.well coordinate, [trimmed fastq files from that well coordinate for that sample]]
			// if splitBarcodeParsing = false: trim.out.fastq -> [sample name, [trimmed fastq files for that sample]]
			trim.out.fastq.dump(tag: 'trim.out.fastq')
			// Run bsbolt(aligner) on trimmed fastq files
			align(genome.bsbolt_index, trim.out.fastq)

			// if splitBarcodeParsing = true: align.out.bam -> [sample name.well coordinate, sample name.well coordinate.bam]
			// if splitBarcodeParsing = false: align.out.bam -> [sample name, sample name.bam]
			align.out.bam.dump(tag: 'align.out.bam')
			// Run scDedup for removing duplicate reads from bam files produced by bsbolt
			dedup(align.out.bam, params.minUniqCount, genome.bsbolt_chrs)
			
			// Trimming and mapping stats only exist when starting from fastq files or a runfolder
			trimming_and_mapping = true
		}
		annot = dedup.out.annot.filter ({ it[1].size() > 0 }).join(dedup.out.bam)
		// if splitBarcodeParsing = true: annot -> [sample name.well coordinate, sample name.well coordinate.uniqcount.annot, sample name.well coordinate.dedup.nstr.bam]
		// if splitBarcodeParsing = false: annot -> [sample name, sample name.uniqcount.annot, sample name.dedup.nstr.bam]
		annot.dump(tag: 'annot')
		
		// Extract methylation stats from deduplicated bam files
		extract(annot)
		
		// Sort bed files
		sortExtractCG(extract.out.metGC)
		// if splitBarcodeParsing = true: sortExtractCG.out.metGC -> [sample name.well coordinate, sample name.well coordinate.CG.chroms.sort]
		// if splitBarcodeParsing = false: sortExtractCG.out.metGC -> [sample name, sample name.CG.chroms.sort]
		sortExtractCG.out.metGC.dump(tag: 'sortExtractCG.out.metGC')
		
		// Sort bed files
		sortExtractCH(extract.out.metCH)
		// if splitBarcodeParsing = true: sortExtractCH.out.metGC -> [sample name.well coordinate, sample name.well coordinate.CH.chroms.sort]
		// if splitBarcodeParsing = false: sortExtractCH.out.metGC -> [sample name, sample name.CH.chroms.sort]
		sortExtractCH.out.metCH.dump(tag: 'sortExtractCH.out.metCH')
		
		// Collect all cell_stats.tsv files for a sample
		complexTabs = dedup.out.complexity.map({it[1]}).map { file ->
			def ns = file.getName().toString().tokenize('.').get(0)
			return tuple(ns, file)
		}.groupTuple()
		// Collect all fragment_hist.tsv files for a sample
		fragTabs = dedup.out.fragHist.map({it[1]}).map { file ->
			def ns = file.getName().toString().tokenize('.').get(0)
			return tuple(ns, file)
		}.groupTuple()
		readsFrags = complexTabs.join(fragTabs)
		// readFrags -> [sample name, [list of cell_stats.tsv files], [list of fragment_hist.tsv files]]
		readsFrags.dump(tag: 'readFrags')
		sampleExtractStats = extract.out.report.map({it[1]}).map { file ->
			def ns = file.getName().toString().tokenize('.')
			return tuple(ns.get(0), file)
		}.groupTuple()
		// sampleExtractStats -> [sample name, [list of cellInfo.txt files]]
		sampleExtractStats.dump(tag: 'sampleExtractStats')
		dedupStats = dedup.out.stats.map({it[1]}).map { file ->
				def ns = file.getName().toString().tokenize('.')
				return tuple(ns.get(0), file)
		}.groupTuple()
		dedupStats.dump(tag: 'dedupStats')
		threshold = samples.map{[it.sample, toIntOr0(it.threshold)]}
		
		if (trimming_and_mapping) {
			// Collect all bsbolt log files for a sample
			sampleMapStats = align.out.log.map({it[1]}).map { file ->
				def ns = file.getName().toString().tokenize('.')
				return tuple(ns.get(0), file)
			}.groupTuple()
			// sampleMapStats -> [sample name, [list of bsbolt.log files]]
			sampleMapStats.dump(tag: 'sampleMapStats')
			// Collect all trim_galore report files for a sample
			sampleTrimStats = trim.out.report.map({it[1]}).map { file ->
				def ns = file.getName().toString().tokenize('.')
				return tuple(ns.get(0), file)
			}.groupTuple()
			// sampleTrimStats -> [sample name, [list of trimming_report_R2.txt files]]
			sampleTrimStats.dump(tag: 'sampleTrimStats')
			report_stats = readsFrags.join(sampleExtractStats).join(dedupStats).join(sampleMapStats).join(sampleTrimStats)
			// report_stats -> [sample name, [list of cell_stats.tsv files], [list of fragment_hist.tsv files], [list of cellInfo.txt files], [list of bsbolt.log files], [list of read2 trimming report files]]
			report_stats.dump(tag: 'report_stats')
			demuxMetrics = merge_demux.out.barcode_metrics
			// Merge per well coordinate files and format them to feed into reporting
			// if splitBarcodeParsing is false only do the latter
			// Having separate merge_stats process facilitates rerunning reporting by relying on only the outputs of this process
			merge_stats(report_stats)
			
			merged_report_stats = merge_stats.out.merged_stats
			merged_report_stats = merged_report_stats.join(threshold)
			// Generate metrics that are necessary to produce reports
			generateMetrics(merged_report_stats, libJson.getParent())

			report_input = generateMetrics.out.complex.join(generateMetrics.out.cell_map_methyl_stats).join(merge_stats.out.fragment_hist).join(merge_stats.out.dedup_stats).join(merge_stats.out.trimming_and_mapping_stats)
		}
		else {
			report_stats = readsFrags.join(sampleExtractStats).join(dedupStats)
			// Empty lists indicate absence of trimming and mapping report/log files
			report_stats = report_stats.map {
				tuple(it[0], it[1], it[2], it[3], it[4], [], [])
			}
			// report_stats -> [sample name, [list of cell_stats.tsv files], [list of fragment_hist.tsv files], [list of cellInfo.txt files], [], []]
			report_stats.dump(tag: 'report_stats')
			demuxMetrics = samples.map{tuple(it.libName, [])}
			// Merge per well coordinate files and format them to feed into reporting
			merge_stats(report_stats)
			
			merged_report_stats = merge_stats.out.merged_stats
			merged_report_stats = merged_report_stats.join(threshold)
			// Generate metrics that are necessary to produce reports
			generateMetrics(merged_report_stats, libJson.getParent())
			
			report_input = generateMetrics.out.complex.join(generateMetrics.out.cell_map_methyl_stats).join(merge_stats.out.fragment_hist).join(merge_stats.out.dedup_stats)
			// Empty lists indicate absence of trimming and mapping report/log files
			report_input = report_input.map {
				tuple(it[0], it[1], it[2], it[3], it[4], [], [])
			}
		}
		dedupTmp = dedup.out.bam.map({it[1]}).map {file ->
			def ns = file.getName().toString().tokenize('.')
			return tuple(ns[0], file, ns[0] + '.' + ns[1])
		}
		if (params.runTssEnrich) {
			tssEnrichIn = generateMetrics.out.cell_map_methyl_stats.cross(dedupTmp).map{
				tuple(it[0][0], it[0][1], it[1][1], it[1][2])
			}
			tssEnrichIn.dump(tag: 'tssEnrichIn')
			tssEnrich(tssEnrichIn, genome.tssWin, genome.backgroundWin)
			tssEnrichBySample = tssEnrich.out.enrich.filter ({ it[1].size() > 1 }).map({it[1]}).map { file ->
				def ns = file.getName().toString().tokenize('.')
				return tuple(ns.get(0), file)
			}.groupTuple()
		}
		else {
			tssEnrichBySample = samples.map{tuple(it.sample, [])}
		}
		report_input = report_input.join(tssEnrichBySample)
		tssEnrichBySample.dump(tag: 'tssEnrichBySample')
		// Generate per sample html report
		generateSampleReport(report_input, libJson.getParent(), trimming_and_mapping, false)

		sample_name_lib_name = samples.map{tuple(
			it.sample, it.libName
		)}
		allMetStats = sample_name_lib_name.join(generateSampleReport.out.cell_map_methyl_stats).groupTuple(by: 1)
		allMetStats = allMetStats.map {
			tuple(it[1], it[2])
		}
		// allMetStats -> [libName, [passingCellsMapMethylStats.csv files for all samples]]
		allMetStats.dump(tag: 'allMetStats')
		allComplexStats = sample_name_lib_name.join(generateMetrics.out.complex).groupTuple(by: 1)
		allComplexStats = allComplexStats.map {
			tuple(it[1], it[2])
		}
		// allComplexStats -> [libName, [allCells.csv files for all samples]]
		allComplexStats.dump(tag: 'allComplexStats')
		passing_cell_met_stats = sample_name_lib_name.join(generateSampleReport.out.passing_cell_methyl_stats).groupTuple(by: 1)
		passing_cell_met_stats = passing_cell_met_stats.map {
			tuple(it[1], it[2])
		}
		// passing_cell_met_stats -> [libName, [passingCellSummaryStats.csv files for all samples]]
		passing_cell_met_stats.dump(tag: 'passing_cell_met_stats')
		library_stats = allMetStats.join(allComplexStats).join(passing_cell_met_stats).join(demuxMetrics)
		// Generate reports for all samples from all libraries
		combinedSampleReport(library_stats, libJson.getParent())
		
		//Different logic required for tokenization to generate input necessary for feeding into matrix generation processes when we split in bcParser vs when we don't
		if (params.splitBarcodeParsing) {
			// Get sample name.well coordinate, filter cells output for a well coordinate and join it with sorted bed files for that well coordinate
			cgPass = generateMetrics.out.annot.filter ({ it[1].size() > 0 }).flatMap({it[1]}).map {file ->
				def ns = file.getName().toString().tokenize('.')
				return tuple(ns[0] + '.' + ns[1], file)
			}.groupTuple().join(sortExtractCG.out.metGC)
			chPass = generateMetrics.out.annot.filter ({ it[1].size() > 0 }).flatMap({it[1]}).map {file ->
				def ns = file.getName().toString().tokenize('.')
				return tuple(ns[0] + '.' + ns[1], file)
			}.groupTuple().join(sortExtractCH.out.metCH)
		}
		else {
			cgPassTmp = generateMetrics.out.annot.filter ({ it[1].size() > 0 }).flatMap({it[1]}).map {file ->
					def ns = file.getName().toString().tokenize('.')
					return tuple(ns[0], file)
			}
			// Need to do a cross here since cgPassTmp will contain per well coordinate files, however sortExtract.out will be on a per sample basis since we did not split in bc_parser
			cgPass = sortExtractCG.out.metGC.cross(cgPassTmp)
			// Crossing will result in a nested list which is unpacked in the map statement
			cgPass = cgPass.map {
				def ns = it[1][1].getName().toString().tokenize('.')
				return tuple(ns[0] + '.' + ns[1], it[1][1], it[0][1])
			}
			chPassTmp = generateMetrics.out.annot.filter ({ it[1].size() > 0 }).flatMap({it[1]}).map {file ->
					def ns = file.getName().toString().tokenize('.')
					return tuple(ns[0], file)
			}
			chPass = sortExtractCH.out.metCH.cross(chPassTmp)
			chPass = chPass.map {
				def ns = it[1][1].getName().toString().tokenize('.')
				return tuple(ns[0] + '.' + ns[1], it[1][1], it[0][1])
			}
		}
        if (params.matrixGenerationCG) {
			// cgPass -> [sample name.well coordinate, [barcodes for passing cells for that well coordinate], sorted bed file] 
			cgPass.dump(tag: 'cgPass')
			mtxCG(genome.genomeTiles,cgPass)
			
			// All of the channels necessary for merging the matrix files rely on getting all matrix file which were generated on a per well coordinate basis grouped together based on the sample name
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
		}
		
		if (params.matrixGenerationCH) {
			// chPass -> [sample name.well coordinate, [barcodes for passing cells for that well coordinate], sorted bed file] 
			chPass.dump(tag: 'chPass')
			mtxCH(genome.genomeTilesCh,chPass)
			
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
    }
	// run only reporting
	else {
		if (params.resultDir == null) {
			ParamLogger.throwError(log, "Must specify --resultDir when running reporting only")
		}
		regularizeSamplesCsv(file(params.samples))
		samples = regularizeSamplesCsv.out.splitCsv(header:true, strip:true)
		// Since trimming and mapping stats might or might not be in the output directory, to avoid
		// confusion and reduce dependencies we assume that there aren't trimming and mapping stats
		// Users can refer to original report for trimming and mapping stats
		metric_input = samples.map{tuple(
			it.sample,
			file("${params.resultDir}/concat/${it.sample}.cell_stats.csv", checkIfExists:true),
			file("${params.resultDir}/concat/${it.sample}.cellInfo.csv", checkIfExists:true),
			toIntOr0(it.threshold)
		)}
		frag_hist = samples.map{tuple(
			it.sample,
			file("${params.resultDir}/concat/${it.sample}.fragment_hist.csv", checkIfExists:true)
		)}
		dedup_stats = samples.map{tuple(
			it.sample,
			file("${params.resultDir}/concat/${it.sample}.dedup_stats.csv", checkIfExists:true)
		)}
		tss_enrich = samples.map{tuple(
			it.sample,
			file("${params.resultDir}/concat/${it.sample}.mergedTssEnrich.csv")
		)}
		tss_enrich = tss_enrich.map {
			if (it[1].exists()) {
				return tuple(it[0], it[1])
			}
			else {
				return tuple(it[0], [])
			}
		}
		// Generate metrics that are necessary to produce reports
		generateMetrics(metric_input, libJson.getParent())
		
		report_input = generateMetrics.out.complex.join(generateMetrics.out.cell_map_methyl_stats).join(frag_hist).join(dedup_stats)
		// Empty lists are for trimming and mapping stats
		report_input = report_input.map {
			tuple(it[0], it[1], it[2], it[3], it[4], [], [])
		}
		report_input = report_input.join(tss_enrich)
		generateSampleReport(report_input, libJson.getParent(), false, true)

		sample_name_lib_name = samples.map{tuple(
			it.sample, it.libName
		)}
		// Since bcParser metrics might or might not be in the output directory, to avoid confusion
		// and reduce dependencies we assume that there aren't bc_parser metrics
		// Users can refer to original report for bc_parser metrics
		// Empty list is for bc_parser metrics
		demuxMetrics = samples.map{tuple(
			it.libName,
			[]
		)}.unique()
		allMetStats = sample_name_lib_name.join(generateSampleReport.out.cell_map_methyl_stats).groupTuple(by: 1)
		allMetStats = allMetStats.map {
			tuple(it[1], it[2])
		}
		allComplexStats = sample_name_lib_name.join(generateMetrics.out.complex).groupTuple(by: 1)
		allComplexStats = allComplexStats.map {
			tuple(it[1], it[2])
		}
		passing_cell_met_stats = sample_name_lib_name.join(generateSampleReport.out.passing_cell_methyl_stats).groupTuple(by: 1)
		passing_cell_met_stats = passing_cell_met_stats.map {
			tuple(it[1], it[2])
		}
		library_stats = allMetStats.join(allComplexStats).join(passing_cell_met_stats).join(demuxMetrics)
		// Generate reports for all samples from all libraries
		combinedSampleReport(library_stats, libJson.getParent())
	}
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
		     "Duration": "$workflow.duration",
			 "Exit status": "$workflow.exitStatus",
			 "Error message": "$workflow.errorMessage"
            ]
		   ]
	def params_data = ["Parameters": [:]]
	for (p in params) {
		if (!p.key.contains('-')) {
			params_data."Parameters".put("$p.key", "$p.value")
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
