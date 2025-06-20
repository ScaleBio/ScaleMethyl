nextflow.enable.dsl=2
import groovy.json.JsonSlurper
import groovy.json.JsonOutput

include { INPUT_READS } from './modules/input_reads.nf'
include { ALIGNMENT } from './modules/alignment.nf'
include { INPUT_BAM_READS } from './modules/input_bam_reads.nf'
include { METRICS } from './modules/metrics.nf'
include { REPORTING } from './modules/reporting.nf'
include { MATRIX_GENERATION } from './modules/matrix_generation.nf'
include { DEDUP_AND_EXTRACT } from './modules/dedup_and_extract.nf'

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
	if(params.aligner=="bsbolt") {
		genome.bsbolt_index = expandPath(genome.bsbolt_index, baseDir)
	} else if(params.aligner=="bwa-meth") {
		genome.bwa_index = expandPath(genome.bwa_index, baseDir)
	} else if(params.aligner=="parabricks") {
		genome.parabricks_index = expandPath(genome.parabricks_index, baseDir)
	} else {
		ParamLogger.throwError("Invalid aligner specified in config: ${params.aligner}")
	}
	genome.filter_chrs = expandPath(genome.filter_chrs, baseDir)
	if(genome.genomeTileBlackList) {
		genome.genomeTileBlackList = expandPath(genome.genomeTileBlackList, baseDir)
	}
	if(genome.genomeTiles) {
		if(params.windowTileSizeCG) {
			log.warn("genomeTiles and windowTileSizeCG are both set. Ignoring windowTileSizeCG.")
			params.windowTileSizeCG = null
		}
		genome.genomeTiles = expandPath(genome.genomeTiles, baseDir)
	}
	if(!genome.genomeTiles && !params.windowTileSizeCG) {
		ParamLogger.throwError("No genomeTiles in genome.json or windowTileSizeCG in config set")
	}
	if(genome.genomeTilesCH) {
		if(params.windowTileSizeCH) {
			log.warn("genomeTilesCH and windowTileSizeCH are both set. Ignoring windowTileSizeCH.")
			params.windowTileSizeCH = null
		}
		genome.genomeTilesCH = expandPath(genome.genomeTilesCH, baseDir)
	}
	if(!genome.genomeTilesCH && !params.windowTileSizeCH) {
		ParamLogger.throwError("No genomeTilesCH in genome.json or windowTileSizeCH in config set")
	}
	if(genome.chromSizes) {
		genome.chromSizes = expandPath(genome.chromSizes, baseDir)
	}
	genome.tssWin = expandPath(genome.tssWin, baseDir)
	genome.backgroundWin = expandPath(genome.backgroundWin, baseDir)
	return genome
}

def updateGenomeBedFileCG(genome, bedFilePath) {
	genome.genomeTiles = bedFilePath
	return genome
}

def updateGenomeBedFileCH(genome, bedFilePath) {
	genome.genomeTilesCh = bedFilePath
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
process RegularizeSamplesCsv {
input: 
	path("samples.in.csv")
	val(libStructName)
	path(libStructDir)
output: 
	path("samples.csv")
publishDir file(params.outDir), mode: 'copy'
label 'process_single'
cache 'deep'
script:
	libStruct = "${libStructDir}/${libStructName}"
	opts=""
	if (params.splitFastq) {
		opts += "--splitFastq"
	}
	if (params.merged) {
		opts += " --merge"
	}
"""
	regularize_samples_csv.py samples.in.csv --libraryStruct $libStruct $opts > samples.csv
"""
}

process CreateWindowsBedCG {
input:
	path(chromSizes)
	path(genomeWindowBlackList)
	val(minimumWindowSize)
output:
	path("genome_tiles_CG.bed")
publishDir file(params.outDir) / "samples" / "genome_bin_matrix", mode: 'copy'
script:
if(genomeWindowBlackList)
"""
		create_windows_bed.py --genomeTiles $params.windowTileSizeCG --chromSizes $chromSizes --output genome_tiles_CG.bed --minimumWindowSize $minimumWindowSize --blacklist $genomeWindowBlackList
"""
else
"""
		create_windows_bed.py --genomeTiles $params.windowTileSizeCG --chromSizes $chromSizes --output genome_tiles_CG.bed --minimumWindowSize $minimumWindowSize
"""
}

process CreateWindowsBedCH {
input:
	path(chromSizes)
	path(genomeWindowBlackList)
	val(minimumWindowSize)
output:
	path("genome_tiles_CH.bed")
publishDir file(params.outDir) / "samples" / "genome_bin_matrix", mode: 'copy'
script:
if(genomeWindowBlackList)
"""
		create_windows_bed.py --genomeTiles $params.windowTileSizeCH --chromSizes $chromSizes --output genome_tiles_CH.bed --minimumWindowSize $minimumWindowSize --blacklist $genomeWindowBlackList
"""
else
"""
		create_windows_bed.py --genomeTiles $params.windowTileSizeCH --chromSizes $chromSizes --output genome_tiles_CH.bed --minimumWindowSize $minimumWindowSize
"""
}

process CreateChromSizes {
input:
	path(fastaFile)
output:
	path("chrom.sizes")
publishDir file(params.outDir) / "samples" / "genome_bin_matrix", mode: 'copy'
script:
"""
	create_chrom_sizes.py --fastaFile $fastaFile --output chrom.sizes
"""
}

process CombineBarcodeCSVs {
input:
	tuple val(sample), val(threshold), path(cellInfoFiles, stageAs:'cellInfo?'), path(cellStatFiles, stageAs:'cellStat?'), path(tssEnrichFiles, stageAs:'tssEnrich?')
output:
	tuple val(sample), val(threshold), path("${sample}.cellInfo.csv"), path("${sample}.cell_stats.csv"), path("${sample}.tss_enrich.csv")
publishDir file(params.outDir) / "metrics_for_reporting", mode: 'copy'
script:
"""
	combine_barcode_csv.py --inFiles cellInfo* --outFile ${sample}.cellInfo.csv
	combine_barcode_csv.py --inFiles cellStat* --outFile ${sample}.cell_stats.csv --header
	combine_barcode_csv.py --inFiles tssEnrich* --outFile ${sample}.tss_enrich.csv --header
"""
}

process CombineSummaryCSVs {
input:
	tuple val(sample), val(threshold), path(dedupFiles, stageAs:'dedupFile?'), path(fragFiles, stageAs:'fragFile?')
output:
	tuple val(sample), val(threshold), path("${sample}.dedup_stats.csv"), path("${sample}.fragment_hist.csv")
publishDir file(params.outDir) / "metrics_for_reporting", mode: 'copy'
script:
"""
	combine_summary_csv.py --inDedupFiles dedupFile* --outDedupFile ${sample}.dedup_stats.csv --inFragFiles fragFile* --outFragFile ${sample}.fragment_hist.csv
"""
}

// Concatenate files produced on a per well coordinate basis and format them for the purpose of publishing/reporting
// Only do the latter if splitFastq = false
process MergeStats {
input:
	tuple val(sample), path(cellStats), path(frags), path(cellInfoStats), path(dedupStats), path(mappingStats), path(trimmingStats)
output:
	tuple val(sample), path("metrics_for_reporting/${sample}.cell_stats.csv"), path("metrics_for_reporting/${sample}.cellInfo.csv"), emit: mergedStats
	tuple val(sample), path("metrics_for_reporting/${sample}.fragment_hist.csv"), emit: fragmentHist
	tuple val(sample), path("metrics_for_reporting/${sample}.dedup_stats.csv"), emit: dedupStats
	tuple val(sample), path("metrics_for_reporting/${sample}.mapping_stats.csv"), path("metrics_for_reporting/${sample}.trimming_stats.csv"), emit: trimmingAndMappingStats, optional: true
tag "$sample"
label 'process_single'
publishDir file(params.outDir), mode:'copy'
script:
"""
	merge_stats.py --sampleName ${sample} --outDir metrics_for_reporting --aligner ${params.aligner}
"""
}

workflow {
	// Initialize ParamLogger class in lib/ParamLogger.groovy
	// Pretty prints select parameters for a run and initializes a logger
	ParamLogger.initialise(workflow, params, log)
	libJson = expandPath(params.libStructure, file("${projectDir}/references/"))
	RegularizeSamplesCsv(file(params.samples, checkIfExists:true),libJson.getName(),libJson.getParent())
	samples = RegularizeSamplesCsv.out.splitCsv(header:true, strip:true)
	
	genome = loadGenome(file(params.genome))
	// chrom.size only needs to be created if a windowTileSize is set
	if(params.windowTileSizeCG || params.windowTileSizeCH) {
		// Uses aligner reference fasta to create a chrom.sizes file for matrix generation
		if(params.aligner == "bsbolt") {
			faFile = genome.bsbolt_index.resolve("BSB_ref.fa")
		} else if(params.aligner == "bwa-meth") {
			faFile = genome.bwa_index.resolve(genome.ref_fasta)
		} else if(params.aligner == "parabricks") {
			faFile = genome.parabricks_index.resolve(genome.ref_fasta)
		} else {
			ParamLogger.throwError("Invalid aligner specified in config: ${params.aligner}")
		}
		chromSizeFile = null
		if(genome.chromSizes) {
			chromSizeFile = genome.chromSizes
		}
		else {
			chromSizeFile = CreateChromSizes(faFile)
		}
		
		if(params.windowTileSizeCG) {
			if(genome.genomeWindowBlackList) {
				bedFilePathCG = CreateWindowsBedCG(chromSizeFile,genome.genomeWindowBlackList,params.minimumWindowSize)
			} else {
				bedFilePathCG = CreateWindowsBedCG(chromSizeFile,[],params.minimumWindowSize)
			}
			
			genome = updateGenomeBedFileCG(genome, bedFilePathCG)
		}
		if(params.windowTileSizeCH && params.calculateCH) {
			if(genome.genomeWindowBlackList) {
				bedFilePathCH = CreateWindowsBedCH(chromSizeFile,genome.genomeWindowBlackList,params.minimumWindowSize)
			} else {
				bedFilePathCH = CreateWindowsBedCH(chromSizeFile,[],params.minimumWindowSize)
			}
			genome = updateGenomeBedFileCH(genome, bedFilePathCH)
		}
	}
	
	allCells = null
	metCG = null
	metCH = null
	trimmingAndMapping = null
	metrics = null

	// If Matrix Checkpointing is enabled or reporting only is checked, start from matrix generation
	if (!params.startPostExtraction && !params.reportingOnly) {
		// Indicates starting from bam files
		if (params.bam1Dir != null && params.bam2Dir != null) {
			// Module that merges per TN5 bam files in two aligner output directories
			INPUT_BAM_READS(params.bam1Dir, params.bam2Dir)
			dedupBamInput = INPUT_BAM_READS.out.mergedBam
			// If starting from bam files, trimming and mapping stats do not exist since those steps are not run
			trimmingAndMapping = false
		}
		else {
			// Module that deals with input reads demux, qc, trimming and mapping
			INPUT_READS(samples, RegularizeSamplesCsv.out, libJson, params.fastqDir, params.runFolder)
			ALIGNMENT(INPUT_READS.out.trimFastq, genome)
			// Aligned bam files to be used as input for deduplication and methylation
			dedupBamInput = ALIGNMENT.out.alignedBam
			// Trimming and mapping stats only exist when starting from fastq files or a runfolder
			trimmingAndMapping = true
		}		
		
		// Module that deduplicates aligned bam files and extracts methylation stats
		DEDUP_AND_EXTRACT(dedupBamInput, genome)
		threshold = samples.map{[it.sample, toIntOr0(it.threshold)]}

		if (trimmingAndMapping) {
			// Collect all aligner log files for a sample
			sampleMapStats = ALIGNMENT.out.alignLog.map({it[1]}).map { file ->
				def ns = file.getName().toString().tokenize('.')
				return tuple(ns.get(0), file)
			}.groupTuple()
			// sampleMapStats -> [sample name, [list of aligner log files]]
			sampleMapStats.dump(tag: 'sampleMapStats')
			// Collect all trim_galore report files for a sample
			sampleTrimStats = INPUT_READS.out.trimLog.map({it[1]}).map { file ->
				def ns = file.getName().toString().tokenize('.')
				return tuple(ns.get(0), file)
			}.groupTuple()
			// sampleTrimStats -> [sample name, [list of trimming_report_R2.txt files]]
			sampleTrimStats.dump(tag: 'sampleTrimStats')
			reportStats = DEDUP_AND_EXTRACT.out.readsFrags.join(DEDUP_AND_EXTRACT.out.sampleExtractStats).join(DEDUP_AND_EXTRACT.out.dedupStats).join(sampleMapStats).join(sampleTrimStats)
			// reportStats -> [sample name, [list of cellStats.tsv files], [list of fragment_hist.tsv files], [list of cellInfo.txt files], [list of bsbolt.log files], [list of read2 trimming report files]]
			reportStats.dump(tag: 'reportStats')
			demuxMetrics = INPUT_READS.out.mergedBcParserMetrics
			// Merge per well coordinate files and format them to feed into reporting
			// if splitFastq is false only do the latter
			// Having separate MergeStats process facilitates rerunning reporting by relying on only the outputs of this process
			MergeStats(reportStats)

			trimmingAndMappingStats = MergeStats.out.trimmingAndMappingStats
		}
		else {
			reportStats = DEDUP_AND_EXTRACT.out.readsFrags.join(DEDUP_AND_EXTRACT.out.sampleExtractStats).join(DEDUP_AND_EXTRACT.out.dedupStats)
			// Empty lists indicate absence of trimming and mapping report/log files
			reportStats = reportStats.map {
				tuple(it[0], it[1], it[2], it[3], it[4], [], [])
			}
			// reportStats -> [sample name, [list of cell_stats.tsv files], [list of fragment_hist.tsv files], [list of cellInfo.txt files], [], []]
			reportStats.dump(tag: 'reportStats')
			demuxMetrics = samples.map{tuple(it.libName, [])}.unique()
			// Merge per well coordinate files and format them to feed into reporting
			MergeStats(reportStats)

			trimmingAndMappingStats = null
		}
		
		
		mergedReportStats = MergeStats.out.mergedStats
		mergedReportStats = mergedReportStats.join(threshold)
		
		// Generate metrics from merged stats and use those metrics for generating html reports
		metrics = METRICS(mergedReportStats, MergeStats.out.fragmentHist, MergeStats.out.dedupStats, demuxMetrics, trimmingAndMappingStats, null, DEDUP_AND_EXTRACT.out.dedupBam, libJson, trimmingAndMapping, false, genome, samples)
			
		allCells = METRICS.out.allCells.groupTuple()
		metCG = DEDUP_AND_EXTRACT.out.metCG.groupTuple()
		metCH = DEDUP_AND_EXTRACT.out.metCH.groupTuple()
	} else {
		// Reporting code before post-alignment code
		// Set up variables for reporting
		metricInput=null
		fragHist=null
		dedupStats=null
		tssEnrich=null
		// if merged, use samplesheet to merge metrics
		if(params.merged) {
			// allCellInfos: [sample name, threshold, cellInfo file, cell_stats file, tss_enrich file]
			allCellInfos = samples
				.map{tuple(it.sample,toIntOr0(it.threshold), file("${it.resultDir}/metrics_for_reporting/${it.sample}.cellInfo.csv", checkIfExists:true),file("${it.resultDir}/metrics_for_reporting/${it.sample}.cell_stats.csv", checkIfExists:true),file("${it.resultDir}/metrics_for_reporting/${it.sample}.tss_enrich.csv"))}
				.groupTuple()
				.map{tuple(it[0], it[1][0], it[2], it[3], it[4])}
			// Combine cellInfo, cell_stats, and tss_enrich files
			CombineBarcodeCSVs(allCellInfos)

			// allDedupStats: [sample name, threshold, dedup_stats file, fragment_hist file]
			allDedupStats = samples
				.map{tuple(it.sample,toIntOr0(it.threshold), file("${it.resultDir}/metrics_for_reporting/${it.sample}.dedup_stats.csv", checkIfExists:true), file("${it.resultDir}/metrics_for_reporting/${it.sample}.fragment_hist.csv", checkIfExists:true))}
				.groupTuple()
				.map{tuple(it[0], it[1][0], it[2], it[3])}
			// Combine dedup_stats and fragment_hist files
			CombineSummaryCSVs(allDedupStats)

			// Set reporting variables to combination function outputs
			metricInput=CombineBarcodeCSVs.out.map{tuple(it[0],it[3],it[2],it[1])}
			fragHist=CombineSummaryCSVs.out.map{tuple(it[0],it[3])}
			dedupStats=CombineSummaryCSVs.out.map{tuple(it[0],it[2])}
			tssEnrich=CombineBarcodeCSVs.out.map{tuple(it[0],it[4])}

		// Else, use the output directory to merge metrics
		} else {
			// Set output directory to the previous output directory if it exists, otherwise use the current output directory
			if(params.previousOutDir != null) {
				outDir = params.previousOutDir
			} else {
				outDir = params.outDir
			}

			// metricInput: [sample name, cell_stats file, cellInfo file, threshold]
			metricInput = samples.map{tuple(
				it.sample,
				file("${outDir}/metrics_for_reporting/${it.sample}.cell_stats.csv"),
				file("${outDir}/metrics_for_reporting/${it.sample}.cellInfo.csv"),
				toIntOr0(it.threshold)
			)}

			// fragHist: [sample name, fragment_hist file]
			fragHist = samples.map{tuple(
				it.sample,
				file("${outDir}/metrics_for_reporting/${it.sample}.fragment_hist.csv")
			)}

			// dedupStats: [sample name, dedup_stats file]
			dedupStats = samples.map{tuple(
				it.sample,
				file("${outDir}/metrics_for_reporting/${it.sample}.dedup_stats.csv")
			)}

			// tssEnrich: [sample name, tss_enrich file]
			tssEnrich = samples.map{tuple(
				it.sample,
				file("${outDir}/metrics_for_reporting/${it.sample}.tss_enrich.csv")
			)}
		}
		
		// For the case where original analysis was run with runTssEnrich = false
		tssEnrich = tssEnrich.map {
			if (it[1].exists()) {
				return tuple(it[0], it[1])
			}
			else {
				return tuple(it[0], [])
			}
		}

		// Since bcParser metrics might or might not be in the output directory, to avoid confusion
		// and reduce dependencies we assume that there aren't bc_parser metrics
		// Users can refer to original report for bc_parser metrics
		// Empty list is for bc_parser metrics
		// bcParser metrics will not be there in the case of starting from bam files
		demuxMetrics = samples.map{tuple(
			it.libName,
			[]
		)}.unique()

		// Generate metrics from merged stats and use those metrics for generating html reports
		metrics = METRICS(metricInput, fragHist, dedupStats, demuxMetrics, null, tssEnrich, null, libJson, false, true, genome, samples)
		trimmingAndMapping = false
		// If Matrix Checkpointing is enabled, start from matrix generation
		if(params.startPostExtraction) {
			// if using previous output directory, use the previous allCells.csv and met_ files
			if(params.previousOutDir != null) {
				// allCells: [sample name, list of allCells files]
				allCells = samples
									.map{tuple(it.sample, file("${params.previousOutDir}/samples/${it.sample}.allCells.csv", checkIfExists:true))}
									.groupTuple()
				// contexts: [CG, CH]
				contexts = ""
				if(params.calculateCH) {
					contexts = Channel.of('CG', 'CH')
				} else {
					contexts = Channel.of('CG')
				}

				// subsamples: [sample name, list of parquet files, context]
									// [sample name, result directory, allCells file]
				subsamples = samples.map{tuple(it.sample, params.previousOutDir, file("${params.previousOutDir}/samples/${it.sample}.allCells.csv", checkIfExists:true).splitCsv(header:true, strip:true))}
									// [sample name, result directory, wells from allCells]
									.map{tuple(it[0], it[1], it[2].tgmt_well.unique())}
									// [sample name, result directory, well list]
									.flatMap{sampleName, dir, wells -> wells.collect{well -> tuple(sampleName,dir,well)}}
									// [sample name.well, parquet file path beginning]
									.map{tuple("${it[0]}.${it[2]}","${it[1]}/matrix_generation/${it[0]}.${it[2]}.met_")}
									// [sample name.well, parquet file path with context]
									.combine(contexts)
									// [sample name.well, parquet file path, context]
									.map{tuple(it[0], file("${it[1]}${it[2]}.parquet"), it[2])}
									// Remove non-existent files
									.filter{it[1].exists()}
									// Group by sample name and context [sample name.well, list of parquet files, context]
									.groupTuple(by:[0,2])

				// branch subsamples based on context for matrix generation
				subsamples.branch {
					metCG: it[2]=='CG'
					metCH: it[2]=='CH'
				}.set { mets }
				metCG = mets.metCG.map{tuple(it[0],it[1])}
				metCH = mets.metCH.map{tuple(it[0],it[1])}
			
			// if merged option is selected, use the sample sheet to find the files to merge
			} else if(params.merged) {
				// allCells: [sample name, list of allCells files]
				allCells = samples.map{tuple(it.sample, file("${it.resultDir}/samples/${it.sample}.allCells.csv", checkIfExists:true))}
								  .groupTuple()

				// contexts: [CG, CH]
				contexts = Channel.of('CG', 'CH')

				// subsamples: [sample name, list of parquet files, context]
				subsamples = samples.map{tuple(it.sample, it.resultDir, file("${it.resultDir}/samples/${it.sample}.allCells.csv", checkIfExists:true).splitCsv(header:true, strip:true))}
									.map{tuple(it[0], it[1], it[2].tgmt_well.unique())}
									.flatMap{sampleName, dir, wells -> wells.collect{well -> tuple(sampleName,dir,well)}}
									.map{tuple("${it[0]}.${it[2]}","${it[1]}/matrix_generation/${it[0]}.${it[2]}.met_")}
									.combine(contexts)
									.map{tuple(it[0], file("${it[1]}${it[2]}.parquet"), it[2])}
									.filter{it[1].exists()}
									.groupTuple(by:[0,2])

				// branch subsamples based on context for matrix generation					
				subsamples.branch {
					metCG: it[2]=='CG'
					metCH: it[2]=='CH'
				}.set { mets }
				metCG = mets.metCG.map{tuple(it[0],it[1])}
				metCH = mets.metCH.map{tuple(it[0],it[1])}
			}
		}
	}

	CGmtx = Channel.of([])
	CGmtxBarcodes = Channel.of([])
	
	if(!params.reportingOnly) {
		// Generate binned matrices
		MATRIX_GENERATION(allCells,metCG,metCH, genome)

		// CG matrix per sample
		CGmtx = MATRIX_GENERATION.out.CGmtx
		CGmtxBarcodes = MATRIX_GENERATION.out.CGmtxBarcodes
	}
	REPORTING(metrics.reportInput, libJson, trimmingAndMapping, params.reportingOnly, metrics.passingCellMethylStats, metrics.csvFolder, metrics.combinedCsvFolder, metrics.libraryStats, CGmtx, CGmtxBarcodes)
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
	def paramsData = ["Parameters": [:]]
	for (p in params) {
		if (!p.key.contains('-')) {
			paramsData."Parameters".put("$p.key", "$p.value")
		}
    }
	def referenceData = ["Reference Genome": [:]]
	for (p in genome) {
		referenceData."Reference Genome".put("$p.key", "$p.value")
	}
	def manifestData = ["Manifest": [:]]
	for (p in workflow.manifest.getProperties()) {
		p = p.toString()
		def splitStr = p.split("=")
		if (splitStr[0].equals("name") || splitStr[0].equals("version")) {
			manifestData."Manifest".put(splitStr[0], splitStr[1])
		}
	}
	def jsonStr = JsonOutput.toJson(data+paramsData+referenceData+manifestData)
	def jsonBeauty = JsonOutput.prettyPrint(jsonStr)
	workflowInfo = file("$params.outDir/workflow_info.json")
	workflowInfo.write(jsonBeauty)
}
