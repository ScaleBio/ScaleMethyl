nextflow.enable.dsl=2
import groovy.json.JsonSlurper
import groovy.json.JsonOutput

include { INPUT_READS } from './modules/input_reads.nf'
include { INPUT_BAM_READS } from './modules/input_bam_reads.nf'
include { METRICS_AND_REPORTING } from './modules/metrics_and_reporting.nf'
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
process RegularizeSamplesCsv {
input: 
	path("samples.in.csv")
output: 
	path("samples.csv")
publishDir file(params.outDir), mode: 'copy'
label 'process_single'
cache 'deep'
script:
	opts=""
	if (params.splitFastq) {
		opts += "--splitFastq"
	}
"""
	regularize_samples_csv.py samples.in.csv $opts > samples.csv
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
	merge_stats.py --sampleName ${sample} --outDir metrics_for_reporting
"""
}

workflow {
	// Initialize ParamLogger class in lib/ParamLogger.groovy
	// Pretty prints select parameters for a run and initializes a logger
	ParamLogger.initialise(workflow, params, log)
	libJson = expandPath(params.libStructure, file("${projectDir}/references/"))
	RegularizeSamplesCsv(file(params.samples))
	samples = RegularizeSamplesCsv.out.splitCsv(header:true, strip:true)
	genome = loadGenome(file(params.genome))
	// If it is not a reporting only run
	if (!params.reporting) {
		// Indicates starting from bam files
		if (params.bam1Dir != null && params.bam2Dir != null) {
			// Module that merges per TN5 bam files in two bsbolt output directories
			INPUT_BAM_READS(params.bam1Dir, params.bam2Dir)
			dedupBamInput = INPUT_BAM_READS.out.mergedBam
			// If starting from bam files, trimming and mapping stats do not exist since those steps are not run
			trimmingAndMapping = false
		}
		else {
			// Module that deals with input reads demux, qc, trimming and mapping
			INPUT_READS(samples, RegularizeSamplesCsv.out, libJson, params.fastqDir, params.runFolder, genome)
			dedupBamInput = INPUT_READS.out.alignedBam						
			// Trimming and mapping stats only exist when starting from fastq files or a runfolder
			trimmingAndMapping = true
		}
		// Module that deduplicates aligned bam files and extracts methylation stats
		DEDUP_AND_EXTRACT(dedupBamInput, genome)
		threshold = samples.map{[it.sample, toIntOr0(it.threshold)]}
		
		if (trimmingAndMapping) {
			// Collect all bsbolt log files for a sample
			sampleMapStats = INPUT_READS.out.alignLog.map({it[1]}).map { file ->
				def ns = file.getName().toString().tokenize('.')
				return tuple(ns.get(0), file)
			}.groupTuple()
			// sampleMapStats -> [sample name, [list of bsbolt.log files]]
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
		METRICS_AND_REPORTING(mergedReportStats, MergeStats.out.fragmentHist, MergeStats.out.dedupStats, trimmingAndMappingStats, demuxMetrics, null, DEDUP_AND_EXTRACT.out.dedupBam, libJson, trimmingAndMapping, false, genome, samples)

		// Generate binned matrices
		MATRIX_GENERATION(METRICS_AND_REPORTING.out.allCells, DEDUP_AND_EXTRACT.out.metCG, DEDUP_AND_EXTRACT.out.metCH, genome)
    }
	// run only reporting
	else {
		// Since trimming and mapping stats might or might not be in the output directory, to avoid
		// confusion and reduce dependencies we assume that there aren't trimming and mapping stats
		// Users can refer to original report for trimming and mapping stats
		metricInput = samples.map{tuple(
			it.sample,
			file("${params.resultDir}/metrics_for_reporting/${it.sample}.cell_stats.csv", checkIfExists:true),
			file("${params.resultDir}/metrics_for_reporting/${it.sample}.cellInfo.csv", checkIfExists:true),
			toIntOr0(it.threshold)
		)}
		fragHist = samples.map{tuple(
			it.sample,
			file("${params.resultDir}/metrics_for_reporting/${it.sample}.fragment_hist.csv", checkIfExists:true)
		)}
		dedupStats = samples.map{tuple(
			it.sample,
			file("${params.resultDir}/metrics_for_reporting/${it.sample}.dedup_stats.csv", checkIfExists:true)
		)}
		tssEnrich = samples.map{tuple(
			it.sample,
			file("${params.resultDir}/metrics_for_reporting/${it.sample}.tss_enrich.csv")
		)}
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
		METRICS_AND_REPORTING(metricInput, fragHist, dedupStats, null, demuxMetrics, tssEnrich, null, libJson, false, true, genome, samples)
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
