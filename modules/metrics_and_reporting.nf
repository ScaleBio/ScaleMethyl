/**
* Module to generate metrics, sample reports, and library reports
*
* Process:
*     GenerateMetrics
*     TssEnrich
*     GenerateSampleReport
*     CombinedSampleReport
*/

// Generate per cell metrics that will be used for reporting and also perform thresholding
process GenerateMetrics {
input:
	tuple val(sample), path(cellStats), path(cellInfoStats), val(threshold)
	path(references)
    val(libStructJsonFileName)
output:
	tuple val(sample), path("${sample}.allCells.csv") , emit: allCells
tag "$sample"
label 'process_single'
script:
    libStruct = "$references/$libStructJsonFileName"
	"""
	generate_metrics.py --sample_name ${sample} --top_cell_percentage ${params.topCellPercentile} --threshold $threshold --cell_info ${cellInfoStats} --cell_stats ${cellStats} \
	--min_cell_ratio ${params.minCellRatio} --min_reads ${params.minUniqCount} --min_uniq_total ${params.minUniqTotal} --max_uniq_total ${params.maxUniqTotal} --library_structure_json $libStruct
	"""
}

process TssEnrich {
input:
	tuple val(sample), path(allCells), path(bam), val(sampleWellCoordinate)
	path(tssBed)
	path(backgroundBed)
output:
	tuple val(sample), path("${sampleWellCoordinate}.tss_enrich.csv"), emit: enrich
tag "$sampleWellCoordinate"
label 'process_single'
script:
	"""
	tss_enrich.py --tss_bed $tssBed --bg_bed $backgroundBed --id $sampleWellCoordinate --cells $allCells --bam $bam
	"""
}

// Generate datapane html report for each sample
process GenerateSampleReport {
input:
	tuple val(sample), path("allCells.csv"), path(fragmentHist), path(dedupStats), path(mappingStats), path(trimmingStats), path(tssEnrich) // avoid overwriting input *.allCells.csv
        path(references)
        val(libStructJsonFileName)
	val(trimmingAndMapping)
	val(reportingOnly)
output:
	tuple val(sample), path("${sample}/${sample}.report.html"), emit: report
	tuple val(sample), path("${sample}/csv/${sample}.passingCellSummaryStats.csv"), emit: passingCellMethylStats
        tuple val(sample), path("${sample}.allCells.csv")
	tuple val(sample), path("${sample}.tss_enrich.csv"), emit: mergedTssEnrich, optional: true
	tuple val(sample), path("${sample}/png")
	tuple val(sample), path("${sample}/csv")
publishDir { outDir }, pattern: "${sample}/*.html", mode:'copy'
publishDir { outDir }, pattern: "${sample}/csv", mode:'copy'
publishDir { outDir }, pattern: "${sample}/png", mode:'copy'
publishDir file(params.outDir) / "metrics_for_reporting", pattern: "${sample}.tss_enrich.csv", mode:'copy'
publishDir file(params.outDir) / "samples", pattern: "${sample}.allCells.csv", mode:'copy'
tag "$sample"
label 'process_single'
script:
	opts = ""
	if (trimmingAndMapping) {
		opts = "--trimming_stats ${trimmingStats} --mapping_stats ${mappingStats} "
	}
	if (reportingOnly) {
		opts = opts + "--reporting_only"
	}
    libStruct = "$references/$libStructJsonFileName"
    outDir = file(params.outDir) / "report" / "sample_reports"
	"""
	generate_sample_report.py --sample_name ${sample} --out_dir ${sample} --all_cells allCells.csv --library_structure_json $libStruct \
	--fragment_hist ${fragmentHist} --tss_enrich ${tssEnrich} --dedup_stats ${dedupStats} $opts
	"""
}

// Generate combined report for the entire run
process CombinedSampleReport {
input:
	tuple val(libName), path(allCells), path(passingCellStats), path("demuxMetrics.json")
	path(references)
        val(libStructJsonFileName)
output:
	tuple val(libName), path("${libName}/library.${libName}.report.html"), emit: combinedReport
	tuple val(libName), path("${libName}/library.${libName}.combinedPassingCellStats.csv"), emit: combinedPassingStats
	tuple val(libName), path("${libName}/png")
	tuple val(libName), path("${libName}/csv")
publishDir file(params.outDir) / "report" / "library_report", pattern: "${libName}/*.html", mode:'copy'
publishDir file(params.outDir) / "report" / "library_report", pattern: "${libName}/*.csv", mode:'copy'
publishDir file(params.outDir) / "report" / "library_report", pattern: "${libName}/csv", mode:'copy'
publishDir file(params.outDir) / "report" / "library_report", pattern: "${libName}/png", mode:'copy'
tag "$libName"
label 'process_single'
script:
    libStruct = "$references/$libStructJsonFileName"
	"""
	export TMPDIR=\$(mktemp -p `pwd` -d)
	combined_sample_report.py --library_name ${libName} --out_dir ${libName} --all_cells ${allCells} \
	--passing_cell_stats ${passingCellStats} --library_structure_json $libStruct
	"""
}

workflow METRICS_AND_REPORTING {

take:
    mergedReportStats // Stats for all samples to aid metric generation
    fragmentHist // Fragment histogram stats from scDedup
    dedupStats // Deduplication stats from scDedup
    trimmingAndMappingStats // Trimming and mapping stats from trim_galore and bsbolt
    demuxMetrics // bcParser metrics on a per library basis
    tssEnrichBySample //TSS enrichment stats
    dedupBam // Deduplicated bam files from scDedup
    libJson // Library structure definition file
    trimmingAndMapping // Flag to indicate if trimming and mapping stats are available
    reporting // Flag to indicate if reporting is enabled
    genome // Reference genome
    samples // sample and libName combination

main:
    // Generate metrics that are necessary to produce reports
	GenerateMetrics(mergedReportStats, libJson.getParent(), libJson.getName())

    if (trimmingAndMapping) {
        reportInput = GenerateMetrics.out.allCells.join(fragmentHist).join(dedupStats).join(trimmingAndMappingStats)
    }
    else {
        reportInput = GenerateMetrics.out.allCells.join(fragmentHist).join(dedupStats)
        // Empty lists indicate absence of trimming and mapping report/log files
        reportInput = reportInput.map {
            tuple(it[0], it[1], it[2], it[3], [], [])
        }
    }
    if (!reporting) {
        if (params.runTssEnrich) {
            dedupTmp = dedupBam.map({it[1]}).map {file ->
                def ns = file.getName().toString().tokenize('.')
                return tuple(ns[0], file, ns[0] + '.' + ns[1])
            }
            tssEnrichIn = GenerateMetrics.out.allCells.cross(dedupTmp).map{
                tuple(it[0][0], it[0][1], it[1][1], it[1][2])
            }
            tssEnrichIn.dump(tag: 'tssEnrichIn')
            TssEnrich(tssEnrichIn, genome.tssWin, genome.backgroundWin)
            tssEnrichBySample = TssEnrich.out.enrich.filter ({ it[1].size() > 1 }).map({it[1]}).map { file ->
                def ns = file.getName().toString().tokenize('.')
                return tuple(ns.get(0), file)
            }.groupTuple()
        }
        else {
            tssEnrichBySample = samples.map{tuple(it.sample, [])}
        }
    }
    reportInput = reportInput.join(tssEnrichBySample)
    tssEnrichBySample.dump(tag: 'tssEnrichBySample')
    // Generate per sample html report
    GenerateSampleReport(reportInput, libJson.getParent(), libJson.getName(), trimmingAndMapping, reporting)

    sampleNameLibName = samples.map{tuple(
        it.sample, it.libName
    )}
    allCellsByLibrary = sampleNameLibName.join(GenerateMetrics.out.allCells).groupTuple(by: 1)
    allCellsByLibrary = allCellsByLibrary.map {
        tuple(it[1], it[2])
    }
    // allCellsByLibrary -> [libName, [allCells.csv files for all samples]]
    allCellsByLibrary.dump(tag: 'allCellsByLibrary')
    passingCellMetStats = sampleNameLibName.join(GenerateSampleReport.out.passingCellMethylStats).groupTuple(by: 1)
    passingCellMetStats = passingCellMetStats.map {
        tuple(it[1], it[2])
    }
    // passingCellMetStats -> [libName, [passingCellSummaryStats.csv files for all samples]]
    passingCellMetStats.dump(tag: 'passingCellMetStats')
    libraryStats = allCellsByLibrary.join(passingCellMetStats).join(demuxMetrics)
    // Generate reports for all samples from all libraries
    CombinedSampleReport(libraryStats, libJson.getParent(), libJson.getName())

emit:
    allCells = GenerateMetrics.out.allCells // Information about all cell barcodes
}
