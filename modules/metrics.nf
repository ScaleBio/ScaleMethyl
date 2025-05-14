/**
* Module to generate metrics, sample reports, and library reports
*
* Process:
*     GenerateMetrics
*     CombinedGenerateMetrics
*     TssEnrich
*/

// Generate per cell metrics that will be used for reporting and also perform thresholding
process GenerateMetrics {
input:
	tuple val(sample), path(cellStats), path(cellInfoStats), val(threshold), path(mappingStats), path(trimmingStats)
	path(references)
    val(libStructJsonFileName) 
	val(trimmingAndMapping)
output:
	tuple val(sample), path("${sample}.allCells.csv") , emit: allCells
    tuple val(sample), path("${sample}/csv/${sample}.passingCellSummaryStats.csv"), emit: passingCellMethylStats
    
        tuple val(sample), path("${sample}.allCells.csv")
	tuple val(sample), path("${sample}.tss_enrich.csv"), emit: mergedTssEnrich, optional: true
	tuple val(sample), path("${sample}/csv"), emit: csvFolder
tag "$sample"
label 'process_single'

publishDir file(params.outDir) / "metrics_for_reporting", pattern: "${sample}.tss_enrich.csv", mode:'copy'
publishDir { outDir }, pattern: "${sample}/csv", mode:'copy'
publishDir file(params.outDir) / "samples", pattern: "${sample}.allCells.csv", mode:'copy'
script:
    libStruct = "$references/$libStructJsonFileName"
    opts = ""
	if (trimmingAndMapping) {
		opts = "--trimming_stats ${trimmingStats} --mapping_stats ${mappingStats} "
	}
    outDir = file(params.outDir) / "report" / "sample_reports"
	"""
	generate_metrics.py --sample_name ${sample} --top_cell_percentage ${params.topCellPercentile} --threshold $threshold --cell_info ${cellInfoStats} --cell_stats ${cellStats} \
	--min_cell_ratio ${params.minCellRatio} --min_reads ${params.minUniqCount} --min_uniq_total ${params.minUniqTotal} --max_uniq_total ${params.maxUniqTotal} --library_structure_json $libStruct \
     ${opts}
	"""
}

// Generate combined report for the entire run
process CombinedGenerateMetrics {
input:
	tuple val(libName), path(allCells), path(passingCellStats), path("demuxMetrics.json")
	path(references)
        val(libStructJsonFileName)
    val(workflowVersion)
output:
	tuple val(libName), path("${libName}/${libName}.combinedPassingCellStats.csv"), emit: combinedPassingStats
	tuple val(libName), path("${libName}/csv"), emit: csvFolder
publishDir file(params.outDir) / "report" / "library_report", pattern: "${libName}/*.csv", mode:'copy'
publishDir file(params.outDir) / "report" / "library_report", pattern: "${libName}/csv", mode:'copy'
tag "$libName"
label 'process_single'
script:
    libStruct = "$references/$libStructJsonFileName"
	"""
	
	generate_combined_metrics.py --library_name ${libName} --out_dir ${libName} --all_cells ${allCells} \
	--passing_cell_stats ${passingCellStats} --library_structure_json $libStruct --workflow_version $workflowVersion
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


workflow METRICS {

take:
    mergedReportStats // Stats for all samples to aid metric generation
    fragmentHist // Fragment histogram stats from scDedup
    dedupStats // Deduplication stats from scDedup
    demuxMetrics // Demultiplexing metrics
    trimmingAndMappingStats // Trimming and mapping stats from trim_galore and bsbolt
    tssEnrichBySample //TSS enrichment stats
    dedupBam // Deduplicated bam files from scDedup
    libJson // Library structure definition file
    trimmingAndMapping // Flag to indicate if trimming and mapping stats are available
    reporting // Flag to indicate if reporting is enabled
    genome // Reference genome
    samples // sample and libName combination

main:
    // Generate metrics that are necessary to produce reports
    if(trimmingAndMapping) {
        mergedReportStats=mergedReportStats.join(trimmingAndMappingStats)
    } else {
        mergedReportStats=mergedReportStats.map {
            tuple(it[0],it[1],it[2],it[3], [], [])
        }
    }
	GenerateMetrics(mergedReportStats, libJson.getParent(), libJson.getName(), trimmingAndMapping)

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
    samplesBarcodeAndLibrary = samples.map{tuple(
        it.sample, it.libName, it.barcodes
    )}
    reportInput = reportInput.join(samplesBarcodeAndLibrary)
    tssEnrichBySample.dump(tag: 'tssEnrichBySample')

    sampleNameLibName = samples.map{tuple(
        it.sample, it.libName
    )}
    allCellsByLibrary = sampleNameLibName.join(GenerateMetrics.out.allCells).groupTuple(by: 1)
    allCellsByLibrary = allCellsByLibrary.map {
        tuple(it[1], it[2])
    }
    // allCellsByLibrary -> [libName, [allCells.csv files for all samples]]
    allCellsByLibrary.dump(tag: 'allCellsByLibrary')
    passingCellMetStats = sampleNameLibName.join(GenerateMetrics.out.passingCellMethylStats).groupTuple(by: 1)
    passingCellMetStats = passingCellMetStats.map {
        tuple(it[1], it[2])
    }
    // passingCellMetStats -> [libName, [passingCellSummaryStats.csv files for all samples]]
    passingCellMetStats.dump(tag: 'passingCellMetStats')
    libraryStats = allCellsByLibrary.join(passingCellMetStats).join(demuxMetrics)

    CombinedGenerateMetrics(libraryStats, libJson.getParent(), libJson.getName(), workflow.manifest.version)
    
emit:
    allCells = GenerateMetrics.out.allCells // Information about all cell barcodes
    reportInput = reportInput // Input for the report
    demuxMetrics = mergedReportStats // bcParser metrics on a per library basis
    passingCellMethylStats = GenerateMetrics.out.passingCellMethylStats // Passing cell methylation stats
    csvFolder = GenerateMetrics.out.csvFolder // Folder containing csv files
    combinedCsvFolder = CombinedGenerateMetrics.out.csvFolder // Folder containing csv files
    libraryStats = libraryStats // Library stats
}
