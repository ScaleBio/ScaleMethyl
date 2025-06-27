/**
* Module to generate metrics, sample reports, and library reports
*
* Process:
*     GenerateSampleReport
*     CombinedSampleReport
*/

// Generate datapane html report for each sample
process GenerateSampleReport {
input:
	tuple val(sample), path("allCells.csv"), path(fragmentHist), path(dedupStats), path(mappingStats), path(trimmingStats), path(tssEnrich), val(libraryName), val(sampleBarcodes), path(passingCellStats), path(csvFolder), path(CGmtx), path(CGmtxBarcodes) 
        path(references)
        val(libStructJsonFileName)
	val(trimmingAndMapping)
	val(reportingOnly)
    val(workflowVersion)
	
output:
	tuple val(sample), path("${sample}/${sample}.report.html"), emit: report
	
	tuple val(sample), path("${sample}/png")
	tuple val(sample), path("${sample}/csv/${sample}.report_clusters.tsv"), optional: true
	tuple val(sample), path("${sample}.tss_enrich.csv"), emit: mergedTssEnrich, optional: true

publishDir { outDir }, pattern: "${sample}/*.html", mode:'copy'
publishDir { outDir }, pattern: "${sample}/png", mode:'copy'
publishDir { outDir }, pattern: "${sample}/csv/${sample}.report_clusters.tsv", mode:'copy'
publishDir file(params.outDir) / "metrics_for_reporting", pattern: "${sample}.tss_enrich.csv", mode:'copy'
tag "$sample"
label 'process_medium'
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
	options = CGmtx ? "--mtxFile ${CGmtx} --mtxBarcodes ${CGmtxBarcodes}" : ""
	"""
	export DATAPANE_CDN_BASE="https://d3j2ibc7el7s1k.cloudfront.net/v0.17.0"
	export TMPDIR=\$(mktemp -p `pwd` -d)
	generate_sample_report.py --sample_name ${sample} --out_dir ${sample} --all_cells allCells.csv --library_structure_json $libStruct \
	--fragment_hist ${fragmentHist} --tss_enrich ${tssEnrich} --dedup_stats ${dedupStats} $opts --library_name $libraryName --sample_barcodes "$sampleBarcodes" --workflow_version $workflowVersion \
    --passing_cell_stats ${passingCellStats} --csv_folder ${csvFolder} ${options}
	"""
}

// Generate combined report for the entire run
process CombinedSampleReport {
input:
	tuple val(libName), path(allCells), path(passingCellStats), path("demuxMetrics.json"), path(combinedCsvFolder)
	path(references)
        val(libStructJsonFileName)
    val(workflowVersion)
output:
	tuple val(libName), path("${libName}/${libName}.report.html"), emit: combinedReport
	tuple val(libName), path("${libName}/png")
publishDir file(params.outDir) / "report" / "library_report", pattern: "${libName}/*.html", mode:'copy'
publishDir file(params.outDir) / "report" / "library_report", pattern: "${libName}/png", mode:'copy'
tag "$libName"
label 'process_single'
script:
    libStruct = "$references/$libStructJsonFileName"
	"""
	export DATAPANE_CDN_BASE="https://d3j2ibc7el7s1k.cloudfront.net/v0.17.0"
	export TMPDIR=\$(mktemp -p `pwd` -d)
	combined_sample_report.py --library_name ${libName} --out_dir ${libName} --all_cells ${allCells} \
	--passing_cell_stats ${passingCellStats} --library_structure_json $libStruct --workflow_version $workflowVersion \
    --csv_folder ${combinedCsvFolder}
	"""
}

workflow REPORTING {

take:
    reportInput // Input for the sample report from GenerateMetrics
    libJson // Library structure definition file
    trimmingAndMapping // Flag to indicate if trimming and mapping stats are available
    reportingOnly // Flag to indicate if reportingOnly is set
    passingCellMethylStats // Passing cell methylation stats
    csvFolder // Folder containing csv files for sample plate plots
    combinedCsvFolder // Folder containing csv files for combined plate plots
    libraryStats // Combined library stat csvs from GenerateCombinedMetrics
	CGmtx // CG methylation matrix
	CGmtxBarcodes // Barcodes for CG mtx
main:
	if(!params.reportingOnly && params.windowMatrixOut) {
    	reportInput = reportInput.join(passingCellMethylStats).join(csvFolder).join(CGmtx).join(CGmtxBarcodes)
	} else {
		// CGmtx = CGmtxBarcodes = [] if params.windowMatrixOut is false or params.reportingOnly is true
		reportInput = reportInput.join(passingCellMethylStats).join(csvFolder).merge(CGmtx).merge(CGmtxBarcodes)
	}
    GenerateSampleReport(reportInput, libJson.getParent(), libJson.getName(), trimmingAndMapping, reportingOnly, workflow.manifest.version)

    libraryStats = libraryStats.join(combinedCsvFolder)
    // Generate reports for all samples from all libraries
    CombinedSampleReport(libraryStats, libJson.getParent(), libJson.getName(), workflow.manifest.version,)
}
