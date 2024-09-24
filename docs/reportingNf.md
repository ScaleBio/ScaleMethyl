# Reporting-only workflow

A separate subworkflow, `reporting`, is included to run the downstream steps of the workflow, after alignment and quantification (after extract and matrix generation).\
This includes:
* Cell filtering (unique reads threshold)
* Sample metric and report generation

This sub-workflow can be used to regenerate reports and filtered metric outputs after adjusting the threshold in sample.csv or using `--minUniqTotal` or `--maxUniqTotal` without re-running upstream processes.

## Usage
The workflow can be started by running `nextflow run /PATH/TO/ScaleRna -profile ... --genome ... --samples samples.csv --reporting --resultDir <OutDir from previous pipeline run>

Where `samples.csv` includes a `threshold` column, see documentation on [samples csv](samplesCsv.md) for details.

## Inputs
The `reporting` workflow will read the metrics from the `concat` folder in the output directory. Other than the standard samples csv and reference genome json file, to run the reporting only workflow, the `--resultDir` argument needs to be provided. This points to the output directory of the original run. For the reporting only workflow, it is recommended to provide a new output directory, otherwise the original output directory will be overwritten

## Outputs
The `reporting` workflow produces:
* A `report` and `library_report` folders with all html, csv and png for cells with the new threshold in the new outDir destination. It is recommended that you do not use the same folder as the original run, as you will need to reference your originals for the run summary stats. 
  * Cell metrics (`<sample>.allCells.csv`, `<sample>.passingCellMapMethylStats.csv`, `<sample>.passingCellSummaryStats.csv`)
  * New per-sample QC reports (`report/<sample>/<sample>.report.html`)
    - And `.csv` metric files
  * New combined sample library reports (`library.<libName>.report.html`)
The reports will not have the summary stats and barcode metrics in them, since we do not read in the files containing that information for the reporting only workflow. This is to reduce the number of dependencies, and since those numbers are already available in the original report, and do not change with subsequent workflow executions, the user can always reference the original report for those numbers.
