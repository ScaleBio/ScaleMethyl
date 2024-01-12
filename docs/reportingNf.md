# Reporting-only workflow

A separate subworkflow, `reporting`, is included to run the downstream steps of the workflow, after alignment and quantification (after extract and matrix generation).\
This includes:
* Cell filtering (unique reads threshold)
* Sample metric and report generation

This sub-workflow can be used to regenerate reports and filtered metric outputs after adjusting the threshold in sample.csv or using`--minUniqTotal` or `--maxUniqTotal` without re-running upstream processes.

## Usage
The workflow can be started by running `nextflow run /PATH/TO/ScaleRna -profile ... --genome ... --samples samples.csv --reporting --resultDir <OutDir from previous pipeline run>

Where `samples.csv` includes a `threshold` column, see [samples.csv](docs/samplesCsv.md) for details.

## Inputs
The `reporting` workflow will read the output in `report` and `library_report` from the previous pipeline run specified in `resultDir`. Specify a new `outDir`

## Outputs
The `reporting` workflow produces:
* A `report` and `library_report` folders with all html, csv and png for cells with the new threshold in the new outDir destination. It is recommended that you do not use the same folder as the original run, as you will need to reference your originals for the run summary stats. 
  * Cell metrics (`samples/<sample>.<libName>.allCells.csv`, `<sample>.passingCellMapMethylStats.csv`, `<sample>.passingCellSummaryStats.csv`)
  * New per-sample QC reports (`report/<sample>/<sample>.<libName>.report.html`)
    - And `.csv` metric files
  * New combined sample library reports (`library.<libName>.report.html`)

