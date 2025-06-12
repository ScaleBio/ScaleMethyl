# Merging Multiple Alignments

It could be useful for large kits to split up the ScaleMethyl pipeline for each plate separately to save on compute. This involves running ScaleMethyl on a separate set of fastqs from each plate with matrix generation deactivated (--windowMatrixOut false --covOut false --amethystOut false --allcOut false) and creating a special merger samplesheet (example in [samples_merge.csv](examples/samples_merge.csv) and running ScaleMethyl with --startPostExtraction  and --merged options. The merger samplesheet has the following property: a sample name for each sample across all plates and a resultDir option with the output directory from each of the smaller ScaleMethyl runs. This will also merge the reports and perform any matrix generation that's activated. Use --reportingOnly instead of --startPostExtraction  if you only wish for the combined reports.

### Example

sample | barcodes | libName | resultDir
-- | -- | -- | --
One | 1A01-2D12 | combinedLibrary | resultDirectory1
Two | 2E01-2H12 | combinedLibrary | resultDirectory2



