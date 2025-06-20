# ScaleMethyl QC Reports

## Sample Report
The files called _<SampleName>.report.html_ contain the summary report for a single sample, i.e. all or a subset of TN5 wells from a ScaleMethyl library. It shows read, cell and barcode level summary metrics and plots for library and sample QC.

### Read Stats
**Total Reads**: Number of reads with valid barcodes allocated to the sample\
**Reads Passing Trimming**: Percent of **Total Reads** that pass adapter trimming filters\
**Reads Passing Mapping**: Percent of **Total Reads** that can be mapped to the genome\
**Reads Per Passing Cell**: **Total Reads** divided by the **Number of Passing Cells**
**Saturation**: The overall sequencing saturation level of this sample. _Saturation_ is defined as `1 - (UniqueReads / TotalReads)`

### Cell Statistics
**Note**: All numbers in this table depend on the number of cells called and hence on the _Unique Reads Threshold_. Check the threshold indicated on the _Barcode Rank Plot_.

**Number of Passing Cells**: The number of cell barcodes passing the _Unique Reads Threshold_\
**Reads in Passing Cells**: The percent of total reads belonging to called cells\
**Median Unique Reads**: The median number of unique reads per cell\
**Median CG's Covered**: The median number of covered cytosines in the `CG` context per cell\
**Median CH's Covered**: The median number of covered cytosines in the `CH` context per cell\
**Median CG Methylation**: The median percentage of methylated cytosines in the `CG` sequence context\
**Median CH Methylation**: The median percentage of methylated cytosines in the `CA`, `CT` or `CC` sequence context\
**Median TSS Enrichment**: The median transcription start site enrichment per cells

### Plots
#### Barcode Rank plot
This log/log plot shows the number of unique reads for each cell-barcode, sorted from high to low. The threshold used for cell-calling is indicated by a dashed line.

#### Reads Per Cell
This shows the complexity per cell, i.e. unique reads per cell-barcode, as an absolute number and relative to total reads.The color indicates called cells vs. background barcodes

### Cell Stats Tab
**Fragment Length Histogram**: Computed from the genomic distance of aligned paired-end reads\
**TSS Enrichment**: Ratio of coverage +/-1kb around Transcription Start sites vs. control regions\
**Mitochondrial Reads**: Reads mapping to the mitochondrial chromosome

### Barcodes Tab
The heatmaps on the left show the number of reads with each TN5 barcode. \
The heatmaps on the right shows the number of passing cells with each TN5 barcode.

## Library Report
The file called _library.<LibraryName>.report.html_ contains the summary report at the library level, i.e. sequencing, barcode and demultiplexing information for all samples processed in one run of the ScaleMethyl kit.

### Barcode Metrics
This table gives the barcode matching statistics for the full library. Reads that fail barcode matching are not assigned to any sample and are hence not included in any downstream analysis or metrics.

**Pass**: Reads for which all expected barcodes were found \
**Error**: Reads which were filtered because at least one barcode could not be found or matched against the expected sequences (whitelist). These reads are excluded from all further analysis. For more information about the different types of errors, refer to [the page on detailed barcode stats](detailedBarcodeStats.md)

### Cells Called
Cell-barcodes passing the unique read threshold for each sample

## Barcodes Tab
In the first row, the heatmap on the left shows the number of reads with each i5 barcode, and the heatmap on the right shows the number of passing cells with each i5 barcode. \
In the second row, the heatmap on the left shows the number of reads with each i7 barcode, and the heatmap on the right shows the number of passing cells with each i7 barcode.

## UMAP tab
Reports include a default umap on the 20,000 most variable windows calculated in the MATRIX_GENERATION step with nearest neighbor clustering with knn = min(100, # of passing cells - 1).