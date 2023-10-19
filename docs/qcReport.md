# ScaleMethyl QC Reports

## Sample Report
The files called _<SampleName>.report.html_ contain the summary report for a single sample, i.e. all or a subset of TN5 wells from a ScaleMethyl library. It shows read, cell and barcode level summary metrics and plots for library and sample QC.

### Passing Cell Stats
**Note**: All numbers in this table depend on the number of cells called and hence on the _Unique Reads Threshold_. Check the threshold indicated on the _Knee plot_.

**Number of Passing Cells**: The number of cell barcodes passing the _Unique Reads Threshold_\
**Median Total Reads**: The median number of total reads per passing cell barcode\
**Median Unique Reads**: The median number of unique reads per passing cell barcode\
**Median of Total Covered C's**: The median number of total covered cytosines per passing cell barcode\
**Median of Percent of Unique over Total Reads(%)**: The median percentage of unique reads versus unique over total reads\
**Median of Covered CG's**: The median number of covered cytosines in the CG context per passing cell barcode\
**Median of Covered CH's**: The median number of covered cytosines in the CH context per passing cell barcode\
**Median of CG Methylation Percent(%)**: The median percentage of covered cytosines in the CG context\
**Median of CH Methylation Percent(%)**: The median percentage of covered cytosines in the CH context

### Barcodes Tab
The heatmap on the left shows the number of reads with each TN5 barcode. \
The heatmap on the right shows the number of passing cells with each TN5 barcode.

## Library Report
The file called _library.<LibraryName>.report.html_ contains the summary report at the library level, i.e. sequencing, barcode and demultiplexing information for all samples processed in one run of the ScaleMethyl kit.

### Barcode Read Status
This table gives the barcode matching statistics for the full library. Reads that fail barcode matching are not assigned to any sample and are hence not included in any downstream analysis or metrics.

**Pass**: Reads for which all expected barcodes were found \
**Error**: Reads which were filtered because at least one barcode could not be found or matched against the expected sequences (whitelist). These reads are excluded from all further analysis. For more information about the different types of errors, refer to [the page on detailed barcode stats](detailedBarcodeStats.md)

### Unique Reads Threshold for Cell Calling

Unique reads threshold for each sample in the library

## Barcodes Tab
In the first row, the heatmap on the left shows the number of reads with each i5 barcode, and the heatmap on the right shows the number of passing cells with each i5 barcode. \
In the second row, the heatmap on the left shows the number of reads with each i7 barcode, and the heatmap on the right shows the number of passing cells with each i7 barcode.
