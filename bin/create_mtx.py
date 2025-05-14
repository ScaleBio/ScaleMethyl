#!/usr/bin/env python
"""
Create matrix of features defined by a BED file and met calls from parquet files
"""
import duckdb
import numpy as np
import pandas as pd
from scipy.io import mmwrite
from numpy import unique
from scipy.sparse import coo_array, csc_array
import argparse
from pathlib import Path
from reporting.reporting import Utils


def score_mtx(met_rate: csc_array):
    """
    Calculate score matrix from methylation rate matrix

    Args:
        met_rate: Methylation rate matrix
    """
    # Calculate cell methylation rate (mean per column ignoring zeros)
    colSums = np.matrix(met_rate.sum(axis=0, dtype=np.float32)).getA1()
    row, col = met_rate.nonzero()
    nonzeroEntriesPerCol = np.zeros(met_rate.shape[1])
    uniqueCols = unique(col)

    for i in uniqueCols:
        nonzeroEntriesPerCol[i] = sum(col == i)

    cell_met_rates = np.zeros(met_rate.shape[1], dtype=np.float32)
    cell_met_rates[uniqueCols] = colSums[uniqueCols] / nonzeroEntriesPerCol[uniqueCols]

    diff = np.matrix(met_rate[row, col]).getA1() - cell_met_rates[col]

    # Calculate the score for each non-zero entry
    score_data = np.where(
        diff > 0,
        diff / (1 - (cell_met_rates[col] - 1)),
        diff / (cell_met_rates[col] - 1),
    )
    return coo_array((score_data, (row, col)), shape=met_rate.shape)


def create_mtx(met_calls: list[str], bedfile: Path, passing_bc: list, sample: str, context: str):
    """
    Creates matrices of methylation rates for regions defined by a BED file.

    Using the methylation calls from a Parquet file and a BED file binning regions
    of the genome, create a score matrix for CG context methylation and a percent
    methylated matrix for CH context methylation for each region and cell barcode.
    Matrices are written in sparse Matrix Market array format.

    Args:
        met_calls: Path to the Parquet file containing methylation calls.
        bedfile: Path to the BED file defining regions of interest.
        passing_bc: Passing cell barcodes for this well coordinate.
        sample: Used for constructing the output file names.
        context: Either CG or CH, for CG a score matrix is created, for CH a methylation rate matrix is created.
    """

    # read bed file
    bed = duckdb.sql(
        f"""
    SELECT
        column0 AS chr,
        column1 AS start,
        column2 AS stop,
        chr || '_' || start || '_' || stop as region
    FROM read_csv(
        '{bedfile}',
        header = false,
        delim = '\t'
        );
    """
    )

    # regions will be rows in matrix
    regions = duckdb.sql(
        """
    SELECT region
    FROM bed;
    """
    ).df()["region"]
    region_mapping = pd.Series(regions.index, index=regions)
    # Create a mapping from barcodes to their indices
    col_mapping = {bc: idx for idx, bc in enumerate(passing_bc)}

    mtx = csc_array((len(regions), len(passing_bc)), dtype=np.float32)
    file_names = [file.name for file in met_calls]
    for idx, file in enumerate(file_names):
        # get non-sparse entries of the matrix
        # Adding 1 here to keep the values non-zero so they are not removed when crushing into a sparse matrix
        # This will be subtracted later in the score_mtx function and when outputing the CH context matrix
        result = duckdb.sql(
            f"""
        SELECT
            barcode,
            region,
            1+(SUM(methylated) / (SUM(methylated) + SUM(unmethylated))) as perc_methylated
        FROM read_parquet({file}) m, bed b
        WHERE m.pos >= b.start
        AND m.pos < b.stop
        AND m.chr = b.chr
        GROUP BY barcode, region
        """
        ).df()
        # removing barcodes that are not passing to save memory
        result = result[result["barcode"].isin(passing_bc)]
        # Check if the result DataFrame is empty
        if result.empty:
            continue  # Skip to the next file if no data

        # Map regions to their corresponding indices
        row = np.array(result["region"].map(region_mapping))
        # Map barcodes to their corresponding indices
        col = np.array(result["barcode"].map(col_mapping))
        # Create a temporary COO matrix for the current file's data
        tempmtx = csc_array(
            (result["perc_methylated"].to_numpy(dtype=np.float32), (row, col)),
            shape=(len(regions), len(passing_bc)),
        )
        # Incorporate the new data into the matrix
        mtx += tempmtx

    if context == "CG":
        with open(f"{sample}.{context}.score.mtx", "wb") as f:
            # score_mtx function will subtract the pseudo count by subtracting the average
            mmwrite(f, coo_array(score_mtx(mtx)), precision=3)
    else:
        # Use percent methylated for CH context
        with open(f"{sample}.{context}.mtx", "wb") as f:
            # Subtract 1 to get the original values
            nonzeroIndices = mtx.nonzero()
            mtx[nonzeroIndices] = mtx[nonzeroIndices] - 1
            mmwrite(f, coo_array(mtx), precision=3)
    # Write column and row names
    with open(f"{sample}.barcodes.tsv", "w") as f:
        for bc in passing_bc:
            f.write(f"{bc}\n")
    regions.to_csv(f"{sample}.{context}.features.tsv", sep="\t", index=False, header=False)


def main():
    parser = argparse.ArgumentParser("Create a matrix of met calls")
    parser.add_argument(
        "--bedfile",
        type=Path,
        required=True,
        help="BED file with non-overlapping regions",
    )
    parser.add_argument("--met_calls", type=Path, nargs="+", help="Parquet file with met calls")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--CH", action="store_true", default=False, help="Use CH context met calls")
    parser.add_argument("--all_cells", type=Path, nargs="+", help="File with all cell info")
    args = parser.parse_args()

    create_mtx(
        met_calls=args.met_calls,
        bedfile=args.bedfile,
        passing_bc=Utils.get_passing_cells(args.all_cells, args.sample),
        sample=args.sample,
        context="CG" if not args.CH else "CH",
    )


if __name__ == "__main__":
    main()
