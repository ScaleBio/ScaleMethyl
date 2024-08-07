#!/usr/bin/env python
"""
Create matrix of features defined by a BED file and met calls from parquet files
"""
import duckdb
import numpy as np
import pandas as pd
from scipy.io import mmwrite
from scipy.sparse import coo_array
import argparse
from pathlib import Path
from reporting.reporting import Utils


def score_mtx(met_rate:np.ndarray):
    """
    Calculate score matrix from methylation rate matrix

    Args:
        met_rate: Methylation rate matrix
    """
    cell_met_rate = np.nanmean(met_rate, axis=0)
    diff = met_rate - cell_met_rate
    score = np.where(diff > 0, diff/(1 - cell_met_rate), diff/cell_met_rate)
    return score


def create_mtx(met_calls:Path, bedfile:Path, passing_bc:list, sample:str, context:str):
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
    bed = duckdb.sql(f"""
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
    """)
    
    calls = duckdb.sql(f"""
    SELECT barcode, chr, pos, methylated, unmethylated
    FROM '{met_calls}';
    """)
    
    # get non-sparse entries of the matrix
    perc_methylated_df = duckdb.sql(f"""
    SELECT
        barcode,
        region,
        SUM(methylated) / (SUM(methylated) + SUM(unmethylated)) as perc_methylated
    FROM calls m, bed b
    WHERE m.pos >= b.start
    AND m.pos < b.stop
    AND m.chr = b.chr
    GROUP BY barcode, region
    """).df()
    
    # filter out barcodes that are not passing
    # these are typically less than 1% of positions
    perc_methylated_df = perc_methylated_df[perc_methylated_df['barcode'].isin(passing_bc)]
    
    # regions will be rows in matrix
    regions = duckdb.sql(f"""
    SELECT region
    FROM bed;
    """).df()['region']
    region_mapping = pd.Series(regions.index, index=regions)
    row = np.array(perc_methylated_df['region'].map(region_mapping))

    # passing barcodes will be columns in matrix
    col = np.array(perc_methylated_df['barcode'].map({bc: idx for idx, bc in enumerate(passing_bc)}))
    data = perc_methylated_df['perc_methylated'].to_numpy()
    
    mtx = coo_array((data, (row, col)), shape=(len(regions), len(passing_bc))).todense()
    mask = np.ones_like(mtx, dtype=bool)
    mask[row, col] = False
    # set missing data to NaN
    mtx[mask] = np.nan

    if context == "CG":
        with open(f"{sample}.{context}.score.mtx", 'wb') as f:
            # set missing data to 0, representing average cell methylation rate
            mmwrite(f, coo_array(np.nan_to_num(score_mtx(mtx))), precision=3)
    else:
        # Use percent methylated for CH context
        with open(f"{sample}.{context}.mtx", 'wb') as f:
            # set missing data to 0
            mmwrite(f, coo_array(np.nan_to_num(mtx)), precision=3)

    # Write column and row names
    with open(f"{sample}.barcodes.tsv", 'w') as f:
        for bc in passing_bc:
            f.write(f"{bc}\n")

    regions.to_csv(f"{sample}.{context}.features.tsv", sep='\t', index=False, header=False)

def main():
    parser = argparse.ArgumentParser("Create a matrix of met calls")
    parser.add_argument("--bedfile", type=Path, required=True, help="BED file with non-overlapping regions")
    parser.add_argument("--met_calls", type=Path, help="Parquet file with met calls")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--CH", action="store_true", default=False, help="Use CH context met calls")
    parser.add_argument("--all_cells", type=Path, help="File with all cell info")
    args = parser.parse_args()

    
    create_mtx(
        met_calls=args.met_calls,
        bedfile=args.bedfile,
        passing_bc=Utils.get_passing_cells(args.all_cells, args.sample),
        sample=args.sample,
        context="CG" if not args.CH else "CH"
    )


if __name__ == "__main__":
    main()
