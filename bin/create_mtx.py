#!/usr/bin/env python
"""
Create matrix of features defined by a BED file and met calls from parquet files
"""
import duckdb
import gzip
import numpy as np
from scipy.io import mmwrite
from scipy.sparse import csc_array
import argparse
from pathlib import Path
from reporting.reporting import Utils


def score_mtx(met_rate:np.ndarray):
    """
    Calculate score matrix from methylation rate matrix

    Args:
        met_rate: Methylation rate matrix
    """
    cell_met_rate = np.mean(met_rate, axis=0)
    diff = met_rate - cell_met_rate
    score = np.where(diff > 0, diff/(1 - cell_met_rate), diff/cell_met_rate)
    return score


def create_mtx(met_calls:Path, bedfile:Path, passing_bc:list, sample:str, context:str):
    """
    Creates matrices of methylation counts for regions defined by a BED file.

    Using the methylation calls from a Parquet file and a BED file binning regions
    of the genome, create a matrix of methylation counts (met.mtx) for each region
    and cell barcode, also create a matrix of coverage (cov.mtx) for each region
    and cell barcode. Matrices are written in dense Matrix Market array format.

    Args:
        met_calls: Path to the Parquet file containing methylation calls.
        bedfile: Path to the BED file defining regions of interest.
        passing_bc: Passing cell barcodes for this well coordinate.
        sample: Used for constructing the output file names.
        context: Either CG or CH, for CG a score matrix is created, for CH a methylation rate matrix is created.
    """

    # read bed file
    bed = duckdb.sql(f"""
    SELECT column0 AS chr, column1 AS start, column2 AS stop, chr || '_' || start || '_' || stop as region
    FROM read_csv(
        '{bedfile}',
        header = false,
        delim = '\t'
        );
    """)
    
    num_regions = duckdb.sql("SELECT COUNT(*) FROM bed").fetchall()[0][0]
    met = np.zeros((num_regions, len(passing_bc)), dtype=int)
    cov = np.zeros((num_regions, len(passing_bc)), dtype=int)
    
    for idx, bc in enumerate(passing_bc):
        # get methylation calls for current barcode
        calls = duckdb.sql(f"""
        SELECT chr, pos, met
        FROM read_parquet('{met_calls}')
        WHERE barcode = '{bc}';
        """)

        # get non-zero elements of binned matrix
        cov_bins = duckdb.sql(f"""
        SELECT bed.region, COUNT(*) as count
        FROM calls
        INNER JOIN bed
        USING (chr)
        WHERE calls.pos >= bed.start AND calls.pos < bed.stop
        GROUP BY region
        """)

        # get methylated sums in binned matrix
        met_bins = duckdb.sql(f"""
        SELECT bed.region, COUNT(*) as count
        FROM calls
        INNER JOIN bed
        USING (chr)
        WHERE (calls.pos >= bed.start AND calls.pos < bed.stop) AND (calls.met SIMILAR TO '[XYZ]')
        GROUP BY region
        """)

        # get all elements of binned matrix as a numpy array
        cov_cell = duckdb.sql(f"""
        SELECT COALESCE(cov_bins.count, 0) as count
        FROM bed
        LEFT JOIN cov_bins
        USING (region)
        ORDER BY chr, start;
        """).fetchnumpy()['count']
        
        # update coverage matrix
        cov[:, idx] = cov_cell

        # get all elements of binned matrix as a numpy array
        met_cell = duckdb.sql(f"""
        SELECT COALESCE(met_bins.count, 0) as count
        FROM bed
        LEFT JOIN met_bins
        USING (region)
        ORDER BY chr, start;
        """).fetchnumpy()['count']
        
        # update coverage matrix
        met[:, idx] = met_cell
    
    met_rate = np.where(cov != 0, met / cov, 0)
    if context == "CG":
        with gzip.open(f"{sample}.{context}.score.mtx.gz", 'wb') as f:
            mmwrite(f, csc_array(score_mtx(met_rate)))
    else:
        # Use percent methylated for CH context
        with gzip.open(f"{sample}.{context}.mtx.gz", 'wb') as f:
            mmwrite(f, csc_array(met_rate))

    # Write column and row names
    with open(f"{sample}.barcodes.tsv", 'w') as f:
        for bc in passing_bc:
            f.write(f"{bc}\n")

    duckdb.sql(f"""
    COPY (
        SELECT region
        FROM bed
        ORDER BY chr, start
    ) TO '{sample}.{context}.features.tsv' (HEADER false, DELIMITER '\t');
    """)

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
