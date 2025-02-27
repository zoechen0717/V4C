import cooler
import numpy as np
import pandas as pd
import argparse
from .utils import get_promoter_coords

def extract_v4c(mcool_files, resolutions, coords, genome=None, bed_file=None, flank=50000, balance=True, scale=True, output="extracted_data.tsv"):
    """
    Extracts Virtual 4C contact frequencies from .mcool files at specific resolutions.

    Parameters:
    - mcool_files (list): List of .mcool files.
    - resolutions (list): List of resolutions to extract.
    - coords (str): Genomic coordinates (e.g., "chr17:45878152-46000000").
    - genome (str, optional): Reference genome ("hg38" or "hg19").
    - bed_file (str, optional): Custom BED file for genomic coordinates.
    - flank (int): Number of base pairs to extend upstream and downstream.
    - balance (bool): Whether to use ICE balancing.
    - scale (bool): Whether to normalize values between 0 and 1.
    - output (str): Output file name.

    Returns:
    - Saves extracted contact frequencies to a TSV file.
    """

    chrom, start, end = coords.split(":")[0], int(coords.split(":")[1].split("-")[0]), int(coords.split(":")[1].split("-")[1])

    if genome:
        promoter_coords = get_promoter_coords(genome, chrom, start, end)
    elif bed_file:
        promoter_coords = pd.read_csv(bed_file, sep="\t", header=None, names=["chrom", "start", "end"])
    else:
        promoter_coords = [(chrom, start, end)]

    results = []

    for mcool in mcool_files:
        for res in resolutions:
            c = cooler.Cooler(f"{mcool}::/resolutions/{res}")
            for chrom, start, end in promoter_coords:
                region_start = max(0, start - flank)
                region_end = end + flank
                matrix = c.matrix(balance=balance, sparse=False).fetch((chrom, region_start, region_end))

                row_index = (start // res) - (region_start // res)
                row_values = matrix[row_index, :].tolist()

                if scale:
                    min_val, max_val = np.min(row_values), np.max(row_values)
                    row_values = [(x - min_val) / (max_val - min_val) if max_val > min_val else 0 for x in row_values]

                genomic_coords = np.arange(region_start, region_end, res)
                results.append([mcool, res, chrom, start, end] + row_values)

    df = pd.DataFrame(results)
    df.to_csv(output, sep="\t", index=False)
