import pandas as pd
import numpy as np
import os
from typing import List

def get_promoter_coords(genome, chrom=None, start=None, end=None, gene=None):
    """
    Retrieves promoter coordinates for a given genome build.

    Parameters:
    - genome (str): "hg38" or "hg19"
    - chrom (str, optional): Chromosome
    - start (int, optional): Start position
    - end (int, optional): End position
    - gene (str, optional): Gene name

    Returns:
    - list of tuples [(chrom, start, end)]
    """
    promoter_file = f"genome/{genome}_promoters.bed"
    promoters = pd.read_csv(promoter_file, sep="\t", header=None, names=["chrom", "start", "end", "gene"])

    if gene:
        return promoters[promoters["gene"] == gene][["chrom", "start", "end"]].values.tolist()

    return promoters[(promoters["chrom"] == chrom) & (promoters["start"] <= end) & (promoters["end"] >= start)][["chrom", "start", "end"]].values.tolist()


def normalize_data(data):
    """
    Applies Min-Max normalization to scale data between 0 and 1.

    Parameters:
    - data (list or np.array): Raw contact frequency data

    Returns:
    - np.array: Normalized data
    """
    min_val, max_val = np.min(data), np.max(data)
    if max_val > min_val:
        return (data - min_val) / (max_val - min_val)
    else:
        return np.zeros_like(data)

def validate_mcool_files(mcool_files: List[str]) -> List[str]:
    """
    Validates and filters .mcool files.

    Args:
        mcool_files: List of paths to .mcool files

    Returns:
        List of valid .mcool file paths

    Raises:
        ValueError: If no valid files are found
    """
    valid_files = []
    for file in mcool_files:
        if not os.path.exists(file):
            print(f"Warning: File not found: {file}")
            continue
        if not file.endswith('.mcool'):
            print(f"Warning: File is not a .mcool file: {file}")
            continue
        valid_files.append(file)
    
    if not valid_files:
        raise ValueError("No valid .mcool files found")
    
    return valid_files