import pandas as pd
import numpy as np

def get_promoter_coords(genome, chrom, start, end):
    """
    Retrieves promoter coordinates for a given genome build.

    Parameters:
    - genome (str): "hg38" or "hg19"
    - chrom (str): Chromosome
    - start (int): Start position
    - end (int): End position

    Returns:
    - list of tuples [(chrom, start, end)]
    """
    promoter_file = f"genome/{genome}_promoters.bed"
    promoters = pd.read_csv(promoter_file, sep="\t", header=None, names=["chrom", "start", "end", "gene"])
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

def validate_mcool_files(mcool_files):
    """
    Checks if provided .mcool files exist.

    Parameters:
    - mcool_files (list): List of .mcool file paths

    Returns:
    - list: Valid .mcool file paths
    """
    import os
    valid_files = [f for f in mcool_files if os.path.exists(f)]
    if len(valid_files) != len(mcool_files):
        missing_files = set(mcool_files) - set(valid_files)
        print(f"Warning: Missing files: {missing_files}")
    return valid_files
