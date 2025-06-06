import pytest
import pandas as pd
import numpy as np
from v4c import extract_v4c
from v4c.extract import V4CError, InputValidationError, FileProcessingError
import os

def test_extract_v4c_basic(sample_mcool_file, sample_coords, sample_resolutions, sample_output_file):
    """Test basic functionality of extract_v4c."""
    extract_v4c(
        mcool_files=[sample_mcool_file],
        resolutions=sample_resolutions,
        coords=sample_coords,
        genome="hg38",
        output=sample_output_file
    )
    
    # Verify output file exists and is not empty
    assert os.path.exists(sample_output_file)
    df = pd.read_csv(sample_output_file, sep="\t")
    assert not df.empty

def test_extract_v4c_with_genes(sample_mcool_file, sample_genes, sample_resolutions, sample_output_file):
    """Test extract_v4c with gene names."""
    extract_v4c(
        mcool_files=[sample_mcool_file],
        resolutions=sample_resolutions,
        genes=sample_genes,
        genome="hg38",
        output=sample_output_file
    )
    
    df = pd.read_csv(sample_output_file, sep="\t")
    assert not df.empty
    assert "mcool" in df.columns
    assert "res" in df.columns

def test_extract_v4c_invalid_inputs():
    """Test extract_v4c with invalid inputs."""
    with pytest.raises(InputValidationError):
        extract_v4c([], [5000], coords="chr17:45878152-46000000")
    
    with pytest.raises(InputValidationError):
        extract_v4c(["test.mcool"], [], coords="chr17:45878152-46000000")
    
    with pytest.raises(InputValidationError):
        extract_v4c(["test.mcool"], [5000], genome="invalid_genome")

def test_extract_v4c_invalid_coords():
    """Test extract_v4c with invalid coordinates."""
    with pytest.raises(InputValidationError):
        extract_v4c(["test.mcool"], [5000], coords="invalid_coords")
    
    with pytest.raises(InputValidationError):
        extract_v4c(["test.mcool"], [5000], coords="chr17:46000000-45878152")  # start > end

def test_extract_v4c_missing_files():
    """Test extract_v4c with non-existent files."""
    with pytest.raises(FileProcessingError):
        extract_v4c(["nonexistent.mcool"], [5000], coords="chr17:45878152-46000000")

def test_extract_v4c_scaling(sample_mcool_file, sample_coords, sample_resolutions, sample_output_file):
    """Test extract_v4c with and without scaling."""
    # Test with scaling
    extract_v4c(
        mcool_files=[sample_mcool_file],
        resolutions=sample_resolutions,
        coords=sample_coords,
        genome="hg38",
        scale=True,
        output=sample_output_file
    )
    df_scaled = pd.read_csv(sample_output_file, sep="\t")
    
    # Test without scaling
    extract_v4c(
        mcool_files=[sample_mcool_file],
        resolutions=sample_resolutions,
        coords=sample_coords,
        genome="hg38",
        scale=False,
        output=sample_output_file
    )
    df_unscaled = pd.read_csv(sample_output_file, sep="\t")
    
    # Verify scaling effect
    assert not np.array_equal(df_scaled.iloc[:, 5:], df_unscaled.iloc[:, 5:])

def test_extract_v4c_balancing(sample_mcool_file, sample_coords, sample_resolutions, sample_output_file):
    """Test extract_v4c with and without ICE balancing."""
    # Test with balancing
    extract_v4c(
        mcool_files=[sample_mcool_file],
        resolutions=sample_resolutions,
        coords=sample_coords,
        genome="hg38",
        balance=True,
        output=sample_output_file
    )
    df_balanced = pd.read_csv(sample_output_file, sep="\t")
    
    # Test without balancing
    extract_v4c(
        mcool_files=[sample_mcool_file],
        resolutions=sample_resolutions,
        coords=sample_coords,
        genome="hg38",
        balance=False,
        output=sample_output_file
    )
    df_unbalanced = pd.read_csv(sample_output_file, sep="\t")
    
    # Verify balancing effect
    assert not np.array_equal(df_balanced.iloc[:, 5:], df_unbalanced.iloc[:, 5:])

def test_extract_v4c_multiple_files(sample_mcool_file, sample_coords, sample_resolutions, sample_output_file):
    """Test extract_v4c with multiple input files."""
    extract_v4c(
        mcool_files=[sample_mcool_file, sample_mcool_file],  # Using same file twice for testing
        resolutions=sample_resolutions,
        coords=sample_coords,
        genome="hg38",
        output=sample_output_file
    )
    
    df = pd.read_csv(sample_output_file, sep="\t")
    assert len(df["mcool"].unique()) == 2

def test_extract_v4c_output_format(sample_mcool_file, sample_coords, sample_resolutions, sample_output_file):
    """Test the format of extract_v4c output."""
    extract_v4c(
        mcool_files=[sample_mcool_file],
        resolutions=sample_resolutions,
        coords=sample_coords,
        genome="hg38",
        output=sample_output_file
    )
    
    df = pd.read_csv(sample_output_file, sep="\t")
    required_columns = ["mcool", "res", "chrom", "start", "end"]
    assert all(col in df.columns for col in required_columns)
    assert len(df.columns) > len(required_columns)  # Should have additional contact frequency columns
