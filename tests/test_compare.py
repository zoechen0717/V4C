import pytest
import pandas as pd
import numpy as np
import os
from v4c import compare_v4c
from v4c.extract import V4CError, InputValidationError, FileProcessingError

def test_compare_v4c_basic(sample_dataframe, sample_output_file, tmp_path):
    """Test basic functionality of compare_v4c."""
    # Save sample data
    sample_dataframe.to_csv(sample_output_file, sep="\t", index=False)
    
    # Test comparison
    output_plot = str(tmp_path / "test_compare.png")
    compare_v4c(
        input_files=[sample_output_file],
        ylim=0.4,
        scale=True,
        output_file=output_plot
    )
    
    # Verify plot was created
    assert os.path.exists(output_plot)

def test_compare_v4c_multiple_files(sample_dataframe, sample_output_file, tmp_path):
    """Test compare_v4c with multiple input files."""
    # Create two different datasets
    df1 = sample_dataframe.copy()
    df2 = sample_dataframe.copy()
    df2["contact_freqs"] = df2["contact_freqs"].apply(lambda x: [v * 1.5 for v in x])
    
    # Save datasets
    file1 = str(tmp_path / "data1.tsv")
    file2 = str(tmp_path / "data2.tsv")
    df1.to_csv(file1, sep="\t", index=False)
    df2.to_csv(file2, sep="\t", index=False)
    
    # Test comparison
    output_plot = str(tmp_path / "test_compare.png")
    compare_v4c(
        input_files=[file1, file2],
        ylim=0.4,
        scale=True,
        output_file=output_plot
    )
    
    assert os.path.exists(output_plot)

def test_compare_v4c_invalid_inputs():
    """Test compare_v4c with invalid inputs."""
    with pytest.raises(InputValidationError):
        compare_v4c([], ylim=0.4)
    
    with pytest.raises(InputValidationError):
        compare_v4c(["nonexistent.tsv"], ylim=0.4)
    
    with pytest.raises(InputValidationError):
        compare_v4c(["test.txt"], ylim=0.4)  # Wrong file extension
    
    with pytest.raises(InputValidationError):
        compare_v4c(["test.tsv"], ylim=-1)  # Invalid ylim

def test_compare_v4c_invalid_data(sample_output_file):
    """Test compare_v4c with invalid data format."""
    # Create invalid data
    invalid_data = pd.DataFrame({
        "invalid_column": [1, 2, 3]
    })
    invalid_data.to_csv(sample_output_file, sep="\t", index=False)
    
    with pytest.raises(InputValidationError):
        compare_v4c([sample_output_file], ylim=0.4)

def test_compare_v4c_custom_parameters(sample_dataframe, sample_output_file, tmp_path):
    """Test compare_v4c with custom parameters."""
    # Save sample data
    sample_dataframe.to_csv(sample_output_file, sep="\t", index=False)
    
    # Test with custom parameters
    output_plot = str(tmp_path / "test_compare.png")
    compare_v4c(
        input_files=[sample_output_file],
        ylim=0.8,
        scale=False,
        output_file=output_plot,
        dpi=300,
        figsize=(12, 8),
        colors=["#FF0000", "#00FF00"]
    )
    
    assert os.path.exists(output_plot)

def test_compare_v4c_scaling(sample_dataframe, sample_output_file, tmp_path):
    """Test compare_v4c with and without scaling."""
    # Save sample data
    sample_dataframe.to_csv(sample_output_file, sep="\t", index=False)
    
    # Test with scaling
    output_scaled = str(tmp_path / "test_compare_scaled.png")
    compare_v4c(
        input_files=[sample_output_file],
        ylim=0.4,
        scale=True,
        output_file=output_scaled
    )
    
    # Test without scaling
    output_unscaled = str(tmp_path / "test_compare_unscaled.png")
    compare_v4c(
        input_files=[sample_output_file],
        ylim=0.4,
        scale=False,
        output_file=output_unscaled
    )
    
    assert os.path.exists(output_scaled)
    assert os.path.exists(output_unscaled)

def test_compare_v4c_empty_data(sample_output_file):
    """Test compare_v4c with empty data."""
    # Create empty data
    empty_data = pd.DataFrame(columns=["mcool", "res", "chrom", "start", "end"])
    empty_data.to_csv(sample_output_file, sep="\t", index=False)
    
    with pytest.raises(FileProcessingError):
        compare_v4c([sample_output_file], ylim=0.4)

def test_compare_v4c_file_permissions(sample_dataframe, sample_output_file, tmp_path):
    """Test compare_v4c with file permission issues."""
    sample_dataframe.to_csv(sample_output_file, sep="\t", index=False)
    
    # Try to save to a directory without write permissions
    output_plot = "/root/test_compare.png"  # This should fail on most systems
    with pytest.raises(FileProcessingError):
        compare_v4c(
            input_files=[sample_output_file],
            ylim=0.4,
            output_file=output_plot
        )
