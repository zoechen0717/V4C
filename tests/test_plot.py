import pytest
import pandas as pd
import numpy as np
import os
from v4c import plot_v4c
from v4c.extract import V4CError, InputValidationError, FileProcessingError

def test_plot_v4c_basic(sample_dataframe, sample_output_file, tmp_path):
    """Test basic functionality of plot_v4c."""
    # Save sample data
    sample_dataframe.to_csv(sample_output_file, sep="\t", index=False)
    
    # Test plotting
    output_plot = str(tmp_path / "test_plot.png")
    plot_v4c(
        input_file=sample_output_file,
        ylim=0.4,
        flank=50000,
        output_file=output_plot
    )
    
    # Verify plot was created
    assert os.path.exists(output_plot)

def test_plot_v4c_invalid_inputs():
    """Test plot_v4c with invalid inputs."""
    with pytest.raises(InputValidationError):
        plot_v4c("nonexistent.tsv", ylim=0.4)
    
    with pytest.raises(InputValidationError):
        plot_v4c("test.txt", ylim=0.4)  # Wrong file extension
    
    with pytest.raises(InputValidationError):
        plot_v4c("test.tsv", ylim=-1)  # Invalid ylim
    
    with pytest.raises(InputValidationError):
        plot_v4c("test.tsv", flank=-1)  # Invalid flank

def test_plot_v4c_invalid_data(sample_output_file):
    """Test plot_v4c with invalid data format."""
    # Create invalid data
    invalid_data = pd.DataFrame({
        "invalid_column": [1, 2, 3]
    })
    invalid_data.to_csv(sample_output_file, sep="\t", index=False)
    
    with pytest.raises(InputValidationError):
        plot_v4c(sample_output_file, ylim=0.4)

def test_plot_v4c_custom_parameters(sample_dataframe, sample_output_file, tmp_path):
    """Test plot_v4c with custom parameters."""
    # Save sample data
    sample_dataframe.to_csv(sample_output_file, sep="\t", index=False)
    
    # Test with custom parameters
    output_plot = str(tmp_path / "test_plot.png")
    plot_v4c(
        input_file=sample_output_file,
        ylim=0.8,
        flank=100000,
        output_file=output_plot,
        dpi=300,
        figsize=(12, 8)
    )
    
    assert os.path.exists(output_plot)

def test_plot_v4c_multiple_resolutions(sample_dataframe, sample_output_file, tmp_path):
    """Test plot_v4c with multiple resolutions."""
    # Create data with multiple resolutions
    data = {
        "mcool": ["sample.mcool"] * 2,
        "res": [5000, 10000],
        "chrom": ["chr17"] * 2,
        "start": [45878152] * 2,
        "end": [46000000] * 2,
        "contact_freqs": [
            [0.1, 0.2, 0.3, 0.2, 0.1],
            [0.15, 0.25, 0.35, 0.25, 0.15]
        ]
    }
    df = pd.DataFrame(data)
    df.to_csv(sample_output_file, sep="\t", index=False)
    
    # Test plotting
    output_plot = str(tmp_path / "test_plot.png")
    plot_v4c(
        input_file=sample_output_file,
        ylim=0.4,
        output_file=output_plot
    )
    
    assert os.path.exists(output_plot)

def test_plot_v4c_empty_data(sample_output_file):
    """Test plot_v4c with empty data."""
    # Create empty data
    empty_data = pd.DataFrame(columns=["mcool", "res", "chrom", "start", "end"])
    empty_data.to_csv(sample_output_file, sep="\t", index=False)
    
    with pytest.raises(FileProcessingError):
        plot_v4c(sample_output_file, ylim=0.4)

def test_plot_v4c_invalid_ylim(sample_dataframe, sample_output_file):
    """Test plot_v4c with invalid ylim."""
    sample_dataframe.to_csv(sample_output_file, sep="\t", index=False)
    
    with pytest.raises(InputValidationError):
        plot_v4c(sample_output_file, ylim=0)  # ylim must be positive

def test_plot_v4c_file_permissions(sample_dataframe, sample_output_file, tmp_path):
    """Test plot_v4c with file permission issues."""
    sample_dataframe.to_csv(sample_output_file, sep="\t", index=False)
    
    # Try to save to a directory without write permissions
    output_plot = "/root/test_plot.png"  # This should fail on most systems
    with pytest.raises(FileProcessingError):
        plot_v4c(
            input_file=sample_output_file,
            ylim=0.4,
            output_file=output_plot
        )
