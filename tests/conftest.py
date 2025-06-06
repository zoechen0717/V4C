import os
import pytest
import pandas as pd
import numpy as np
from pathlib import Path

@pytest.fixture
def test_data_dir():
    """Return the path to the test data directory."""
    return Path(__file__).parent / "data"

@pytest.fixture
def sample_mcool_file(test_data_dir):
    """Return the path to a sample mcool file."""
    return str(test_data_dir / "sample.mcool")

@pytest.fixture
def sample_coords():
    """Return sample genomic coordinates."""
    return "chr17:45878152-46000000"

@pytest.fixture
def sample_genes():
    """Return sample gene names."""
    return "MAPT,CRHR1"

@pytest.fixture
def sample_resolutions():
    """Return sample resolutions."""
    return [5000, 10000]

@pytest.fixture
def sample_output_file(tmp_path):
    """Return a temporary output file path."""
    return str(tmp_path / "test_output.tsv")

@pytest.fixture
def sample_dataframe():
    """Return a sample DataFrame for testing."""
    data = {
        "mcool": ["sample1.mcool", "sample2.mcool"],
        "res": [5000, 10000],
        "chrom": ["chr17", "chr17"],
        "start": [45878152, 45878152],
        "end": [46000000, 46000000],
        "contact_freqs": [
            [0.1, 0.2, 0.3, 0.2, 0.1],
            [0.15, 0.25, 0.35, 0.25, 0.15]
        ]
    }
    return pd.DataFrame(data)

@pytest.fixture
def mock_cooler_matrix():
    """Return a mock Hi-C contact matrix."""
    return np.array([
        [1.0, 0.5, 0.3, 0.2, 0.1],
        [0.5, 1.0, 0.4, 0.3, 0.2],
        [0.3, 0.4, 1.0, 0.4, 0.3],
        [0.2, 0.3, 0.4, 1.0, 0.5],
        [0.1, 0.2, 0.3, 0.5, 1.0]
    ]) 