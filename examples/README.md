# V4C Examples

This directory contains examples demonstrating how to use the V4C package.

## Files

### `V4C.ipynb`
Original Jupyter notebook showing the core V4C analysis logic. This notebook served as the foundation for the V4C package development.

### `api_example.py`
Comprehensive Python script demonstrating all V4C API features:
- Basic extraction using gene names
- Coordinate-based extraction
- BED file extraction
- Plotting with custom settings
- Sample comparison
- Complete workflow
- Data analysis and manipulation

### `simple_api_example.py`
Minimal example showing the basic V4C API usage:
- Extract Virtual 4C data
- Generate plots

## Usage

### Command Line Examples
```bash
# Extract data
v4c-extract sample.mcool --resolutions 10000 50000 --genes MYC,GATA6 --genome hg38 --output data.tsv

# Plot results
v4c-plot --input data.tsv --ylim 0.4 --output plots.png

# Compare samples
v4c-compare --inputs sample1.tsv sample2.tsv --ylim 1.0 --scale --output comparison.png
```

### Python API Examples
```python
from v4c.extract import extract_v4c
from v4c.plot import plot_v4c
from v4c.compare import compare_v4c

# Extract data
extract_v4c(
    mcool_files=["sample.mcool"],
    resolutions=[10000, 50000],
    genes="MYC,GATA6",
    genome="hg38",
    flank=500000,
    normalization_method="minmax",
    output="data.tsv"
)

# Plot results
plot_v4c(
    input_file="data.tsv",
    ylim=0.4,
    output_file="plots.png"
)

# Compare samples
compare_v4c(
    input_files=["sample1.tsv", "sample2.tsv"],
    ylim=1.0,
    scale=True,
    output_file="comparison.png"
)
```

## Requirements

- V4C package installed (`pip install -e .`)
- .mcool files with Hi-C data
- Python 3.7 or higher

## Notes

- Replace file paths in examples with your actual .mcool files
- The API examples assume you have valid Hi-C data files
- Check the generated output files for results
