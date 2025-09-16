# V4C - Virtual 4C Analysis Tool

A Python package for extracting and analyzing Virtual 4C contact frequencies from Hi-C data stored in .mcool files.

## Overview

Virtual 4C (V4C) analysis extracts contact frequencies around specific genomic regions (viewpoints) from Hi-C data, providing insights into chromatin interactions and 3D genome organization. This tool supports multiple input formats, normalization methods, and visualization options.

## Features

- Extract Virtual 4C profiles from .mcool files
- Support for gene names, genomic coordinates, and BED files
- Multiple normalization methods (Min-Max scaling, self-normalization)
- Flexible flanking region specification
- High-quality visualization with customizable colors and labels
- Comparison analysis between multiple samples
- Command-line interface for easy integration into workflows

## Installation

### Prerequisites

- Python 3.7 or higher
- Required packages: cooler, numpy, pandas, matplotlib
- Input data: .mcool files (multi-resolution cooler format)

### Input Data Format

V4C requires .mcool files as input. If you have other Hi-C data formats, you can convert them:

**Converting from .hic files to .mcool:**
```bash
pip install hic2cool
hic2cool convert <infile.hic> <outfile.mcool> -r <resolution> -p <nproc>
```

For more details on hic2cool conversion options, visit the [hic2cool repository](https://github.com/4dn-dcic/hic2cool).

**Converting from pairs format:**
You can also convert pairs format to .mcool using cooler tools or the hic2cool package mentioned above.

### Install from source

```bash
git clone <repository-url>
cd V4C
pip install -e .
```

**Note**: Installing in development mode (`pip install -e .`) is recommended to avoid import warnings and enable direct command-line usage.

## Quick Start

### Basic Usage

1. **Extract Virtual 4C data**:
```bash
v4c-extract sample.mcool --resolutions 5000 10000 --genes MYC,DDX11L1 --genome hg38 --output data.tsv
```

2. **Plot individual samples**:
```bash
v4c-plot --input data.tsv --ylim 0.4 --output plots.png
```

3. **Compare multiple samples**:
```bash
v4c-compare --inputs sample1.tsv sample2.tsv --ylim 1.0 --scale --output comparison.png
```

## Command Reference

### v4c-extract

Extract Virtual 4C contact frequencies from .mcool files.

**Required arguments**:
- `mcool_files`: One or more .mcool files
- `--resolutions`: List of resolutions to extract (e.g., 5000 10000 50000)

**Input options** (choose one):
- `--genes`: Comma-separated gene names (requires --genome)
- `--coords`: Genomic coordinates (e.g., 'chr17:45878152-46000000')
- `--bed`: Path to BED file with genomic regions (3 columns: chrom, start, end; or 4+ columns: chrom, start, end, gene_name)

**Optional parameters**:
- `--genome`: Reference genome version (hg38 or hg19)
- `--flank`: Flanking region in base pairs (default: 500000)
- `--normalization`: Normalization method - 'minmax' or 'self' (default: minmax)
- `--use-fixed-center`: Use fixed center position calculation
- `--no-balance`: Disable ICE balancing
- `--no-scale`: Disable normalization
- `--output`: Output file name (default: extracted_data.tsv)

**Examples**:
```bash
# Extract by gene names
v4c-extract sample.mcool --resolutions 5000 10000 --genes MYC,DDX11L1 --genome hg38 --flank 500000 --normalization minmax --output data.tsv

# Extract by coordinates
v4c-extract sample.mcool --resolutions 50000 --coords chr8:127732934-127737934 --genome hg38 --output data.tsv

# Extract from BED file (3 columns: chrom, start, end)
v4c-extract sample.mcool --resolutions 50000 --bed regions.bed --output data.tsv

# Extract from BED file (4 columns: chrom, start, end, gene_name)
v4c-extract sample.mcool --resolutions 50000 --bed genes.bed --output data.tsv
```

### v4c-plot

Plot Virtual 4C contact frequencies.

**Required arguments**:
- `--input`: Input TSV file from v4c-extract

**Optional parameters**:
- `--ylim`: Maximum y-axis value (default: 0.4)
- `--flank`: Flanking region in bp for display
- `--output`: Output file path for saving the plot
- `--dpi`: DPI for the output figure (default: 300)
- `--sample-names`: Custom sample names as JSON dict
- `--colors`: Custom colors as JSON list

**Examples**:
```bash
# Basic plotting
v4c-plot --input data.tsv --ylim 0.4 --output plots.png

# With custom sample names and colors
v4c-plot --input data.tsv --ylim 0.4 --sample-names '{"sample1.mcool": "Control", "sample2.mcool": "Treatment"}' --colors '["#ff0000", "#00ff00"]'

# Basic plotting with gene names from extract
v4c-plot --input data.tsv --ylim 0.4 --output plots.png
```

### v4c-compare

Compare Virtual 4C data from multiple files.

**Required arguments**:
- `--inputs`: Input TSV files from v4c-extract

**Optional parameters**:
- `--ylim`: Maximum y-axis value (default: 1.0)
- `--scale`: Normalize values between 0 and 1
- `--output`: Output file path for saving the plot
- `--dpi`: DPI for the output figure (default: 300)
- `--sample-names`: Custom sample names as JSON dict
- `--colors`: Custom colors as JSON list

**Examples**:
```bash
# Basic comparison
v4c-compare --inputs sample1.tsv sample2.tsv --ylim 1.0 --scale --output comparison.png

# With custom sample names and colors
v4c-compare --inputs sample1.tsv sample2.tsv --ylim 1.0 --scale --sample-names '{"sample1.mcool": "Control", "sample2.mcool": "Treatment"}' --colors '["#ff0000", "#00ff00"]'
```

## Tutorial

### Step 1: Prepare Your Data

Ensure you have:
- .mcool files containing Hi-C data
- Gene names or genomic coordinates of interest
- Optional: BED file with custom genomic regions

**Data Format Conversion:**
If your Hi-C data is in .hic format, convert it to .mcool first:
```bash
pip install hic2cool
hic2cool convert input.hic output.mcool -r 5000,10000,25000,50000 -p 4
```

For pairs format or other conversions, refer to the [hic2cool documentation](https://github.com/4dn-dcic/hic2cool).

### Step 2: Extract Virtual 4C Data

Choose your input method:

**Method 1: Gene names** (recommended for most users)
```bash
v4c-extract sample.mcool --resolutions 5000 10000 --genes MYC,GATA6,KLF5 --genome hg38 --flank 500000 --normalization minmax --output extracted_data.tsv
```

**Method 2: Genomic coordinates**
```bash
v4c-extract sample.mcool --resolutions 50000 --coords chr8:127732934-127737934 --genome hg38 --flank 500000 --output extracted_data.tsv
```

**Method 3: BED file**
```bash
# 3-column BED file (chrom, start, end)
v4c-extract sample.mcool --resolutions 50000 --bed regions.bed --flank 500000 --output extracted_data.tsv

# 4-column BED file (chrom, start, end, gene_name)
v4c-extract sample.mcool --resolutions 50000 --bed genes.bed --flank 500000 --output extracted_data.tsv
```

### Step 3: Visualize Results

**Individual plots** (one plot per sample per gene):
```bash
v4c-plot --input extracted_data.tsv --ylim 0.4 --flank 500000 --output individual_plots.png
```

**Comparison plots** (one plot per gene comparing samples):
```bash
v4c-compare --inputs sample1.tsv sample2.tsv --ylim 1.0 --scale --output comparison_plots.png
```

### Step 4: Customize Visualization

**Custom colors and sample names**:
```bash
v4c-plot --input data.tsv --ylim 0.4 --colors '["#ff0000", "#00ff00", "#0000ff"]' --sample-names '{"sample1.mcool": "Control", "sample2.mcool": "Treatment"}' --output custom_plots.png
```

## Output Files

### TSV Format

The extracted data is saved in TSV format with the following columns:
- `mcool`: Input file path
- `res`: Resolution used
- `chrom`: Chromosome
- `start`: Start position
- `end`: End position
- `gene_name`: Gene name (if available)
- Contact frequency columns: Genomic coordinates as column names

### Plot Features

- **Individual plots**: One plot per sample per gene/region
- **Comparison plots**: One plot per gene/region comparing all samples
- **Genomic coordinates**: X-axis shows actual genomic positions
- **Scientific notation**: Large coordinates displayed in scientific format
- **Coordinate range**: Annotated region boundaries
- **Custom colors**: User-defined color schemes
- **Sample labels**: Customizable sample names in legends

## Normalization Methods

### Min-Max Normalization (default)
Scales contact frequencies to range [0, 1] using the formula:
```
normalized_value = (value - min_value) / (max_value - min_value)
```

### Self-Normalization
Divides all contact frequencies by the contact frequency at the viewpoint (TSS):
```
normalized_value = value / viewpoint_value
```

## Troubleshooting

### Common Issues

1. **Import warnings**: Install in development mode with `pip install -e .`
2. **File not found**: Check file paths and ensure .mcool files exist
3. **Wrong file format**: V4C requires .mcool files. Convert .hic files using [hic2cool](https://github.com/4dn-dcic/hic2cool)
4. **Out of bounds errors**: Some genomic regions may be outside Hi-C data coverage
5. **Memory issues**: Use higher resolutions (e.g., 50000 instead of 5000) for large datasets

### Getting Help

- Check command help: `v4c-extract --help`, `v4c-plot --help`, `v4c-compare --help`
- Verify input file formats and parameters
- Ensure all required dependencies are installed

## Dependencies

- cooler: Hi-C data processing
- numpy: Numerical computations
- pandas: Data manipulation
- matplotlib: Plotting and visualization

## License

[Add your license information here]

## Citation

If you use this tool in your research, please cite:

[Add citation information here]