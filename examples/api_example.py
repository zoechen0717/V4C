#!/usr/bin/env python3
"""
V4C API Example
===============

This example demonstrates how to use V4C as a Python package instead of command-line tools.
This is useful for integrating V4C into larger analysis pipelines or custom scripts.

Requirements:
- V4C package installed
- .mcool files with Hi-C data
- Python 3.7+
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

# Add the parent directory to the path to import v4c
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from v4c.extract import extract_v4c
from v4c.plot import plot_v4c
from v4c.compare import compare_v4c


def example_basic_extraction():
    """
    Example 1: Basic extraction using gene names
    """
    print("=== Example 1: Basic Extraction ===")
    
    # Define input parameters
    mcool_files = ["/path/to/your/sample.mcool"]  # Replace with actual file path
    resolutions = [5000, 10000, 50000]
    genes = "MYC,GATA6,KLF5"
    genome = "hg38"
    flank = 500000
    normalization_method = "minmax"
    output_file = "api_extracted_data.tsv"
    
    try:
        # Extract Virtual 4C data
        extract_v4c(
            mcool_files=mcool_files,
            resolutions=resolutions,
            genes=genes,
            genome=genome,
            flank=flank,
            normalization_method=normalization_method,
            output=output_file
        )
        
        print(f"✓ Data extracted successfully to {output_file}")
        
        # Read and display the results
        df = pd.read_csv(output_file, sep='\t')
        print(f"✓ Extracted data shape: {df.shape}")
        print(f"✓ Columns: {list(df.columns[:10])}...")  # Show first 10 columns
        
    except Exception as e:
        print(f"✗ Error in extraction: {e}")


def example_coordinate_extraction():
    """
    Example 2: Extraction using genomic coordinates
    """
    print("\n=== Example 2: Coordinate-based Extraction ===")
    
    # Define input parameters
    mcool_files = ["/path/to/your/sample.mcool"]  # Replace with actual file path
    resolutions = [10000, 50000]
    coords = "chr8:127732934-127737934"  # MYC gene region
    genome = "hg38"
    flank = 200000
    normalization_method = "self"
    output_file = "api_coordinate_data.tsv"
    
    try:
        # Extract Virtual 4C data
        extract_v4c(
            mcool_files=mcool_files,
            resolutions=resolutions,
            coords=coords,
            genome=genome,
            flank=flank,
            normalization_method=normalization_method,
            output=output_file
        )
        
        print(f"✓ Coordinate data extracted successfully to {output_file}")
        
    except Exception as e:
        print(f"✗ Error in coordinate extraction: {e}")


def example_bed_file_extraction():
    """
    Example 3: Extraction using BED file
    """
    print("\n=== Example 3: BED File Extraction ===")
    
    # Create a sample BED file
    bed_content = """chr8	127732934	127737934	MYC
chr18	22166892	22171892	GATA6
chr13	73052476	73057476	KLF5"""
    
    bed_file = "sample_regions.bed"
    with open(bed_file, 'w') as f:
        f.write(bed_content)
    
    # Define input parameters
    mcool_files = ["/path/to/your/sample.mcool"]  # Replace with actual file path
    resolutions = [50000]
    bed_file_path = bed_file
    flank = 100000
    normalization_method = "minmax"
    output_file = "api_bed_data.tsv"
    
    try:
        # Extract Virtual 4C data
        extract_v4c(
            mcool_files=mcool_files,
            resolutions=resolutions,
            bed_file=bed_file_path,
            flank=flank,
            normalization_method=normalization_method,
            output=output_file
        )
        
        print(f"✓ BED file data extracted successfully to {output_file}")
        
        # Clean up
        os.remove(bed_file)
        
    except Exception as e:
        print(f"✗ Error in BED file extraction: {e}")
        if os.path.exists(bed_file):
            os.remove(bed_file)


def example_plotting():
    """
    Example 4: Plotting extracted data
    """
    print("\n=== Example 4: Plotting ===")
    
    # Assume we have extracted data
    input_file = "api_extracted_data.tsv"
    
    if not os.path.exists(input_file):
        print(f"✗ Input file {input_file} not found. Run extraction examples first.")
        return
    
    try:
        # Plot with default settings
        plot_v4c(
            input_file=input_file,
            ylim=0.4,
            flank=500000,
            output_file="api_plots.png",
            dpi=300
        )
        
        print("✓ Plots generated successfully")
        
        # Plot with custom colors and sample names
        custom_colors = ["#ff0000", "#00ff00", "#0000ff"]
        sample_names = {
            "sample.mcool": "Control Sample"
        }
        
        plot_v4c(
            input_file=input_file,
            ylim=0.5,
            flank=500000,
            output_file="api_custom_plots.png",
            dpi=300,
            colors=custom_colors,
            sample_names=sample_names
        )
        
        print("✓ Custom plots generated successfully")
        
    except Exception as e:
        print(f"✗ Error in plotting: {e}")


def example_comparison():
    """
    Example 5: Comparing multiple samples
    """
    print("\n=== Example 5: Sample Comparison ===")
    
    # Assume we have multiple extracted data files
    input_files = ["sample1_data.tsv", "sample2_data.tsv"]
    
    # Check if files exist
    missing_files = [f for f in input_files if not os.path.exists(f)]
    if missing_files:
        print(f"✗ Missing files: {missing_files}")
        print("Create multiple sample data files first.")
        return
    
    try:
        # Compare with default settings
        compare_v4c(
            input_files=input_files,
            ylim=1.0,
            scale=True,
            output_file="api_comparison.png",
            dpi=300
        )
        
        print("✓ Comparison plots generated successfully")
        
        # Compare with custom settings
        custom_colors = ["#ff0000", "#00ff00"]
        sample_names = {
            "sample1.mcool": "Treatment",
            "sample2.mcool": "Control"
        }
        
        compare_v4c(
            input_files=input_files,
            ylim=1.0,
            scale=True,
            output_file="api_custom_comparison.png",
            dpi=300,
            colors=custom_colors,
            sample_names=sample_names
        )
        
        print("✓ Custom comparison plots generated successfully")
        
    except Exception as e:
        print(f"✗ Error in comparison: {e}")


def example_workflow():
    """
    Example 6: Complete workflow
    """
    print("\n=== Example 6: Complete Workflow ===")
    
    # This example shows a complete analysis workflow
    # Note: Replace file paths with actual .mcool files
    
    mcool_files = [
        "/path/to/control.mcool",
        "/path/to/treatment.mcool"
    ]
    
    # Step 1: Extract data for both samples
    for i, mcool_file in enumerate(mcool_files):
        if not os.path.exists(mcool_file):
            print(f"✗ File not found: {mcool_file}")
            continue
            
        output_file = f"sample_{i+1}_data.tsv"
        
        try:
            extract_v4c(
                mcool_files=[mcool_file],
                resolutions=[10000, 50000],
                genes="MYC,GATA6",
                genome="hg38",
                flank=500000,
                normalization_method="minmax",
                output=output_file
            )
            print(f"✓ Extracted data for sample {i+1}")
            
        except Exception as e:
            print(f"✗ Error extracting sample {i+1}: {e}")
    
    # Step 2: Generate individual plots
    for i in range(len(mcool_files)):
        input_file = f"sample_{i+1}_data.tsv"
        if os.path.exists(input_file):
            try:
                plot_v4c(
                    input_file=input_file,
                    ylim=0.4,
                    output_file=f"sample_{i+1}_plots.png"
                )
                print(f"✓ Generated plots for sample {i+1}")
            except Exception as e:
                print(f"✗ Error plotting sample {i+1}: {e}")
    
    # Step 3: Compare samples
    input_files = [f"sample_{i+1}_data.tsv" for i in range(len(mcool_files))]
    existing_files = [f for f in input_files if os.path.exists(f)]
    
    if len(existing_files) >= 2:
        try:
            compare_v4c(
                input_files=existing_files,
                ylim=1.0,
                scale=True,
                output_file="workflow_comparison.png",
                colors=["#ff0000", "#00ff00"],
                sample_names={
                    "control.mcool": "Control",
                    "treatment.mcool": "Treatment"
                }
            )
            print("✓ Generated comparison plots")
        except Exception as e:
            print(f"✗ Error in comparison: {e}")


def example_data_analysis():
    """
    Example 7: Data analysis and manipulation
    """
    print("\n=== Example 7: Data Analysis ===")
    
    # Load extracted data
    input_file = "api_extracted_data.tsv"
    
    if not os.path.exists(input_file):
        print(f"✗ Input file {input_file} not found.")
        return
    
    try:
        # Read the data
        df = pd.read_csv(input_file, sep='\t')
        print(f"✓ Loaded data: {df.shape}")
        
        # Basic statistics
        print(f"✓ Number of genes: {df['gene_name'].nunique()}")
        print(f"✓ Number of resolutions: {df['res'].nunique()}")
        print(f"✓ Resolutions: {sorted(df['res'].unique())}")
        
        # Get contact frequency columns
        contact_cols = [col for col in df.columns if col not in 
                       ['mcool', 'res', 'chrom', 'start', 'end', 'gene_name']]
        
        print(f"✓ Number of contact bins: {len(contact_cols)}")
        
        # Calculate statistics for each gene
        for gene in df['gene_name'].unique():
            if pd.isna(gene) or gene == '':
                continue
                
            gene_data = df[df['gene_name'] == gene]
            print(f"\n--- {gene} ---")
            
            for _, row in gene_data.iterrows():
                contact_values = row[contact_cols].values
                print(f"  Resolution {row['res']}: mean={contact_values.mean():.4f}, "
                      f"max={contact_values.max():.4f}, min={contact_values.min():.4f}")
        
        # Find peak contact frequencies
        for gene in df['gene_name'].unique():
            if pd.isna(gene) or gene == '':
                continue
                
            gene_data = df[df['gene_name'] == gene]
            for _, row in gene_data.iterrows():
                contact_values = row[contact_cols].values
                peak_idx = contact_values.argmax()
                peak_value = contact_values.max()
                peak_position = contact_cols[peak_idx]
                
                print(f"  {gene} (res {row['res']}): Peak at position {peak_position} "
                      f"with value {peak_value:.4f}")
        
    except Exception as e:
        print(f"✗ Error in data analysis: {e}")


def main():
    """
    Main function to run all examples
    """
    print("V4C API Examples")
    print("================")
    print("This script demonstrates how to use V4C as a Python package.")
    print("Note: Replace file paths with actual .mcool files before running.")
    print()
    
    # Run examples
    example_basic_extraction()
    example_coordinate_extraction()
    example_bed_file_extraction()
    example_plotting()
    example_comparison()
    example_workflow()
    example_data_analysis()
    
    print("\n=== Summary ===")
    print("✓ All examples completed!")
    print("✓ Check the generated files:")
    print("  - api_extracted_data.tsv")
    print("  - api_coordinate_data.tsv")
    print("  - api_bed_data.tsv")
    print("  - api_plots.png")
    print("  - api_custom_plots.png")
    print("  - api_comparison.png")
    print("  - api_custom_comparison.png")
    print("  - workflow_comparison.png")
    
    print("\nFor more information, see the V4C documentation.")


if __name__ == "__main__":
    main()
