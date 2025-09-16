import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from typing import List, Optional, Union
from .extract import V4CError, InputValidationError, FileProcessingError

def validate_plot_inputs(input_file: str, ylim: float, flank: int) -> None:
    """
    Validates input parameters for plot_v4c function.
    
    Args:
        input_file: Path to input TSV file
        ylim: Maximum y-axis value
        flank: Flanking region in bp
        
    Raises:
        InputValidationError: If any input parameter is invalid
    """
    if not os.path.exists(input_file):
        raise InputValidationError(f"Input file not found: {input_file}")
    
    if not input_file.endswith('.tsv'):
        raise InputValidationError(f"Input file must be a TSV file: {input_file}")
    
    if ylim <= 0:
        raise InputValidationError(f"ylim must be positive: {ylim}")
    
    if flank <= 0:
        raise InputValidationError(f"flank must be positive: {flank}")

def plot_v4c(input_file: str, 
             ylim: float = 0.4, 
             flank: int = 50000,
             output_file: Optional[str] = None,
             dpi: int = 300,
             figsize: tuple = (10, 6),
             sample_names: Optional[dict] = None,
             colors: Optional[List[str]] = None) -> None:
    """
    Plots Virtual 4C contact frequencies.

    Args:
        input_file: TSV file containing extracted V4C data
        ylim: Maximum y-axis value
        flank: Flanking region (bp) to display
        output_file: Optional path to save the plot
        dpi: DPI for the output figure
        figsize: Figure size (width, height) in inches
        sample_names: Optional dict mapping mcool filenames to custom sample names
        colors: Optional list of colors for different samples

    Returns:
        None: Displays or saves the plot

    Raises:
        InputValidationError: If input parameters are invalid
        FileProcessingError: If file processing fails
        V4CError: For other V4C-related errors
    """
    try:
        # Validate inputs
        validate_plot_inputs(input_file, ylim, flank)
        
        
        # Read and validate data
        try:
            df = pd.read_csv(input_file, sep="\t")
            required_columns = ["mcool", "res", "chrom", "start", "end"]
            if not all(col in df.columns for col in required_columns):
                raise InputValidationError(f"Input file missing required columns: {required_columns}")
            
            # Check if gene_name column exists
            has_gene_name = "gene_name" in df.columns
        except Exception as e:
            raise FileProcessingError(f"Error reading input file: {str(e)}")

        # Set up colors
        if colors is None:
            colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        
        # Create plots for each sample and gene combination
        # Group by mcool file and coordinate (gene)
        for (mcool_file, chrom, start, end), coord_group in df.groupby(['mcool', 'chrom', 'start', 'end']):
            plt.figure(figsize=figsize)
            
            # Plot each sample for this coordinate
            for i, (_, row) in enumerate(coord_group.iterrows()):
                try:
                    # Skip metadata columns (mcool, res, chrom, start, end, gene_name)
                    coords = list(map(float, row[6:].values))
                    # Create sample label from mcool filename
                    mcool_filename = os.path.basename(row['mcool'])
                    
                    # Check if user provided custom sample names
                    if sample_names and mcool_filename in sample_names:
                        sample_label = sample_names[mcool_filename]
                    else:
                        # Use default naming: extract prefix from mcool filename
                        sample_name = mcool_filename.replace('.mcool', '')
                        # Extract cell type and genome info if available
                        if 'genome1' in sample_name.lower():
                            # Extract cell type (e.g., "Astro" from "Astro_merged_chr17.genome1")
                            cell_type = sample_name.split('_')[0]
                            sample_label = f"{cell_type} Genome1"
                        elif 'genome2' in sample_name.lower():
                            # Extract cell type (e.g., "Astro" from "Astro_merged_chr17.genome2")
                            cell_type = sample_name.split('_')[0]
                            sample_label = f"{cell_type} Genome2"
                        else:
                            sample_label = sample_name
                    
                    # Create genomic coordinates for x-axis
                    resolution = row['res']
                    # Ensure genomic_coords has the same length as coords
                    genomic_coords = np.linspace(start - flank, end + flank, len(coords))
                    
                    # Use provided color or default (always use first color since each group has one sample)
                    color = colors[0] if colors else None
                    
                    plt.plot(genomic_coords, coords, 
                            label=sample_label,
                            color=color,
                            alpha=0.7,
                            linewidth=2)
                except Exception as e:
                    raise FileProcessingError(f"Error plotting data: {str(e)}")

            # Customize plot
            plt.xlabel("Genomic Position (bp)", fontsize=12)
            plt.ylabel("Hi-C Contact Frequency", fontsize=12)
            
            # Format x-axis with scientific notation for large numbers
            ax = plt.gca()
            ax.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))
            
            # Add coordinate range annotation
            coord_range = f"{chrom}:{start-flank:,}-{end+flank:,}"
            plt.text(0.02, 0.98, f"Region: {coord_range}", 
                    transform=plt.gca().transAxes, 
                    fontsize=10, 
                    verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            # Create title with gene name if available, otherwise use coordinate
            coord_key = (chrom, start, end)
            
            # Use gene_name from data if available, otherwise use coordinate
            if has_gene_name and coord_group['gene_name'].iloc[0] and str(coord_group['gene_name'].iloc[0]).strip():
                gene_name = coord_group['gene_name'].iloc[0]
                title = f"Virtual 4C - {gene_name}"
            else:
                coord_label = f"{chrom}:{start}-{end}"
                title = f"Virtual 4C - {coord_label}"
            plt.title(title, fontsize=14)
            plt.legend()
            plt.ylim(0, ylim)
            
            # Save or show plot
            if output_file:
                try:
                    # Create unique filename for each coordinate
                    base_name = output_file.replace('.png', '').replace('.pdf', '')
                    
                    # Use gene name and sample info in filename
                    if has_gene_name and coord_group['gene_name'].iloc[0] and str(coord_group['gene_name'].iloc[0]).strip() and str(coord_group['gene_name'].iloc[0]) != 'nan':
                        gene_name = coord_group['gene_name'].iloc[0]
                        # Extract sample name from mcool file path
                        sample_name = os.path.basename(mcool_file).replace('.mcool', '')
                        coord_suffix = f"_{gene_name}_{sample_name}"
                    else:
                        sample_name = os.path.basename(mcool_file).replace('.mcool', '')
                        coord_suffix = f"_{chrom}_{start}_{end}_{sample_name}"
                    unique_output = f"{base_name}{coord_suffix}.png"
                    plt.savefig(unique_output, dpi=dpi, bbox_inches='tight')
                    plt.close()
                except Exception as e:
                    raise FileProcessingError(f"Error saving plot to {output_file}: {str(e)}")
            else:
                plt.show()
                plt.close()

    except V4CError:
        raise
    except Exception as e:
        raise V4CError(f"Unexpected error in plot_v4c: {str(e)}")

def main():
    """Command-line interface for plot_v4c"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Plot Virtual 4C contact frequencies")
    parser.add_argument("--input", required=True, help="Input TSV file from v4c-extract")
    parser.add_argument("--ylim", type=float, default=0.4, help="Maximum y-axis value")
    parser.add_argument("--flank", type=int, default=50000, help="Flanking region in bp")
    parser.add_argument("--output", help="Output file path for saving the plot")
    parser.add_argument("--dpi", type=int, default=300, help="DPI for the output figure")
    parser.add_argument("--sample-names", help="Custom sample names as JSON dict (e.g., '{\"file1.mcool\": \"Sample1\", \"file2.mcool\": \"Sample2\"}')")
    parser.add_argument("--colors", help="Custom colors as JSON list (e.g., '[\"#ff0000\", \"#00ff00\", \"#0000ff\"]')")
    
    args = parser.parse_args()
    
    try:
        # Parse sample names if provided
        sample_names = None
        if args.sample_names:
            import json
            sample_names = json.loads(args.sample_names)
        
        # Parse colors if provided
        colors = None
        if args.colors:
            import json
            colors = json.loads(args.colors)
        
        plot_v4c(args.input, args.ylim, args.flank, args.output, args.dpi, sample_names=sample_names, colors=colors)
    except V4CError as e:
        print(f"Error: {str(e)}")
        exit(1)
    except Exception as e:
        print(f"Unexpected error: {str(e)}")
        exit(1)

if __name__ == "__main__":
    main()
