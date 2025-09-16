import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from typing import List, Optional, Union
from .extract import V4CError, InputValidationError, FileProcessingError

def validate_compare_inputs(input_files: List[str], ylim: float, scale: bool) -> None:
    """
    Validates input parameters for compare_v4c function.
    
    Args:
        input_files: List of input TSV files
        ylim: Maximum y-axis value
        scale: Whether to normalize values
        
    Raises:
        InputValidationError: If any input parameter is invalid
    """
    if not input_files:
        raise InputValidationError("No input files provided")
    
    for file in input_files:
        if not os.path.exists(file):
            raise InputValidationError(f"Input file not found: {file}")
        if not file.endswith('.tsv'):
            raise InputValidationError(f"Input file must be a TSV file: {file}")
    
    if ylim <= 0:
        raise InputValidationError(f"ylim must be positive: {ylim}")

def compare_v4c(input_files: List[str], 
                 ylim: float = 1.0, 
                 scale: bool = True,
                 output_file: Optional[str] = None,
                 dpi: int = 300,
                 figsize: tuple = (12, 8),
                 colors: Optional[List[str]] = None,
                 sample_names: Optional[dict] = None) -> None:
    """
    Compares Virtual 4C data from multiple .mcool files.

    Args:
        input_files: List of TSV files containing extracted V4C data
        ylim: Maximum y-axis value
        scale: Whether to normalize values between 0 and 1
        output_file: Optional path to save the plot
        dpi: DPI for the output figure
        figsize: Figure size (width, height) in inches
        colors: Optional list of colors for different datasets
        sample_names: Optional dict mapping mcool filenames to custom sample names

    Returns:
        None: Displays or saves the comparison plot

    Raises:
        InputValidationError: If input parameters are invalid
        FileProcessingError: If file processing fails
        V4CError: For other V4C-related errors
    """
    try:
        # Validate inputs
        validate_compare_inputs(input_files, ylim, scale)
        
        
        # Read and validate data
        dfs = []
        has_gene_name = False
        for file in input_files:
            try:
                df = pd.read_csv(file, sep="\t")
                required_columns = ["mcool", "res", "chrom", "start", "end"]
                if not all(col in df.columns for col in required_columns):
                    raise InputValidationError(f"Input file missing required columns: {required_columns}")
                # Check if any file has gene_name column
                if "gene_name" in df.columns:
                    has_gene_name = True
                dfs.append(df)
            except Exception as e:
                raise FileProcessingError(f"Error reading input file {file}: {str(e)}")

        # Create comparison plot for each coordinate
        # Group all data by coordinate first
        all_coords = set()
        for df in dfs:
            for _, row in df.iterrows():
                all_coords.add((row['chrom'], row['start'], row['end']))
        
        # Use default colors if not provided
        if colors is None:
            colors = plt.cm.tab10.colors[:len(input_files)]
        
        # Create a plot for each coordinate
        print(f"Found {len(all_coords)} coordinates to process: {all_coords}")
        for i, coord in enumerate(all_coords):
            chrom, start, end = coord
            print(f"Processing coordinate {i+1}/{len(all_coords)}: {chrom}:{start}-{end}")
            plt.figure(figsize=figsize)
            
            # Plot each dataset for this coordinate
            for df, color in zip(dfs, colors):
                # Filter data for this coordinate
                coord_data = df[(df['chrom'] == chrom) & (df['start'] == start) & (df['end'] == end)]
                
                for _, row in coord_data.iterrows():
                    try:
                        # Skip metadata columns (mcool, res, chrom, start, end, gene_name)
                        coords = list(map(float, row[6:].values))
                        if scale:
                            min_val, max_val = min(coords), max(coords)
                            coords = [(x - min_val) / (max_val - min_val) if max_val > min_val else 0 for x in coords]
                        
                        # Create genomic coordinates for x-axis
                        resolution = row['res']
                        # Calculate flank for this specific row
                        num_contacts = len(coords)
                        row_flank = (num_contacts * resolution) // 2
                        # Ensure genomic_coords has the same length as coords
                        genomic_coords = np.linspace(start - row_flank, end + row_flank, len(coords))
                        
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
                        
                        plt.plot(genomic_coords, coords, 
                                label=sample_label,
                                color=color,
                                alpha=0.7,
                                linewidth=2)
                    except Exception as e:
                        print(f"Warning: Error plotting data for {row['mcool']} at {coord}: {str(e)}")
                        continue

            # Customize plot
            plt.xlabel("Genomic Position (bp)", fontsize=12)
            plt.ylabel("Hi-C Contact Frequency", fontsize=12)
            
            # Format x-axis with scientific notation for large numbers
            ax = plt.gca()
            ax.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))
            
            # Add coordinate range annotation
            # Calculate flank from the data (assuming it's the same for all samples)
            if len(dfs) > 0:
                sample_data = dfs[0]
                coord_data = sample_data[(sample_data['chrom'] == chrom) & (sample_data['start'] == start) & (sample_data['end'] == end)]
                if not coord_data.empty:
                    # Calculate flank from the number of contact points
                    num_contacts = len(coord_data.columns) - 6  # 6 metadata columns
                    # Assuming resolution is consistent, estimate flank
                    resolution = coord_data['res'].iloc[0]
                    estimated_flank = (num_contacts * resolution) // 2
                    coord_range = f"{chrom}:{start-estimated_flank:,}-{end+estimated_flank:,}"
                else:
                    coord_range = f"{chrom}:{start:,}-{end:,}"
            else:
                coord_range = f"{chrom}:{start:,}-{end:,}"
            
            plt.text(0.02, 0.98, f"Region: {coord_range}", 
                    transform=plt.gca().transAxes, 
                    fontsize=10, 
                    verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            # Create title with gene name if available, otherwise use coordinate
            coord_key = (chrom, start, end)
            
            # Priority 1: Use gene_name from data if available
            gene_name_from_data = None
            if has_gene_name:
                for df in dfs:
                    coord_data = df[(df['chrom'] == chrom) & (df['start'] == start) & (df['end'] == end)]
                    if not coord_data.empty and coord_data['gene_name'].iloc[0] and str(coord_data['gene_name'].iloc[0]).strip():
                        gene_name_from_data = coord_data['gene_name'].iloc[0]
                        break
            
            if gene_name_from_data:
                plt.title(f"Virtual 4C Comparison - {gene_name_from_data}", fontsize=14)
            else:
                coord_label = f"{chrom}:{start}-{end}"
                plt.title(f"Virtual 4C Comparison - {coord_label}", fontsize=14)
            
            plt.legend()
            plt.ylim(0, ylim)
            
            # Save or show plot
            if output_file:
                try:
                    # Create unique filename for each coordinate
                    base_name = output_file.replace('.png', '').replace('.pdf', '')
                    
                    # Use gene name in filename if available, otherwise use coordinate
                    if gene_name_from_data:
                        coord_suffix = f"_{gene_name_from_data}"
                    else:
                        coord_suffix = f"_{chrom}_{start}_{end}"
                    unique_output = f"{base_name}{coord_suffix}.png"
                    print(f"  Saving plot to: {unique_output}")
                    plt.savefig(unique_output, dpi=dpi, bbox_inches='tight')
                    plt.close()
                    print(f"  Plot saved successfully")
                except Exception as e:
                    raise FileProcessingError(f"Error saving plot to {output_file}: {str(e)}")
            else:
                plt.show()
                plt.close()

    except V4CError:
        raise
    except Exception as e:
        raise V4CError(f"Unexpected error in compare_v4c: {str(e)}")

def main():
    """Command-line interface for compare_v4c"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Compare Virtual 4C data from multiple files")
    parser.add_argument("--inputs", required=True, nargs="+", help="Input TSV files from v4c-extract")
    parser.add_argument("--ylim", type=float, default=1.0, help="Maximum y-axis value")
    parser.add_argument("--scale", action="store_true", help="Normalize values between 0 and 1")
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
        
        compare_v4c(args.inputs, args.ylim, args.scale, args.output, args.dpi, sample_names=sample_names, colors=colors)
    except V4CError as e:
        print(f"Error: {str(e)}")
        exit(1)
    except Exception as e:
        print(f"Unexpected error: {str(e)}")
        exit(1)

if __name__ == "__main__":
    main()
