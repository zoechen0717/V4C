import matplotlib.pyplot as plt
import pandas as pd
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
                colors: Optional[List[str]] = None) -> None:
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
        for file in input_files:
            try:
                df = pd.read_csv(file, sep="\t")
                required_columns = ["mcool", "res", "chrom", "start", "end"]
                if not all(col in df.columns for col in required_columns):
                    raise InputValidationError(f"Input file missing required columns: {required_columns}")
                dfs.append(df)
            except Exception as e:
                raise FileProcessingError(f"Error reading input file {file}: {str(e)}")

        # Create comparison plot
        plt.figure(figsize=figsize)
        
        # Use default colors if not provided
        if colors is None:
            colors = plt.cm.tab10.colors[:len(input_files)]
        
        # Plot each dataset
        for df, color in zip(dfs, colors):
            for _, row in df.iterrows():
                try:
                    coords = list(map(float, row[5:].values))
                    if scale:
                        min_val, max_val = min(coords), max(coords)
                        coords = [(x - min_val) / (max_val - min_val) if max_val > min_val else 0 for x in coords]
                    
                    plt.plot(coords, 
                            label=f"{os.path.basename(row['mcool'])} - Res {row['res']}",
                            color=color,
                            alpha=0.7)
                except Exception as e:
                    raise FileProcessingError(f"Error plotting data for {row['mcool']}: {str(e)}")

        # Customize plot
        plt.xlabel("Genomic Position (bp)", fontsize=12)
        plt.ylabel("Hi-C Contact Frequency", fontsize=12)
        plt.title("Virtual 4C Comparison", fontsize=14)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.ylim(0, ylim)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        # Save or show plot
        if output_file:
            try:
                plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
                plt.close()
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
    
    args = parser.parse_args()
    
    try:
        compare_v4c(args.inputs, args.ylim, args.scale, args.output, args.dpi)
    except V4CError as e:
        print(f"Error: {str(e)}")
        exit(1)
    except Exception as e:
        print(f"Unexpected error: {str(e)}")
        exit(1)

if __name__ == "__main__":
    main()
