import matplotlib.pyplot as plt
import pandas as pd
import os
from typing import Optional, Union
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
            figsize: tuple = (10, 6)) -> None:
    """
    Plots Virtual 4C contact frequencies.

    Args:
        input_file: TSV file containing extracted V4C data
        ylim: Maximum y-axis value
        flank: Flanking region (bp) to display
        output_file: Optional path to save the plot
        dpi: DPI for the output figure
        figsize: Figure size (width, height) in inches

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
        except Exception as e:
            raise FileProcessingError(f"Error reading input file: {str(e)}")

        # Create plots for each mcool file
        for mcool_file in df["mcool"].unique():
            plt.figure(figsize=figsize)
            subset = df[df["mcool"] == mcool_file]
            
            # Plot each resolution
            for _, row in subset.iterrows():
                try:
                    coords = list(map(float, row[5:].values))
                    plt.plot(coords, 
                            label=f"Resolution {row['res']}",
                            alpha=0.7)
                except Exception as e:
                    raise FileProcessingError(f"Error plotting data for {mcool_file}: {str(e)}")

            # Customize plot
            plt.xlabel("Genomic Position (bp)", fontsize=12)
            plt.ylabel("Hi-C Contact Frequency", fontsize=12)
            plt.title(f"Virtual 4C - {os.path.basename(mcool_file)}", fontsize=14)
            plt.legend()
            plt.ylim(0, ylim)
            plt.grid(True, alpha=0.3)
            
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
    
    args = parser.parse_args()
    
    try:
        plot_v4c(args.input, args.ylim, args.flank, args.output, args.dpi)
    except V4CError as e:
        print(f"Error: {str(e)}")
        exit(1)
    except Exception as e:
        print(f"Unexpected error: {str(e)}")
        exit(1)

if __name__ == "__main__":
    main()
