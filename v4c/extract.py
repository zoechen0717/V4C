import cooler
import numpy as np
import pandas as pd
import argparse
import os
from typing import List, Optional, Union
from .utils import get_promoter_coords, validate_mcool_files
import sys

class V4CError(Exception):
    """Base exception class for V4C errors"""
    pass

class InputValidationError(V4CError):
    """Exception raised for input validation errors"""
    pass

class FileProcessingError(V4CError):
    """Exception raised for file processing errors"""
    pass

def validate_inputs(mcool_files: List[str], 
                   resolutions: List[int], 
                   coords: Optional[str], 
                   genes: Optional[str], 
                   genome: Optional[str], 
                   bed_file: Optional[str],
                   normalization_method: str = "minmax") -> None:
    """
    Validates input parameters for extract_v4c function.
    
    Args:
        mcool_files: List of .mcool files
        resolutions: List of resolutions
        coords: Genomic coordinates
        genes: Comma-separated gene names
        genome: Reference genome version
        bed_file: Path to BED file
        normalization_method: Normalization method ("minmax" or "self")
        
    Raises:
        InputValidationError: If any input parameter is invalid
    """
    # Validate mcool files
    if not mcool_files:
        raise InputValidationError("No mcool files provided")
    
    # Validate resolutions
    if not resolutions:
        raise InputValidationError("No resolutions provided")
    if not all(isinstance(r, int) and r > 0 for r in resolutions):
        raise InputValidationError("Resolutions must be positive integers")
    
    # Validate genome version
    if genome and genome not in ["hg38", "hg19"]:
        raise InputValidationError(f"Invalid genome version: {genome}. Must be 'hg38' or 'hg19'")
    
    # Validate coordinates format
    if coords:
        try:
            chrom, pos = coords.split(":")
            start, end = map(int, pos.split("-"))
            if start >= end:
                raise InputValidationError("Start position must be less than end position")
        except ValueError:
            raise InputValidationError("Invalid coordinates format. Expected format: 'chr:start-end'")
    
    # Validate bed file
    if bed_file and not os.path.exists(bed_file):
        raise InputValidationError(f"BED file not found: {bed_file}")
    
    # Validate normalization method
    if normalization_method not in ["minmax", "self"]:
        raise InputValidationError(f"Invalid normalization method: {normalization_method}. Must be 'minmax' or 'self'")

def extract_v4c(mcool_files: List[str], 
                resolutions: List[int], 
                coords: Optional[str] = None, 
                genes: Optional[str] = None, 
                genome: Optional[str] = None, 
                bed_file: Optional[str] = None, 
                flank: int = 500000, 
                balance: bool = True, 
                scale: bool = True, 
                normalization_method: str = "minmax",
                use_fixed_center: bool = False,
                output: str = "extracted_data.tsv") -> None:
    """
    Extracts Virtual 4C contact frequencies from .mcool files at specific resolutions.

    Args:
        mcool_files: List of .mcool files
        resolutions: List of resolutions to extract
        coords: Genomic coordinates (e.g., "chr17:45878152-46000000")
        genes: Comma-separated gene names to extract promoter coordinates
        genome: Reference genome ("hg38" or "hg19")
        bed_file: Custom BED file for genomic coordinates
        flank: Number of base pairs to extend upstream and downstream (default: 500000)
        balance: Whether to use ICE balancing
        scale: Whether to normalize values between 0 and 1
        normalization_method: Normalization method - "minmax" or "self" (default: "minmax")
        use_fixed_center: Whether to use fixed center position calculation like original code (default: False)
        output: Output file name

    Returns:
        None: Saves extracted contact frequencies to a TSV file

    Raises:
        InputValidationError: If input parameters are invalid
        FileProcessingError: If file processing fails
        V4CError: For other V4C-related errors
    """
    try:
        # Validate inputs
        validate_inputs(mcool_files, resolutions, coords, genes, genome, bed_file, normalization_method)
        
        # Validate mcool files exist
        valid_mcool_files = validate_mcool_files(mcool_files)
        if not valid_mcool_files:
            raise FileProcessingError("No valid mcool files found")

        promoter_coords = []
        gene_mapping = {}  # Map coordinates to gene names

        if genes and genome:
            gene_list = genes.split(",")  # Split comma-separated gene names
            for gene in gene_list:
                gene = gene.strip()
                coords = get_promoter_coords(genome, gene=gene)
                promoter_coords.extend(coords)
                # Map each coordinate to the gene name
                for coord in coords:
                    gene_mapping[tuple(coord)] = gene
        elif genome and coords:
            chrom, start, end = coords.split(":")[0], int(coords.split(":")[1].split("-")[0]), int(coords.split(":")[1].split("-")[1])
            # For coordinate input, use the center point as the TSS
            center = (start + end) // 2
            promoter_coords = [(chrom, center, center)]
        elif bed_file:
            try:
                # Read BED file and determine number of columns
                bed_df = pd.read_csv(bed_file, sep="\t", header=None)
                
                if len(bed_df.columns) >= 4:
                    # BED file has 4+ columns (chrom, start, end, gene_name, ...)
                    bed_df.columns = ["chrom", "start", "end", "gene_name"] + [f"col_{i}" for i in range(4, len(bed_df.columns))]
                    promoter_coords = bed_df[["chrom", "start", "end"]].values.tolist()
                    # Map coordinates to gene names
                    for i, row in bed_df.iterrows():
                        coord_tuple = (row["chrom"], row["start"], row["end"])
                        gene_mapping[coord_tuple] = row["gene_name"]
                else:
                    # BED file has only 3 columns (chrom, start, end)
                    bed_df.columns = ["chrom", "start", "end"]
                    promoter_coords = bed_df[["chrom", "start", "end"]].values.tolist()
                    # No gene names available, leave gene_mapping empty
                    
            except Exception as e:
                raise FileProcessingError(f"Error reading BED file: {str(e)}")
        else:
            raise InputValidationError("Either --coords, --genes with --genome, or --bed must be specified.")

        results = []

        for mcool in valid_mcool_files:
            for res in resolutions:
                try:
                    c = cooler.Cooler(f"{mcool}::/resolutions/{res}")
                    for chrom, start, end in promoter_coords:
                        # Calculate region boundaries
                        if use_fixed_center:
                            # Use original code logic: fixed 500kb flanking from TSS
                            row_bin_start = (start // res) * res
                            region_start = max(0, row_bin_start - flank)
                            region_end = row_bin_start + flank
                        else:
                            # Use current logic: flanking from TSS region
                            region_start = max(0, start - flank)
                            region_end = end + flank
                        
                        matrix = c.matrix(balance=balance, sparse=False).fetch((chrom, region_start, region_end))
                        np.nan_to_num(matrix, copy=False)

                        # Calculate row index
                        if use_fixed_center:
                            # Use original code logic: fixed center position
                            row_index = (flank // res)
                        else:
                            # Use current logic: relative to TSS position
                            row_index = (start // res) - (region_start // res)
                        
                        row_values = matrix[row_index, :].tolist()

                        # Apply normalization
                        if scale:
                            if normalization_method == "self":
                                # Original code self-normalization
                                if row_values[row_index] > 0:
                                    row_values = [x / row_values[row_index] for x in row_values]
                                else:
                                    row_values = [0] * len(row_values)
                            else:  # minmax normalization
                                min_val, max_val = np.min(row_values), np.max(row_values)
                                row_values = [(x - min_val) / (max_val - min_val) if max_val > min_val else 0 for x in row_values]

                        genomic_coords = np.arange(region_start, region_end, res)
                        # Add gene name if available
                        gene_name = gene_mapping.get((chrom, start, end), "")
                        results.append([mcool, res, chrom, start, end, gene_name] + row_values)
                except Exception as e:
                    raise FileProcessingError(f"Error processing file {mcool} at resolution {res}: {str(e)}")

        try:
            df = pd.DataFrame(results)
            # Add column names
            if len(results) > 0:
                num_contact_cols = len(results[0]) - 6  # 6 metadata columns (added gene_name)
                # Generate genomic coordinate column names
                # Get the first row to determine the genomic region
                first_row = results[0]
                chrom, start, end, resolution = first_row[2], first_row[3], first_row[4], first_row[1]
                
                # Calculate the genomic coordinates for each contact column
                contact_cols = []
                for i in range(num_contact_cols):
                    # Calculate the genomic position for this contact bin
                    # The contact data spans from (start - flank) to (end + flank)
                    total_region_start = start - flank
                    bin_position = total_region_start + (i * resolution)
                    contact_cols.append(str(bin_position))
                
                column_names = ["mcool", "res", "chrom", "start", "end", "gene_name"] + contact_cols
                # Ensure column names match DataFrame columns
                if len(column_names) == len(df.columns):
                    df.columns = column_names
                else:
                    # Fallback: use proper column names
                    df.columns = ["mcool", "res", "chrom", "start", "end", "gene_name"] + [f"contact_{i}" for i in range(len(df.columns) - 6)]
            df.to_csv(output, sep="\t", index=False)
        except Exception as e:
            raise FileProcessingError(f"Error saving results to {output}: {str(e)}")

    except V4CError:
        raise
    except Exception as e:
        raise V4CError(f"Unexpected error in extract_v4c: {str(e)}")

def main():
    """Command line interface for extract_v4c function."""
    parser = argparse.ArgumentParser(description="Extract Virtual 4C contact frequencies from .mcool files")
    parser.add_argument("mcool_files", nargs="+", help="One or more .mcool files")
    parser.add_argument("--resolutions", nargs="+", type=int, required=True, help="List of resolutions to extract")
    parser.add_argument("--coords", help="Genomic coordinates (e.g., 'chr17:45878152-46000000')")
    parser.add_argument("--genes", help="Comma-separated gene names")
    parser.add_argument("--genome", choices=["hg38", "hg19"], help="Reference genome version")
    parser.add_argument("--bed", help="Path to BED file")
    parser.add_argument("--flank", type=int, default=500000, help="Number of base pairs to extend upstream and downstream (default: 500000)")
    parser.add_argument("--no-balance", action="store_true", help="Disable ICE balancing")
    parser.add_argument("--no-scale", action="store_true", help="Disable normalization between 0 and 1")
    parser.add_argument("--normalization", choices=["minmax", "self"], default="minmax", help="Normalization method: 'minmax' or 'self' (default: minmax)")
    parser.add_argument("--use-fixed-center", action="store_true", help="Use fixed center position calculation like original code")
    parser.add_argument("--output", default="extracted_data.tsv", help="Output file name")

    args = parser.parse_args()

    try:
        extract_v4c(
            mcool_files=args.mcool_files,
            resolutions=args.resolutions,
            coords=args.coords,
            genes=args.genes,
            genome=args.genome,
            bed_file=args.bed,
            flank=args.flank,
            balance=not args.no_balance,
            scale=not args.no_scale,
            normalization_method=args.normalization,
            use_fixed_center=args.use_fixed_center,
            output=args.output
        )
    except V4CError as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
