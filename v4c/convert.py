import os
import subprocess
import argparse

def convert_hic_to_mcool(hic_file, output_dir="output", resolutions=None):
    """
    Converts a .hic file to a .mcool file using hic2cool.

    Parameters:
    - hic_file (str): Path to the input .hic file.
    - output_dir (str): Directory to save the .mcool file.
    - resolutions (list, optional): List of resolutions to extract.

    Returns:
    - Path to the converted .mcool file.
    """
    if resolutions is None:
        resolutions = [1000, 5000, 10000, 25000, 50000, 100000]  # Default resolutions

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Define output file
    mcool_file = os.path.join(output_dir, os.path.basename(hic_file).replace(".hic", ".mcool"))

    # Build hic2cool command
    cmd = [
        "hic2cool", "convert",
        "--input", hic_file,
        "--output", mcool_file,
        "--resolutions"
    ] + list(map(str, resolutions))

    try:
        subprocess.run(cmd, check=True)
        print(f"Conversion complete: {mcool_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error during conversion: {e}")

    return mcool_file

# Command-line execution
def main():
    parser = argparse.ArgumentParser(description="Convert .hic to .mcool")
    parser.add_argument("--hic", type=str, required=True, help="Input .hic file")
    parser.add_argument("--output", type=str, default="output", help="Output directory")
    parser.add_argument("--resolutions", nargs="+", type=int, help="List of resolutions")

    args = parser.parse_args()
    convert_hic_to_mcool(args.hic, args.output, args.resolutions)

if __name__ == "__main__":
    main()
