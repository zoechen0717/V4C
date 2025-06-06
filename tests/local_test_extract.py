import os
import pandas as pd
from v4c import extract_v4c

def test_extract_v4c():
    """
    Local test for extract_v4c function.
    """
    mcool_files = ["/Users/zoechen/Documents/Work/Corces/mapq30/mcool/Oligo_merged_chr17.genome1.balanced.mapq30.dedup.mcool", "/Users/zoechen/Documents/Work/Corces/mapq30/mcool/Oligo_merged_chr17.genome2.mapq30.dedup.mcool"]
    resolutions = [5000, 10000]
    coords = None
    genes = "MAPT,CRHR1"
    genome = "hg38"
    bed_file = None
    flank = 50000
    balance = True
    scale = True
    output_file = "test_output.tsv"

    # Run extraction
    extract_v4c(mcool_files, resolutions, coords, genes, genome, bed_file, flank, balance, scale, output_file)

    # Check output
    assert os.path.exists(output_file), "Output file was not created."
    df = pd.read_csv(output_file, sep="\t")
    assert not df.empty, "Output file is empty."
    print("Local test passed successfully!")

if __name__ == "__main__":
    test_extract_v4c()
