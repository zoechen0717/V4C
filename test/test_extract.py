import v4c

def test_extract():
    mcool_files = ["test.mcool"]
    resolutions = [5000]
    coords = "chr17:45878152-46000000"
    output_file = "test_output.tsv"

    v4c.extract_v4c(mcool_files, resolutions, coords, genome="hg38", flank=50000, balance=True, scale=True, output=output_file)

    with open(output_file, "r") as f:
        lines = f.readlines()
    assert len(lines) > 1  # Ensure output is not empty
