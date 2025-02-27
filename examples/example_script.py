# example_script.py
import v4c

# Define inputs
mcool_files = ["sample1.mcool", "sample2.mcool"]
resolutions = [5000, 10000]
coords = "chr17:45878152-46000000"
genome = "hg38"
flank = 50000
balance = True
scale = True
output_file = "extracted_data.tsv"

# Step 1: Extract V4C data
print("Extracting Hi-C contact frequencies...")
v4c.extract_v4c(mcool_files, resolutions, coords, genome, None, flank, balance, scale, output_file)

# Step 2: Plot the extracted data
print("Generating Virtual 4C plot...")
v4c.plot_v4c(output_file, ylim=0.4, flank=50000)

# Step 3: Compare multiple datasets
print("Comparing multiple datasets...")
v4c.compare_v4c(["extracted_data1.tsv", "extracted_data2.tsv"], ylim=0.4, sca
