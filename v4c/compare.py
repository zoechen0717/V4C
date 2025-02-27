def compare_v4c(input_files, ylim=1, scale=True):
    """
    Compares Virtual 4C data from multiple .mcool files.

    Parameters:
    - input_files (list): List of TSV files containing extracted V4C data.
    - ylim (float): Maximum y-axis value.
    - scale (bool): Whether to normalize values between 0 and 1.
    """

    plt.figure(figsize=(10, 6))

    for file in input_files:
        df = pd.read_csv(file, sep="\t")

        for _, row in df.iterrows():
            coords = list(map(float, row[5:].values))
            if scale:
                min_val, max_val = min(coords), max(coords)
                coords = [(x - min_val) / (max_val - min_val) if max_val > min_val else 0 for x in coords]

            plt.plot(coords, label=row['mcool'])

    plt.xlabel("Genomic Position (bp)")
    plt.ylabel("Hi-C Contact Frequency")
    plt.title("Virtual 4C Comparison")
    plt.legend()
    plt.ylim(0, ylim)
    plt.show()
