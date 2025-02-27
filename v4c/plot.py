import matplotlib.pyplot as plt
import pandas as pd

def plot_v4c(input_file, ylim=0.4, flank=50000):
    """
    Plots Virtual 4C contact frequencies.

    Parameters:
    - input_file (str): TSV file containing extracted V4C data.
    - ylim (float): Maximum y-axis value.
    - flank (int): Flanking region (bp) to display.
    """

    df = pd.read_csv(input_file, sep="\t")

    for mcool_file in df["mcool"].unique():
        plt.figure(figsize=(10, 6))
        subset = df[df["mcool"] == mcool_file]

        for _, row in subset.iterrows():
            coords = list(map(float, row[5:].values))
            plt.plot(coords, label=f"Resolution {row['res']}")

        plt.xlabel("Genomic Position (bp)")
        plt.ylabel("Hi-C Contact Frequency")
        plt.title(f"Virtual 4C - {mcool_file}")
        plt.legend()
        plt.ylim(0, ylim)
        plt.show()
