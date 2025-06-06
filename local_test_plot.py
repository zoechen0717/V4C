import os
import matplotlib.pyplot as plt
from v4c import plot_v4c

def test_plot_v4c():
    """
    Local test for plot_v4c function.
    """
    input_file = "test_output.tsv"
    ylim = 1
    flank = 50000

    # Run plotting function
    plot_v4c(input_file, ylim, flank)

    # Check if plot was generated (Matplotlib does not create files by default, so this checks for execution errors)
    print("Local test for plot_v4c executed successfully!")

if __name__ == "__main__":
    test_plot_v4c()
