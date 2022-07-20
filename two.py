from pathlib import Path
import matplotlib.pyplot as plt
import fathon
from fathon import fathonUtils as fu
import numpy as np


def get_ycoords(sequence):
    """Get y coordinates from sequence for the time series"""

    ycoords = [0]
    for element in sequence:
        try:
            if element == 'A' or element == 'G' or element == 'a' or element == 'g':
                ycoords.append(ycoords[-1] + 1)
            elif element == 'C' or element == 'T' or element == 'c' or element == 't':
                ycoords.append(ycoords[-1] - 1)
        except Exception:
            raise
    return ycoords


def plot_time_series(ycoords, filename):
    """
    Generate plot of time series from a list of y coordinates

    :param ycoords:
    :return:
    """

    # Plot time series
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.cla()

    p, = ax.plot(ycoords)
    plt.title(f"{filename} Time Series")
    plt.grid()
    fig.savefig(f"{filename}.png")
    p.remove()

    a = ycoords  # What the hell is this ahah why not just keep ycoords

    # Perform DFA analysis

    a = fu.toAggregated(a)  # zero-mean cumulative sum

    # initialize dfa object
    pydfa = fathon.DFA(a)

    # Compute fluctuation function and Hurst exponent
    wins = fu.linRangeByStep(10, 2000)
    n, F = pydfa.computeFlucVec(wins, revSeg=True,
                                polOrd=3)  # What is n and F and why don't you ever use them again?
    H, H_intercept = pydfa.fitFlucVec()

    # Compute Hurst exponent in different ranges
    limits_list = np.array([[15, 2000], [200, 1000]], dtype=int)
    list_H, list_H_intercept = pydfa.multiFitFlucVec(
        limits_list)  # Same question here
    print(f'Hurst exponent: {H} , H intercept: {H_intercept}')

    # Plot Hurst curve
    x_axis = np.linspace(-10, 10, 20)
    y = H * x_axis + H_intercept

    ax.cla()
    p, = ax.plot(x_axis, y)
    plt.title(f"{filename} Hurst curve")
    plt.grid()

    # Save result as PNG
    fig.savefig(filename + "_hurst_curve.png")
    p.remove()


def read_genome(path_to_directory):
    """
    This function reads the input FASTA file and converts
    the data into a time series to perform a detrended
    fluctuational analysis (DFA) of the input data.

    :param path_to_directory:
    :return:
    """

    directory = Path(path_to_directory)
    # Get paths of FASTA files to read
    file_paths = [path for path in directory.glob("*.fasta")]
    for path in file_paths:
        print(f"File being read: {path}\n")
        with open(path, "r") as f:
            # Isolate nucleotide sequence
            genome = [line for line in f if line[0] != ">"]
            genome = "".join(genome)

        ycoords = get_ycoords(genome)
        plot_time_series(ycoords, path)


read_genome(r"C:\Users\legregam\Documents\Projets\Swurl\files")
