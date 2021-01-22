import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from gadata import gadata


def get_fs_data(shot, tag="FS04"):
    """
    """

    # Load the gadata object. Assuming localhost is linked to atlas. See README
    # for instructions on how to do that.
    ga_obj = gadata(tag, shot)

    # Just return the time, signal.
    return ga_obj.xdata, ga_obj.zdata


def calc_freq(time, signal, tag="FS04", time_window=[2000, 4500], distance=200, height=1e16, verbal=True):
    """
    """

    # Restrict to data in time window.
    window = np.where(np.logical_and(time>=time_window[0], time<=time_window[1]))[0]
    time   = time[window]
    signal = signal[window]

    # Use scipy function. distance is how many points from a peak we go before
    # looking for the next peak. height is the required height of the peak to
    # be considered. This variables change dpeending on shot, so trial and
    # error is likely required.
    peaks = find_peaks(signal, distance=distance, height=height)[0]

    # Plot to see if peaks are being detected corrrectly.
    fig, ax = plt.subplots()
    ax.plot(time, signal, linestyle='-', color='k')
    ax.plot(time[peaks], signal[peaks], lw=0, marker="*", color="r")
    ax.set_xlabel("Time (ms)", fontsize=16)
    ax.set_ylabel(tag, fontsize=16)
    fig.tight_layout()
    fig.show()

    if verbal:
        print("Zoom in on plot to ensure only ELM peaks are being detected.")
        print("Tips:")
        print(" -Use different height (currently {:}) such that detection".format(height))
        print("  height sufficiently separates ELMs from the base signal.")
        print(" -Change distance (currently {:}), which is the amount of ".format(distance))
        print("  points surrounding the peak to ignore. Set this such that once")
        print("  you go distance amount of points away from the peak, you are below ")
        print("  what is set for height.\n")

    # Now calculate the frequency. Simply number of peaks divided by time width.
    freq = len(peaks) / (time_window[1] - time_window[0]) / 1000
    print("ELM frequency between {:}-{:} ms: {:.2e} Hz".format(time_window[0], time_window[1], freq))
