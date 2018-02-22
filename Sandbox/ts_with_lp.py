import sys
sys.path.append("/home/shawn/DIII-D/ORNL-Fusion/Collector-Probes")

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot
import openpyxl as xl
import MDSplus as mds
import scipy.interpolate as scinter
import matplotlib.pyplot as plt
import get_ts as ts

# Fucntion for getting cells and putting into normal numpy array.
def returnArray(sheet, lowRange, highRange):
    cells = sheet[lowRange:highRange]
    cells = np.transpose(cells)
    cells = np.reshape(cells, cells.size)
    values = np.array([cell.value for cell in cells])
    return values

# Open up the LP Excel sheet.
wb = xl.load_workbook("link_to_lp_with_fits.xlsx", data_only=True)
sheet = wb.get_sheet_by_name("LP Data")

# We will use 167192 plunge #2 since it goes a little further in.
lp_r  = returnArray(sheet, "K3", "K55") / 100.0 - 0.016   # cm to m
lp_te = returnArray(sheet, "M3", "M55")

if True:
    # Now lets map these R's to psin's, so we can compare with the TS.
    # Grab time 2960 since it's about the middle time for the LP data.
    MDSconn = mds.Connection("localhost")
    gfile = ts.load_gfile_mds(167192, 2960, connection=MDSconn)

    # Create grid of R's and Z's.
    Rs, Zs = np.meshgrid(gfile['R'], gfile['Z'])

    # Coordinates of magnetic axis.
    Z_axis = gfile['ZmAxis']
    R_axis = gfile['RmAxis']

    # Get R's and Z's of the lcfs. Just the right half of it.
    Zes = np.copy(gfile['lcfs'][:, 1][13:-12])
    Res = np.copy(gfile['lcfs'][:, 0][13:-12])

    # Interpolate to create function where you give a Z, and it gives you
    # the R of the lcfs.
    f_Rs = scinter.interp1d(Zes, Res, assume_sorted=False)

    # Only R's on the right half.
    Rs_trunc = Rs > R_axis

    # Interpolation functions of psin(R, Z) and R(psin, Z).
    f_psin = scinter.Rbf(Rs[Rs_trunc], Zs[Rs_trunc], gfile['psiRZn'][Rs_trunc])

    # No we have the function we want, psin(R, Z). Let's get those psins.
    lp_psins = np.array([f_psin(r, -0.18) for r in lp_r])

    # Plot it just to check it looks good.
    #plt.plot(lp_psins, lp_te)
    #plt.show()

    # Great looks good. Let's get the TS data now.
    ts_dict = ts.run_script(167192, "core", tmin=2930, tmax=2990, tstep=20)
    ts_psins = ts_dict["psins"]["avg_psins"]
    ts_te = ts_dict["psins"]["avg_Tes"]

    def exp_fit(x, a, b, c):
        return a * np.exp(-b * (x-1.0)) + c

    last_sol_index = np.argmax(ts_psins < 1.0)
    popt, pcov = curve_fit(exp_fit, ts_psins[:last_sol_index], ts_te[:last_sol_index], maxfev=5000)
    fit_to_ts = exp_fit(np.arange(0.9, 1.5, 0.01), *popt)

    # And now let's plot them together.
    plt.rcParams.update({'font.size': 22})
    plt.plot(lp_psins, lp_te, label="Langmuir")
    plt.plot(ts_psins, ts_te, ".", ms=15.0, label="Thomson")
    plt.plot(np.arange(0.9, 1.5, 0.01), fit_to_ts, "--", label="Thomson Fit")
    plt.title("Shot 167192 LP and TS Data")
    plt.legend()
    plt.xlabel(r"$\mathrm{\phi_n}$")
    plt.ylabel("Te (eV)")
    plt.xlim([1.0, 1.4])
    plt.ylim([0, 40])
    plt.show()

# Let's do a fit to all the LP data too.
sheet2 = wb.get_sheet_by_name("Fit Te's")
all_lp_rmin = returnArray(sheet2, "N2", "N348")
all_lp_te = returnArray(sheet2, "B2", "B348")

def exp_fit(x, a, b, c):
    return a * np.exp(-b * x) + c

popt, pcov = curve_fit(exp_fit, all_lp_rmin, all_lp_te, maxfev=5000)
fit_to_lp = exp_fit(np.arange(2.1, 15.8, 0.01), *popt)

plt.plot(all_lp_rmin, all_lp_te, ".", ms=8.0, label="Langmuir")
plt.plot(np.arange(2.1, 15.8, 0.01), fit_to_lp, "--", label="Fit to Langmuir")
plt.xlabel(r"$\mathrm{R-R_{sep} (cm)}$")
plt.ylabel("Te (eV)")
plt.title("Plunging Langmuir Probe")
plt.legend()
plt.show()
