from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import openpyxl as xl
import MDSplus as mds
from scipy.optimize import curve_fit

import sys
sys.path.append("/home/shawn/DIII-D/d3d_tools")
import r_to_psin


# Take an Excel sheet and a range of cells (i.e. "H3":"H64") and return an array of the values.
def putIntoArray(sheet, minRange, maxRange):
    cells = sheet[minRange:maxRange]
    cells = np.array(cells)
    cells = np.reshape(cells, cells.size)
    values = [cell.value for cell in cells]
    values = np.transpose(values)
    return values


# Open the plunging LP workbook and grab the sheet with all the data.
wb = xl.load_workbook("recLPdata.xlsx", data_only=True)
sheet = wb.get_sheet_by_name("Sheet1")
LPradii = putIntoArray(sheet, "J3", "J55")
LPTes = putIntoArray(sheet, "L3", "L55")
LPdens = putIntoArray(sheet, "K3", "K55")
LPtimes = putIntoArray(sheet, "I3", "I55")
# Convert to radii from cm to m.
LPradii = [radius/100.0 for radius in LPradii]

LPpsins = np.zeros(len(LPradii))
MDSplusConn = mds.Connection("localhost")
# Convert LP data to psi_n.
for index in range(0, len(LPradii)):
    LPradius = LPradii[index]
    LPtime = LPtimes[index]
    tmp_psin = r_to_psin.RtoPsin(radius=LPradius, z_loc=-0.18, shot=167192,
                                 time=LPtime, MDSplusConn=MDSplusConn)
    LPpsins[index] = tmp_psin


filename = "/home/shawn/DIII-D/OMFITprofiles_data/TS_167192.nc"
rootgrp = Dataset(filename, "r", format="NETCDF4")

if False:
    print("\nVariables in file:")
    for key in rootgrp.variables.keys():
        print("  " + key)

# 0-39 are the core. 40-47 are divertor. 48-53 are tangential.
psi_n = np.array(rootgrp["psi_n"][:39])
T_e = np.array(rootgrp["T_e"][:39])
n_e = np.array(rootgrp["n_e"][:39])

# Reorganize these variables into easy to plot (i.e. value vs. psi_n) arrays.
org_psi_n = np.zeros((psi_n.shape[1], psi_n.shape[0]))
org_T_e = np.zeros((T_e.shape[1], T_e.shape[0]))
org_n_e = np.zeros((n_e.shape[1], n_e.shape[0]))

for value_at_time in range(0, psi_n.shape[1]):
    for chord in range(0, psi_n.shape[0]):
        tmp_value = psi_n[chord][value_at_time]
        org_psi_n[value_at_time][chord] = tmp_value

for value_at_time in range(0, T_e.shape[1]):
    for chord in range(0, T_e.shape[0]):
        tmp_value = T_e[chord][value_at_time]
        org_T_e[value_at_time][chord] = tmp_value

for value_at_time in range(0, n_e.shape[1]):
    for chord in range(0, n_e.shape[0]):
        tmp_value = n_e[chord][value_at_time]
        org_n_e[value_at_time][chord] = tmp_value

for index in range(0, len(org_psi_n)):
    print("Index: " + str(index))
    arr1d_org_psi_n = org_psi_n[index]
    arr1d_org_T_e = org_T_e[index]
    arr1d_org_n_e = org_n_e[index]
    idx = np.argsort(arr1d_org_psi_n)
    arr1d_org_psi_n = arr1d_org_psi_n[idx]
    arr1d_org_T_e = arr1d_org_T_e[idx]
    arr1d_org_n_e = arr1d_org_n_e[idx]
    org_psi_n[index] = arr1d_org_psi_n
    org_T_e[index] = arr1d_org_T_e
    org_n_e[index] = arr1d_org_n_e

# Create array of the average values at each psin from TS.
# Each chord remains at about the same psin value, so use avg values from each chord.
avg_TS_psins = np.array([])
avg_TS_temps = np.array([])
avg_TS_ne = np.array([])
for chord in range(0,psi_n.shape[0]):
    avg_psin_at_chord = np.mean(psi_n[chord])
    avg_temp_at_chord = np.mean(T_e[chord])
    avg_ne_at_chord = np.mean(n_e[chord])
    avg_TS_psins = np.append(avg_TS_psins, avg_psin_at_chord)
    avg_TS_temps = np.append(avg_TS_temps, avg_temp_at_chord)
    avg_TS_ne = np.append(avg_TS_ne, avg_ne_at_chord)

plt.plot(avg_TS_psins, avg_TS_temps, label="Averages")

# Now let's do an exp. fit to the avg data. First find first index where psin > 1.0.

# This will return the index where avg_TS_psins first goes below 1.0. For ex, if
# the array has [1.06, 1.03, 0.98, ...] it will return the index of 0.98, which is 2
# in this case. Start exp fit from this point.
first_sol_index = np.argmax(avg_TS_psins<1.0)

# Fit needs to be shifted to the start of our SOL.
def exp_fit(x, a, b, c):
    return a * np.exp(-b * (x-avg_TS_psins[first_sol_index])) + c

popt1, pcov1 = curve_fit(exp_fit, avg_TS_psins[:first_sol_index], avg_TS_temps[:first_sol_index], maxfev=5000)
plt.plot(np.arange(1.0, 1.4, 0.01), exp_fit(np.arange(1.0, 1.4, 0.01), *popt1), "g", label="Exp. Fit")

#for plot_index in range(1, len(org_psi_n[1])):
#    psi_arr = org_psi_n[plot_index]
#    T_e_arr = org_T_e[plot_index]
#    print("Plot #" + str(plot_index))
#    plt.plot(psi_arr, T_e_arr, "r", alpha=0.025, linewidth=10)

for plot_index in range(1, len(org_psi_n[1])):
    psi_arr = org_psi_n[plot_index]
    T_e_arr = org_T_e[plot_index]
    plt.plot(psi_arr, T_e_arr, "r", alpha=0.05, linewidth=10)

plt.plot(LPpsins, LPTes, label="LP")

plt.xlim([1.0, 1.4])
#plt.xlim([1.0, 1.075])
plt.ylim([0,50])
plt.ylabel(r"$\mathrm{T_e\ (eV)}$")
plt.xlabel(r"$\mathrm{\psi_N}$")
plt.title("TS Data vs. Plunging LP Shot #167192")
plt.legend()
plt.show()

# Do the same but for densities.
for index in range(0, len(avg_TS_ne)):
    avg_TS_ne[index] = avg_TS_ne[index] * 10**(-18)

plt.plot(avg_TS_psins, avg_TS_ne, label="Averages")
first_sol_index = np.argmax(avg_TS_psins<1.0)
guess = (1, 1, 1)
popt2, pcov2 = curve_fit(exp_fit, avg_TS_psins[:first_sol_index], avg_TS_ne[:first_sol_index], guess, maxfev=5000)
plt.plot(np.arange(1.0, 1.4, 0.01), exp_fit(np.arange(1.0, 1.4, 0.01), *popt2), "g", label="Exp. Fit")
#plt.plot(np.arange(0, 1.4, 0.01), exp_fit(np.arange(0, 1.4, 0.01), *popt2), "g", label="Exp. Fit")
for plot_index in range(1, len(org_psi_n[1])):
    psi_arr = org_psi_n[plot_index]
    n_e_arr = org_n_e[plot_index] * 10**(-18)
    plt.plot(psi_arr, n_e_arr, "r", alpha=0.05, linewidth=10)
plt.plot(LPpsins[:len(LPdens)], LPdens, label="LP")
plt.xlim([1.0, 1.4])
plt.ylim([0, 15])
plt.ylabel(r"$\mathrm{n_e\ (10^{18}\ m^{-3})}$")
plt.xlabel(r"$\mathrm{\psi_N}$")
plt.title("TS Data vs. Plunging LP Shot #167192")
plt.legend()
plt.show()
