import sys
# Will need to change this to the directory right above Sandbox for you.
sys.path.append('/home/shawn/DIII-D/ORNL-Fusion/Collector-Probes')

# For pointing to the netCDF files.
from Tkinter import Tk
from tkFileDialog import askopenfilename

# OMFITprofiles data in netCDF file format.
from netCDF4 import Dataset

import numpy as np
from scipy.optimize import curve_fit
from scipy import interpolate
import matplotlib.pyplot as plt

import get_lp as lp
import pull_data_dp as pull
import ProbeClass as Probe
import shawn_model as model


def getOMFITprofilesFiles():
    # You will need the netCDF4 file from OMFITprofiles for the TS data (one for each shot).
    print "Specify the shot and NetCDF file(s) from OMFITprofiles."
    Tk().withdraw()
    netcdf_files = {}
    while(True):
        shot = raw_input("Shot: ")
        filename = askopenfilename()
        netcdf_files[str(shot)] = filename
        ans = raw_input("Another file? (y/n): ")
        if ans is not "y":
            break

    return netcdf_files

def loadOMFITdata(shots, plotData=False):
    print "Make sure the OMFITprofiles files are stored in the current directory as \"OMFITprofiles_[shot]_FIT.nc\"."
    all_avg_te = {}
    for shot in shots:
        # Get the filename and open the file in rootgrp.
        filename = "OMFITprofiles_" + str(shot) + "_FIT.nc"
        rootgrp = Dataset(filename, "r", format="NETCDF4")

        # Get the netCDF variables.
        psin_var = rootgrp["psi_n"]
        te_var    = rootgrp["T_e"]
        times_var  = rootgrp["time"]

        # Get the actual data in psin and time variables.
        psin = psin_var[:]
        times = times_var[:]

        # Compute the average temp for each psin coordinate.
        bad_shots = 0
        te_sum = np.zeros(len(psin))
        for time_index in range(0, len(times)):
            # Get te profile at a time
            temp_prof = te_var[:][time_index]

            # Some shots have sketchy data, so skip them.
            if shot == 167405:
                bad_shots = 2
                if times[time_index] == 3700 or times[time_index] == 3800:
                    continue
            if shot == 167406:
                bad_shots = 1
                if times[time_index] == 3500:
                    continue
            if shot == 167407:
                bad_shots = 3
                if times[time_index] == 2800 or times[time_index] == 2900 or times[time_index] == 4200:
                    continue
            if shot == 167408:
                bad_shots = 4
                if times[time_index] == 3400 or times[time_index] == 4200 or times[time_index] == 4300 or times[time_index] == 4400:
                    continue

            for psin_index in range(0, len(psin)):
                # Get the te value at a psin location of the profile. Add to running sum for each location.
                temp_prof_loc = temp_prof[psin_index]
                te_sum[psin_index] = te_sum[psin_index] + temp_prof_loc

        avg_te = np.array([])
        for value in te_sum:
            # Make sure to subtract the amount of bad shots since they aren't in the sum.
            avg  = float(value) / float(len(times)-bad_shots)
            avg_te = np.append(avg_te, avg)

        all_avg_te[str(shot)] = avg_te
        if plotData:
            plt.plot(psin, avg_te, ":", label=str(shot))

    all_avg_te["psin"] = psin

    total_avg = np.zeros(len(psin))
    for index in range(0, len(psin)):
        tmp_sum_at_psin = 0
        for shot_key in all_avg_te.keys():
            shot_data = all_avg_te[shot_key]
            tmp_sum_at_psin = tmp_sum_at_psin + shot_data[index]
        total_avg[index] = tmp_sum_at_psin / float(len(shots))

    total = np.array([psin, total_avg])

    if plotData:
        plt.plot(psin, total_avg, label="Total Average")
        plt.xlim(1, 1.11)
        plt.ylim(0, 50)
        plt.legend()
        plt.title("TS Te Averaged Profile")
        plt.show()

    return total



def avgTeFromLPs(lp_dict, timeStart, timeEnd, plotData=False):

    all_psins = np.array([])
    all_temps = np.array([])
    for shot in lp_dict.keys():
        lp_set = lp_dict[shot]
        for probe in lp_set.keys():
            print "Working on " + probe + " in set " + shot + "..."
            working_lp = lp_set[probe]
            times = np.array(working_lp["time"])
            psins = np.array(working_lp["psin"])
            temps = np.array(working_lp["temp"])

            # Find indices that define range timeStart-timeEnd.
            tmp = timeStart
            data = True
            while(data):
                try:
                    indexStart, = np.where(times==tmp)[0]
                    break
                except:
                    tmp = tmp + 1
                    if tmp > (timeStart + 10):
                        data=False
                        print "No data for " + probe + "."

            tmp = timeEnd
            while(data):
                try:
                    indexEnd, = np.where(times==tmp)[0]
                    break
                except:
                    tmp = tmp - 1

            if(data):
                avg_psin = np.mean(psins[indexStart:indexEnd])
                avg_temp = np.mean(temps[indexStart:indexEnd])
                all_psins = np.append(all_psins, avg_psin)
                all_temps = np.append(all_temps, avg_temp)

    #plt.plot(all_psins, all_temps, ".")

    # Create exponential fit for LP data.
    def exp_fit(x, a, b, c):
        return a * np.exp(-b * x) + c

    # Put into one 2d array.
    psin_vs_temp = np.array([all_psins, all_temps])

    # Sort by order of the first column (the psins).
    psin_vs_temp = zip(*psin_vs_temp)
    psin_vs_temp.sort(key=lambda x: x[0])
    psin_vs_temp = zip(*psin_vs_temp)


    guess = (1, 1, 1)
    popt, pcov = curve_fit(exp_fit, psin_vs_temp[0], psin_vs_temp[1], guess, maxfev=6000)
    fit_temp = exp_fit(np.array(psin_vs_temp[0]), *popt)
    if plotData:
        plt.plot(psin_vs_temp[0], psin_vs_temp[1], ".")
        plt.plot(psin_vs_temp[0], fit_temp)
        plt.show()

    return psin_vs_temp


def runScript(aNumber=28, bNumber=10, cNumber=10, timeStart=2500, timeEnd=4500):
    # The TS data.
    #netcdf_files = getOMFITprofilesFiles()

    # Get the probe objects as a list.
    probes = Probe.get_multiple(aNumber=aNumber, bNumber=bNumber, cNumber=cNumber,
                                 startTime=timeStart, endTime=timeEnd)
    aProbe = probes[0]
    bProbe = probes[1]
    cProbe = probes[2]

    # The shots the probes were in for.
    shots = aProbe.r2d2DICT["shots"]

    # Load the TS data averaged over all the shots at each psin.
    ts_data = loadOMFITdata(shots)

    # Get the divertor LP data.
    lp_dict = {}
    for shot in shots:
        lp_dict[str(shot)] = lp.get_dict_of_lps(shot)

    lp_data = avgTeFromLPs(lp_dict, timeStart, timeEnd)

    # Get the psiN values for each probe.
    au_psiN = aProbe.atlasDICT["psiN_U"]
    ad_psiN = aProbe.atlasDICT["psiN_D"]

    # Estimate of time of shot.
    timeOfShot = 5000.0
    numberOfShots = 4

    # Rough estimate of the net flux as defined in "shawn_model".
    au_net = aProbe.r2d2DICT["w_areal_U"] / (timeOfShot * numberOfShots)
    ad_net = aProbe.r2d2DICT["w_areal_D"] / (timeOfShot * numberOfShots)

    fig, ax1 = plt.subplots()
    ax1.plot(lp_data[0], lp_data[1], 'r.', label="LP Data")
    ax1.plot(ts_data[0], ts_data[1], 'r-', label="TS Data")
    ax1.set_ylabel(r"$\mathrm{T_e}$ (eV)", color="r")
    ax1.set_xlim([1.0, 1.4])
    ax1.set_ylim([0, 50])
    ax1.tick_params("x", colors="r")

    ax2 = ax1.twinx()
    ax2.plot(au_psiN, au_net, 'b', label="AU Net")
    ax1.set_xlim([1.0, 1.4])
    ax2.set_xlabel(r"$\psi_{\mathrm{N}}$")
    ax2.set_ylabel(r"Flux $\mathrm{m^{-2}s^{-1}}$", color="b")
    ax2.tick_params("y", colors="b")

    #plt.plot(au_psiN, au_net, label="AU Net")
    #plt.plot(lp_data[0], lp_data[1], label="LP Data")
    #plt.plot(ts_data[0], ts_data[1], label="TS Data")
    plt.title("TS/LP Te Data with Net Flux")
    plt.legend()
    plt.show()

    return [lp_data, ts_data, au_psiN, au_net]

    # Get the loss fluxes from the model.
    #au_loss = model.loss_flux(#Need parameters extrapolated from TS and/or divertor LP)
