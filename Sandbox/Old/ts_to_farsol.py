from netCDF4 import Dataset
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def netcdf_to_dict(filename):
    """
    Will return a python dictionary of the psin, ne and Te data from a netCDF4 file.

    filename: The absolute path to the netCDF file.
    """

    rootgrp = Dataset(filename, "r", format="NETCDF4")

    # 0-39 are the core. 40-47 are divertor. 48-53 are tangential.
    psi_n = np.array(rootgrp["psi_n"][:39])
    T_e = np.array(rootgrp["T_e"][:39])
    n_e = np.array(rootgrp["n_e"][:39])

    # Reorganize these variables into easy to plot (i.e. value vs. psi_n) arrays.
    # May be unnecessary here, but could come in hand in later updates, so I include it.
    org_psi_n = np.zeros((psi_n.shape[1], psi_n.shape[0]))
    org_T_e = np.zeros((T_e.shape[1], T_e.shape[0]))
    org_n_e = np.zeros((n_e.shape[1], n_e.shape[0]))

    # Organize psins...
    for value_at_time in range(0, psi_n.shape[1]):
        for chord in range(0, psi_n.shape[0]):
            tmp_value = psi_n[chord][value_at_time]
            org_psi_n[value_at_time][chord] = tmp_value

    # Organize Te...
    for value_at_time in range(0, T_e.shape[1]):
        for chord in range(0, T_e.shape[0]):
            tmp_value = T_e[chord][value_at_time]
            org_T_e[value_at_time][chord] = tmp_value

    # Organize ne...
    for value_at_time in range(0, n_e.shape[1]):
        for chord in range(0, n_e.shape[0]):
            tmp_value = n_e[chord][value_at_time]
            org_n_e[value_at_time][chord] = tmp_value

    # Then sort them according to psin low->high.
    for index in range(0, len(org_psi_n)):
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

    # Put into dictionary and return.
    sliced_ts_dict = {}
    sliced_ts_dict["psin"] = org_psi_n
    sliced_ts_dict["Te"]   = org_T_e
    sliced_ts_dict["ne"]   = org_n_e
    sliced_ts_dict["psin_avg"] = avg_TS_psins
    sliced_ts_dict["Te_avg"]   = avg_TS_temps
    sliced_ts_dict["ne_avg"]   = avg_TS_ne

    return sliced_ts_dict

def fit_from_avg(sliced_ts_dict, min_psin=1.0, max_psin=1.4):
    """
    After running netcdf_to_dict, passing this that returned dict will attempt
    to fit an exponential to the data outside the LCFS.

    sliced_ts_dict: Dictionary returned from netcdf_to_dict.
    min_psin:       Where the expoenential fit will start. Default is LCFS.
    max_psin:       How far out the exponential fit will go. Only really affects
                      the plots and not the fit itself.
    """

    psin_avg = sliced_ts_dict["psin_avg"]
    te_avg = sliced_ts_dict["Te_avg"]
    ne_avg = sliced_ts_dict["ne_avg"]

    # Need to remove nans or the fit doesn't work.
    nan_indices = []
    for index in range(0, len(psin_avg)):
        tmp_psin = psin_avg[index]
        tmp_te   = te_avg[index]
        tmp_ne   = ne_avg[index]

        # If any of them are nan, record what index it is.
        if np.isnan(tmp_psin) or np.isnan(tmp_te) or np.isnan(tmp_ne):
            nan_indices.append(index)

    # For each nan index remove that element.
    psin_avg = np.delete(psin_avg, nan_indices)
    te_avg   = np.delete(te_avg, nan_indices)
    ne_avg   = np.delete(ne_avg, nan_indices)

    # This will return the index where avg_TS_psins first goes below 1.0 (min_psin). For ex, if
    # the array has [1.06, 1.03, 0.98, ...] it will return the index of 0.98, which is 2
    # in this case. Start exp fit from this point.
    first_sol_index = np.argmax(psin_avg<min_psin)

    # Fit needs to be shifted to the start of our SOL.
    def exp_fit(x, a, b, c):
        return a * np.exp(-b * (x-psin_avg[first_sol_index])) + c

    # Fit for Te.
    popt_te, pcov_te = curve_fit(exp_fit, psin_avg[:first_sol_index], te_avg[:first_sol_index], maxfev=5000)
    fit_in_sol_te = exp_fit(np.arange(min_psin, max_psin, 0.01), *popt_te)

    # Fit for ne. Multiply by 10^-19 so the fit uses smaller numbers (works better that way).
    fit_ne_avg = ne_avg * 10**(-19)
    guess = (1, 1e-6, 1)
    popt_ne, pcov_ne = curve_fit(exp_fit, psin_avg[:first_sol_index], fit_ne_avg[:first_sol_index], maxfev=5000, p0=guess)
    fit_in_sol_ne = exp_fit(np.arange(min_psin, max_psin, 0.01), *popt_ne) * 10 **(19)

    # Update dictionary and return.
    sliced_ts_dict["exp_psin"] = np.arange(min_psin, max_psin, 0.01)
    sliced_ts_dict["exp_te"]   = fit_in_sol_te
    sliced_ts_dict["exp_ne"]   = fit_in_sol_ne
    sliced_ts_dict["min_psin"] = min_psin
    sliced_ts_dict["max_psin"] = max_psin

    return sliced_ts_dict

def plot_exp_over_sliced(sliced_ts_dict, shot):
    """
    Plot the exponential fit over the sliced data to eyeball how good a fit
    it is.
    """

    lower_bound = sliced_ts_dict["min_psin"]
    upper_bound = sliced_ts_dict["max_psin"]

    te_upper = np.max(sliced_ts_dict["exp_te"]) * 1.2
    ne_upper = np.max(sliced_ts_dict["exp_ne"]) * 1.2

    plt.figure(1)
    plt.suptitle(r"$\mathrm{Exponential\ Fits\ of\ T_e\ and\ n_e\ Shot\ \#}$" + str(shot))
    plt.subplot(211)
    plt.plot(sliced_ts_dict["exp_psin"], sliced_ts_dict["exp_te"], "r")
    plt.plot(sliced_ts_dict["psin_avg"], sliced_ts_dict["Te_avg"], "r.")
    plt.axis([lower_bound, upper_bound, 0, te_upper])
    plt.xlabel(r"$\mathrm{\psi_n}$")
    plt.ylabel(r"$\mathrm{T_e\ (eV)}$")

    plt.subplot(212)
    plt.plot(sliced_ts_dict["exp_psin"], sliced_ts_dict["exp_ne"])
    plt.plot(sliced_ts_dict["psin_avg"], sliced_ts_dict["ne_avg"], "b.")
    plt.axis([lower_bound, upper_bound, 0, ne_upper])
    plt.xlabel(r"$\mathrm{\psi_n}$")
    plt.ylabel(r"$\mathrm{n_e\ (m^{-3})}$")

    plt.show()

def runScript(filename, shot, min_psin=1.0, max_psin=1.4, return_dict=True):
    """
    Runs the above functions in order. Normally good enough.
    """

    sliced_dict = netcdf_to_dict(filename)
    sliced_dict = fit_from_avg(sliced_dict, min_psin, max_psin)
    plot_exp_over_sliced(sliced_dict, shot)

    if return_dict:
        return sliced_dict
