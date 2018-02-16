import sys
sys.path.append("\home\shawn\DIII-D\ORNL-Fusion\Collector-Probes")

import ts_to_farsol as ts
import ProbeClass as Probe
from scipy import interpolate
import numpy as np



# Mass of Dueterium in eV s^2 m^-2
mass_deut = 2.01 * 931.49 * 10**6 / ((3*10**8)**2.0)
diff_coeff = 1.0
# Probe widths in m
aSize = 3.0 / 100.0
bSize = 1.0 / 100.0
cSize = 0.5 / 100.0

def calc_pert_lengths(ts_filename, ts_shot, pList=None, aNumber=None, bNumber=None, cNumber=None, MDStunnel=True, startTime=2500, endTime=5000):
    data_dict = {}

    if pList==None:
        # Get the list of Probe objects.
        pList = Probe.get_multiple(aNumber, bNumber, cNumber, MDStunnel, startTime, endTime)

    # Make sure TS shot that we have the netcDF file for is in the probe.
    if ts_shot not in pList[0].r2d2DICT["shots"]:
        print "Warning: Probes were not in for TS shot."

    # Get the TS data, where it is an exp. fit out to where the probes are.
    ts_extrap_data = ts.runScript(ts_filename, ts_shot, max_psin=1.5)
    ts_psin = ts_extrap_data["exp_psin"]
    ts_te   = ts_extrap_data["exp_te"]

    # Create a function of these values.
    f_ts_te = interpolate.interp1d(ts_psin, ts_te)

    # Won't distinguish perturbation lengths between sides, so just use AD psins.
    aProbe = pList[0]
    ad_psin = aProbe.atlasDICT["psiN_D"]
    data_dict["ad_psin"] = ad_psin
    try:
        bProbe = pList[1]
        bd_psin = bProbe.atlasDICT["psiN_D"]
        data_dict["bd_psin"] = bd_psin
    except:
        print "No B Probe."
    try:
        cProbe = pList[2]
        cd_psin = cProbe.atlasDICT["psiN_D"]
        data_dict["cd_psin"] = cd_psin
    except:
        print "No C Probe"

    # For each probe get the sound speed calculated from the TS Te data.
    sound_speeds = np.array([])
    for index in range(0, len(ad_psin)):
        tmp_psin = ad_psin[index]
        tmp_te = f_ts_te(tmp_psin)
        tmp_cs = (tmp_te*2 / mass_deut)**0.5
        sound_speeds = np.append(sound_speeds, tmp_cs)

    # Calculate the perturbation lengths.
    lengths_A = sound_speeds * aSize**2 / (8.0 * diff_coeff)
    lengths_B = sound_speeds * bSize**2 / (8.0 * diff_coeff)
    lengths_C = sound_speeds * cSize**2 / (8.0 * diff_coeff)


    #data_dict["sound_speeds"] = sound_speeds    
    data_dict["pert_lengths_A"] = lengths_A
    data_dict["pert_lengths_B"] = lengths_B
    data_dict["pert_lengths_C"] = lengths_C
    data_dict["pList"] = pList

    return data_dict


    #return pert_lengths
