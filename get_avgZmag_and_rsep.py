import matplotlib.pyplot as plt
import numpy as np
import EFIT.load_gfile_d3d as loadg
import MDSplus as mds
import scipy.interpolate as scinter
import math

def return_avg_Z(shots, startTime=2500, endTime=5000, timeStep=500, server='localhost'):
    print "Z Location of A Probe: -0.18 m \nZ Location of B Probe: -0.1546 m \nZ Location of C Probe: -0.2054 m"


    zArray = np.array([])
    rSepArrayA = np.array([])
    rSepArrayB = np.array([])
    rSepArrayC = np.array([])

    MDSplusConn = mds.Connection(server)

    for shot in shots:
        print "Shot: " + str(shot)
        time = startTime
        while time <= endTime:

            print "Time: " + str(time)
            parmDICT = loadg.read_g_file_mds(shot, time,
                                             connection=MDSplusConn,
                                             write2file=False,
                                             tree='EFIT01')
            Zaxis = parmDICT['ZmAxis']
            Zes = np.copy(parmDICT['lcfs'][:, 1][13:-12])
            Res = np.copy(parmDICT['lcfs'][:, 0][13:-12])
            f_Rs = scinter.interp1d(Zes, Res, assume_sorted=False)

            # R_Sep for each z location of the three probes.
            rSep = {}
            rSep['a'] = f_Rs(-0.18)
            rSep['b'] = f_Rs(-0.1546)
            rSep['c'] = f_Rs(-0.2054)

            print "Z Axis: " + str(Zaxis) + "\n"
            zArray = np.append(zArray, Zaxis)
            rSepArrayA = np.append(rSepArrayA, rSep['a'])
            rSepArrayB = np.append(rSepArrayB, rSep['b'])
            rSepArrayC = np.append(rSepArrayC, rSep['c'])

            time += timeStep

    avgZ = np.mean(zArray)
    std_dev = np.std(zArray)
    num_of_samples = len(shots) * int((endTime-startTime)/timeStep)
    err = std_dev / math.sqrt(num_of_samples)

    avgRsepA = np.mean(rSepArrayA)
    avgRsepB = np.mean(rSepArrayB)
    avgRsepC = np.mean(rSepArrayC)
    avgRsepA_err = np.std(rSepArrayA) / math.sqrt(num_of_samples)
    avgRsepB_err = np.std(rSepArrayB) / math.sqrt(num_of_samples)
    avgRsepC_err = np.std(rSepArrayC) / math.sqrt(num_of_samples)

    returned_info = {"Average Z Magnetic Axis":avgZ, "Z Error":err,
                     "Average Rsep at A Probe":avgRsepA, "A Error":avgRsepA_err,
                     "Average Rsep at B Probe":avgRsepB, "B Error":avgRsepB_err,
                     "Average Rsep at C Probe":avgRsepC, "C Error":avgRsepC_err}

    return returned_info
