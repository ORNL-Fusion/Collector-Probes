import matplotlib.pyplot as plt
import numpy as np
import EFIT.load_gfile_d3d as loadg
import MDSplus as mds

def return_avg_Z(shots, startTime=2500, endTime=5000, timeStep=500, server='localhost'):
    zArray = np.array([])

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
            print "Z Axis: " + str(Zaxis) + "\n"
            zArray = np.append(zArray, Zaxis)

            time += timeStep

    avgZ = np.mean(zArray)
    err = np.std(zArray)

    zInfo = {"avgZ":avgZ, "err":err}

    return zInfo
