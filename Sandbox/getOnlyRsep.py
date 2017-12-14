import MDSplus as mds
import EFIT.load_gfile_d3d as loadg
import numpy as np
import scipy.interpolate as scinter



def getRsep(startTime=2500, endTime=5000, server='atlas.gat.com'):

    MDSplusCONN = mds.Connection(server)

    # Lines from Zeke's code.
    parmDICT = loadg.read_g_file_mds(shot, time,
                                     connection=MDSplusCONN,
                                     write2file=False)
    Rs, Zs = np.meshgrid(parmDICT['R'], parmDICT['Z'])
    Zes = np.copy(parmDICT['lcfs'][:, 1][13:-12])
    Res = np.copy(parmDICT['lcfs'][:, 0][13:-12])
    f_Rs = scinter.interp1d(Zes, Res, assume_sorted=False)

    # R_Sep for each z location of the three probes.
    rSep = {}
    rSep['a'] = f_Rs(-0.18)
    rSep['b'] = f_Rs(-0.1546)
    rSep['c'] = f_Rs(-0.2054)

    print "Rsep: " + str(rSep['a'])
