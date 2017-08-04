import numpy as np
import pylab as plt
import scipy.interpolate as scinter

import DIIID.plot_vv_coils as pvc
import EFIT.load_gfile_d3d as loadg


def plot_EFIT_helper(shot, time, probe_tip, levs=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
                     tree='EFIT01', SOLlevs=[1.026, 1.049, 1.074, 1.1, 1.125], write2file=False,
                     plotIT=False, server='atlas.gat.com'):

    # gNAM = 'g'+str(shot)+'.0'+str(time)

    parmDICT = loadg.read_g_file_mds(shot, time, tree=tree, Server=server)
    Rs, Zs = np.meshgrid(parmDICT['R'], parmDICT['Z'])
    Z_axis = parmDICT['ZmAxis']
    R_axis = parmDICT['RmAxis']
    # Have to resrict the Zes and Res range becuase they are "double valued". Trying to keep
    # everything on LFS.
    Zes = np.copy(parmDICT['lcfs'][:, 1][13:-12])
    Res = np.copy(parmDICT['lcfs'][:, 0][13:-12])
    f_Rs = scinter.interp1d(Zes, Res, assume_sorted=False)

    # translate up to outbard midplane
    Rs_trunc = Rs > R_axis
    f_psiN = scinter.interp2d(Rs[Rs_trunc], Zs[Rs_trunc], parmDICT['psiRZn'][Rs_trunc])
    f_Romp = scinter.interp2d(parmDICT['psiRZn'][Rs_trunc], Zs[Rs_trunc], Rs[Rs_trunc])

    # remap each array from probe location to mag. axis midplane
    psiN_A = f_psiN([2.26,2.27], -0.18)
    #psiN_B = f_psiN(Rs, -0.1546)
    #psiN_C = f_psiN(Rs, -0.2054)

    R_omp_AU = f_Romp(psiN_AU, Z_axis)

    print f_Rs(Z_axis)
    print psiN_A
    print R_omp_A
    print "A probe R_sep: {:.6}".format(f_Rs(-0.18))
    print "B probe R_sep: {:.6}".format(f_Rs(-0.1546))
    print "C probe R_sep: {:.6}".format(f_Rs(-0.2054))

    if plotIT:
        plt.figure(figsize=(7, 10), facecolor='w')
        plt.subplot(111, aspect='equal')
        plt.axis([0.9, 2.5, -1.5, 1.5])

        plt.plot([probe_tip, 2.365], [-0.18, -0.18], 'r', linewidth=7)
        plt.plot([probe_tip+0.01, 2.365], [-0.1546, -0.1546], 'r', linewidth=3)
        plt.plot([probe_tip+0.01, 2.365], [-0.2054, -0.2054], 'r', linewidth=2)

        plt.plot(parmDICT['wall'][:, 0], parmDICT['wall'][:, 1], 'k', linewidth=3)
        plt.plot(pvc.RVVIN, pvc.ZVVIN, 'k', linewidth=1)
        plt.plot(pvc.RVVOUT, pvc.ZVVOUT, 'k', linewidth=1)
        plt.contour(Rs, Zs, parmDICT['psiRZn'], levs, colors='b', label='EFIT', linestyles='dashed')
        plt.contour(Rs, Zs, parmDICT['psiRZn'], SOLlevs, colors='b', linestyles='solid')
        plt.contour(Rs, Zs, parmDICT['psiRZn'], [0.9996], colors='b', linewidths=2)
        plt.plot(parmDICT['lcfs'][:, 0], parmDICT['lcfs'][:, 1], 'b', linewidth=2)
        plt.plot(parmDICT['RmAxis'], parmDICT['ZmAxis'], 'b', linewidth=1, marker='+', markersize=10)
        plt.show()

    dict = {'R_sepA': f_Rs(-0.18), 'R_sepB': f_Rs(-0.1546), 'R_sepC': f_Rs(-0.2054)}
    return dict

def find_closest(A, target):
    # A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx
