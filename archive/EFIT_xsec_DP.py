import numpy as np
import pylab as plt
import scipy.interpolate as scinter

import MDSplus as mds
import DIIID.plot_vv_coils as pvc
import EFIT.load_gfile_d3d as loadg


def plot_EFIT_helper(shot, time, probe_tip, levs=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
                     tree='EFIT01', SOLlevs=[1.026, 1.049, 1.074, 1.1, 1.125], write2file=False,
                     plotIT=False, server='atlas.gat.com'):

    # gNAM = 'g'+str(shot)+'.0'+str(time)
    MDSplusCONN = mds.Connection(server)
    parmDICT = loadg.read_g_file_mds(shot, time, tree=tree, connection=MDSplusCONN)
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
    f_psiN = scinter.Rbf(Rs[Rs_trunc], Zs[Rs_trunc], parmDICT['psiRZn'][Rs_trunc],
                         function='linear')
    f_Romp = scinter.interp2d(parmDICT['psiRZn'][Rs_trunc], Zs[Rs_trunc], Rs[Rs_trunc])
    # remap each array from probe location to mag. axis midplane
    alpha_A = 0.521 / 100.0
    alpha_B = 0.0323 / 100.0
    alpha_C = 0.0168 / 100.0
    lamb = 1.27 / 100.0
    psiN_A = f_psiN(probe_tip+alpha_A, -0.18)
    psiN_B = f_psiN(probe_tip+lamb+alpha_B, -0.1546)
    psiN_C = f_psiN(probe_tip+lamb+alpha_C, -0.2054)

    R_omp_A = f_Romp(psiN_A, Z_axis)
    R_omp_B = f_Romp(psiN_B, Z_axis)
    R_omp_C = f_Romp(psiN_C, Z_axis)
    print "R_sep_OMP: {:.6}".format(f_Rs(parmDICT['ZmAxis']))
    print "A probe R_sep_OMP: {:.6}".format(R_omp_A[0])
    print "B probe R_sep_OMP: {:.6}".format(R_omp_B[0])
    print "C probe R_sep_OMP: {:.6}".format(R_omp_C[0])

    if plotIT:
        plt.figure(figsize=(7, 10), facecolor='w')
        plt.subplot(111, aspect='equal')
        plt.axis([0.9, 2.5, -1.5, 1.5])

        plt.plot([probe_tip, 2.365], [-0.18, -0.18], 'r', linewidth=12)
        plt.plot([probe_tip+0.0127, 2.365], [-0.1546, -0.1546], 'orange', linewidth=4)
        plt.plot([probe_tip+0.0127, 2.365], [-0.2054, -0.2054], 'green', linewidth=2)

        plt.plot(parmDICT['wall'][:, 0], parmDICT['wall'][:, 1], 'k', linewidth=3)
        plt.plot(pvc.RVVIN, pvc.ZVVIN, 'k', linewidth=1)
        plt.plot(pvc.RVVOUT, pvc.ZVVOUT, 'k', linewidth=1)
        plt.contour(Rs, Zs, parmDICT['psiRZn'], levs, colors='b', label='EFIT',
                    linestyles='dashed')
        plt.contour(Rs, Zs, parmDICT['psiRZn'], SOLlevs, colors='b', linestyles='solid')
        plt.contour(Rs, Zs, parmDICT['psiRZn'], [0.9996], colors='b', linewidths=2)
        plt.plot(parmDICT['lcfs'][:, 0], parmDICT['lcfs'][:, 1], 'b', linewidth=2)
        plt.plot(parmDICT['RmAxis'], parmDICT['ZmAxis'], 'b', linewidth=1, marker='+',
                 markersize=10)
        plt.plot(R_omp_A, parmDICT['ZmAxis'], marker='o', markersize=7,
                 markerfacecolor='None', mec='red', mew=2)
        plt.plot(R_omp_B, parmDICT['ZmAxis'], marker='s', markersize=7,
                 markerfacecolor='None', mec='orange', mew=2)
        plt.plot(R_omp_C, parmDICT['ZmAxis'], marker='^', markersize=7,
                 markerfacecolor='None', mec='green', mew=2)
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
