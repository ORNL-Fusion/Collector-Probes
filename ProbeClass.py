# File:   ProbeClass.py
# Author: Shawn Zamperini
# Email:  zamp@utk.edu
# Date:   8/23/17
#
# Description:
# Pulls data from the dp_probes tree on r2d2 and puts it into lists
# in the class. Then uses this data to access EFIT on atlas and get R - Rsep and
# R - Rsep_omp. These can then be plotted using the plot functions or the data
# can be used in whatever way seen fit.

import pull_data_dp as pull
import get_Rsep as get

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio


class Probe():
    """
    Contains information about probes. Must run from collector probe MDS+ repository first
    (r2d2.gat.com), then from the EFIT MDS+ repo (atlas.gat.com) to fill out the variables.
    """

    def __init__(self, letter, number):
        self.letter = letter
        self.number = number

    def r2d2(self, server='r2d2.gat.com'):
        conn    = pull.thin_connect(self.number, server=server)
        shots   = pull.pull_shots(conn,  self.letter + 'D')
        r_probe = pull.pull_rprobe(conn, 'A')

        # Relevant lists for each probe, 'U' and 'D' faces.
        locations_U    = []
        w_areal_U      = []
        w_areal_err_U  = []
        locations_D    = []
        w_areal_D      = []
        w_areal_err_D  = []

        # this is needed b/c 167206 has no EFITs (as of 08/22/2017) -- EAU
        ind = np.argwhere(shots == 167206)
        shots = np.delete(shots, ind)
        # temporary until we get the EFIT04s in the database.
        ind = np.argwhere(shots == 167220)
        shots = np.delete(shots, ind)

        # Data from R2D2.
        for run in range(1, 1000):
            # U-face data.
            try:
                print self.letter + "U Run: " + str(run)

                loc = pull.pull_rbs_loc(conn, self.letter + 'U', run) / 10.0
                areal = pull.pull_rbs_areal(conn, self.letter + 'U', run)
                areal_err = pull.pull_rbs_areal_err(conn, self.letter + 'U', run)
                locations_U.append(loc)
                w_areal_U.append(areal)
                w_areal_err_U.append(areal_err)
            except:
                break

        for run in range(1, 1000):
            # D-face data.
            try:
                print self.letter + "D Run: " + str(run)
                loc = pull.pull_rbs_loc(conn, self.letter + 'D', run) / 10.0
                areal = pull.pull_rbs_areal(conn, self.letter + 'D', run)
                areal_err = pull.pull_rbs_areal_err(conn, self.letter + 'D', run)
                locations_D.append(loc)
                w_areal_D.append(areal)
                w_areal_err_D.append(areal_err)
            except:
                break

        self.r2d2DICT = {'shots': shots, 'r_probe': r_probe,
                         'locations_U': np.array(locations_U),
                         'locations_D': np.array(locations_D),
                         'w_areal_U': np.array(w_areal_U),
                         'w_areal_D': np.array(w_areal_D),
                         'w_areal_err_U': np.array(w_areal_err_U),
                         'w_areal_err_D': np.array(w_areal_err_D)}

    def atlas(self, server='atlas.gat.com', EFIT='EFIT01', startTime=2500, endTime=5000, step=500):
        # inputs for final dictionary.
        self.EFIT_tstart = startTime
        self.EFIT_tend = endTime
        self.EFIT_step = step

        rminrsep_U         = []
        rminrsep_err_U     = []
        rminrsep_omp_U     = []
        rminrsep_omp_err_U = []
        rminrsep_D         = []
        rminrsep_err_D     = []
        rminrsep_omp_D     = []
        rminrsep_omp_err_D = []

        print "Analyzing " + self.letter + "U" + str(self.number) + " data..."
        avg_dict_U = get.avg_Rsep_all(self.r2d2DICT['shots'],
                                      self.r2d2DICT['r_probe'],
                                      self.r2d2DICT['locations_U'],
                                      server=server,
                                      Etree=EFIT,
                                      startTime=startTime,
                                      endTime=endTime,
                                      step=step)
        print "Analyzing " + self.letter + "D" + str(self.number) + " data..."
        avg_dict_D = get.avg_Rsep_all(self.r2d2DICT['shots'],
                                      self.r2d2DICT['r_probe'],
                                      self.r2d2DICT['locations_D'],
                                      server=server,
                                      Etree=EFIT,
                                      startTime=startTime,
                                      endTime=endTime,
                                      step=step)

        probe_name = self.letter + 'U'
        for loc in self.r2d2DICT['locations_U']:
            # print "U Loc: " + str(loc)
            rminrsep_U.append(avg_dict_U[probe_name.lower()][str(loc)])
            rminrsep_err_U.append(avg_dict_U[probe_name.lower() + '_err'][str(loc)])
            rminrsep_omp_U.append(avg_dict_U[probe_name.lower() + '_omp'][str(loc)])
            rminrsep_omp_err_U.append(avg_dict_U[probe_name.lower() + '_omp_err'][str(loc)])

        probe_name = self.letter + 'D'
        # for key in avg_dict[probe_name.lower()].keys():
        #    print key
        for loc in self.r2d2DICT['locations_D']:
            # print "D Loc: " + str(loc)
            rminrsep_D.append(avg_dict_D[probe_name.lower()][str(loc)])
            rminrsep_err_D.append(avg_dict_D[probe_name.lower() + '_err'][str(loc)])
            rminrsep_omp_D.append(avg_dict_D[probe_name.lower() + '_omp'][str(loc)])
            rminrsep_omp_err_D.append(avg_dict_D[probe_name.lower() + '_omp_err'][str(loc)])

        self.atlasDICT = {'EFIT tree': EFIT, "EFIT start time": startTime,
                          'EFIT end time': endTime, 'EFIT time step': step,
                          'rminrsep_U': np.array(rminrsep_U),
                          'rminrsep_err_U': np.array(rminrsep_err_U),
                          'rminrsep_omp_U': np.array(rminrsep_omp_U),
                          'rminrsep_omp_err_U': np.array(rminrsep_omp_err_U),
                          'rminrsep_D': np.array(rminrsep_D),
                          'rminrsep_err_D': np.array(rminrsep_err_D),
                          'rminrsep_omp_D': np.array(rminrsep_omp_D),
                          'rminrsep_omp_err_D': np.array(rminrsep_omp_err_D)}

    def plot_norm(self, limit=6, newFIG=True):
        if newFIG:
            plt.figure()
        plt.errorbar(x=self.atlasDICT['rminrsep_U'][limit:], y=self.r2d2DICT['w_areal_U'][limit:],
                     xerr=self.atlasDICT['rminrsep_err_U'][limit:], yerr=self.r2d2DICT['w_areal_err_U'][limit:],
                     label=self.letter + 'U' + str(self.number))
        plt.errorbar(x=self.atlasDICT['rminrsep_D'][limit:], y=self.r2d2DICT['w_areal_D'][limit:],
                     xerr=self.atlasDICT['rminrsep_err_D'][limit:], yerr=self.r2d2DICT['w_areal_err_D'][limit:],
                     label=self.letter + 'D' + str(self.number))
        plt.legend(loc='upper right')
        plt.xlabel("R - Rsep (cm)")
        plt.ylabel('W Areal Density (x10^16 cm^-2)')
        plt.title('W Areal Density for ' + self.letter + 'U/' + self.letter + 'D Probes from ' + str(self.atlasDICT['EFIT start time']) + ' to ' + str(self.atlasDICT['EFIT end time']))
        plt.show()

    def plot_omp(self, limit=6, newFIG=True):
        if newFIG:
            plt.figure()
        plt.errorbar(x=self.atlasDICT['rminrsep_omp_U'][limit:],
                     y=self.r2d2DICT['w_areal_U'][limit:],
                     xerr=self.atlasDICT['rminrsep_omp_err_U'][limit:],
                     yerr=self.r2d2DICT['w_areal_err_U'][limit:],
                     label=self.letter + 'U' + str(self.number) + ' omp', fmt='o')
        plt.errorbar(x=self.atlasDICT['rminrsep_omp_D'][limit:],
                     y=self.r2d2DICT['w_areal_D'][limit:],
                     xerr=self.atlasDICT['rminrsep_omp_err_D'][limit:],
                     yerr=self.r2d2DICT['w_areal_err_D'][limit:],
                     label=self.letter + 'D' + str(self.number) + ' omp', fmt='o')
        plt.legend(loc='upper right')
        plt.xlabel("R - Rsep_omp (cm)")
        plt.ylabel("W Areal Density (x10^16 cm^-2)")
        plt.title('W Areal Density for ' + self.letter + 'U/' + self.letter + 'D Probes at OMP from ' + str(self.atlasDICT['EFIT start time']) + ' to ' + str(self.atlasDICT['EFIT end time']))
        plt.show()

    def to_matlab(self):
        tmp_dict = {}
        arr_loc_U   = np.array(self.r2d2DICT['locations_U'])
        arr_sep_U   = np.array(self.rminrsep_U)
        arr_omp_U   = np.array(self.rminrsep_omp_U)
        arr_areal_U = np.array(self.r2d2DICT['w_areal_U'])
        tmp_dict['locations_U']       = arr_loc_U
        tmp_dict['rminrsep_U']        = arr_sep_U
        tmp_dict['rminrsep_omp_U']    = arr_omp_U
        tmp_dict['w_areal_density_U'] = arr_areal_U

        arr_loc_D   = np.array(self.r2d2DICT['locations_D'])
        arr_sep_D   = np.array(self.rminrsep_D)
        arr_omp_D   = np.array(self.rminrsep_omp_D)
        arr_areal_D = self.r2d2DICT['w_areal_D']
        tmp_dict['locations_D']       = arr_loc_D
        tmp_dict['rminrsep_D']        = arr_sep_D
        tmp_dict['rminrsep_omp_D']    = arr_omp_D
        tmp_dict['w_areal_density_D'] = arr_areal_D

        filename = self.letter.upper() + str(self.number) + 'data.mat'
        sio.savemat(filename, tmp_dict)


def get_multiple(aNumber=None, bNumber=None, cNumber=None, MDStunnel=False, startTime=2500,
                 endTime=5000, step=500, efitTree='EFIT01', toHDF5=False):
    """
    Allows filling out of multiple probe classes at once for probes that
    were inserted together. Returns a list of the up to three probes requested,
    each of class Probe so as to preserve each individual function in the
    Probe class (such as to_matlab).
    """
    # Create probes and put into list.
    pList = []
    if aNumber:
        aProbe = Probe('A', int(aNumber))
        pList.append(aProbe)
    if bNumber:
        bProbe = Probe('B', int(bNumber))
        pList.append(bProbe)
    if cNumber:
        cProbe = Probe('C', int(cNumber))
        pList.append(cProbe)

    # Get the R2D2 data.
    raw_input("SSH into R2D2. Press enter when finished...")
    if MDStunnel:
        server = 'localhost'
    else:
        server = 'r2d2.gat.com'
    for p in pList:
        p.r2d2(server=server)

    # Get the Atlas data.
    raw_input("SSH into Atlas. Press enter when finished...")
    for p in pList:
        if MDStunnel:
            server = 'localhost'
        else:
            server = 'atlas.gat.com'
        p.atlas(server=server, startTime=startTime, endTime=endTime, step=step,
                EFIT=efitTree)

    # Give a warning if the shots don't match up.
    # Case if only two probes are given.
    if len(pList) == 2:
        if (pList[0].shots.all() != pList[1].shots.all()):
            print("Error \n-------")
            print(pList[0].letter + " probe shots do not match " + pList[1].letter + " probe shots.")
            print(pList[0].letter + " shots: " + str(pList[0].shots))
            print(pList[1].letter + " shots: " + str(pList[1].shots) + "\n")

    # Case if three probes are given.
    if len(pList) == 3:
        if (pList[0].shots.all() != pList[1].shots.all()):
            print("Error \n-------")
            print(pList[0].letter + " probe shots do not match " + pList[1].letter + " probe shots.")
            print(pList[0].letter + " shots: " + str(pList[0].shots))
            print(pList[1].letter + " shots: " + str(pList[1].shots) + "\n")
        if (pList[0].shots.all() != pList[2].shots.all()):
            print("Error \n-------")
            print(pList[0].letter + " probe shots do not match " + pList[2].letter + " probe shots.")
            print(pList[0].letter + " shots: " + str(pList[0].shots))
            print(pList[2].letter + " shots: " + str(pList[2].shots) + "\n")
        if (pList[1].shots.all() != pList[2].shots.all()):
            print("Error \n-------")
            print(pList[1].letter + " probe shots do not match " + pList[2].letter + " probe shots.")
            print(pList[1].letter + " shots: " + str(pList[1].shots))
            print(pList[2].letter + " shots: " + str(pList[2].shots) + "\n")

    if toHDF5:
        """
        Warning: make sure you have the hickle package installed via:
                 pip install hickle
        """
        import Misc.hickle.hickle as hkl


    # Return list. Could have up to three probes, but will still be in order or A, B then C.
    return pList

#
#     if self.Anumber is not None:
#         suffix = 'A'+str(self.Anumber)
#     if self.Bnumber is not None:
#         suffix = 'B'+str(self.Bnumber)
#     if self.Cnumber is not None:
#         suffix = 'C'+str(self.Cnumber)
#     if self.Anumber is not None:
#         suffix = 'A'+str(self.Anumber)+'_B'+str(self.Bnumber)
#     if self.Anumber is not None and self.Bnumber is not None and self.Cnumber is not None:
#         suffix = 'A'+str(self.Anumber)+'_B'+str(self.Bnumber)+'_C'+str(self.Cnumber)
#
#     fnam = 'CPdata_mdsplus_'+suffix+'.h5'
#     print fnam
#     hkl.dump(self.r2d2DICT, fnam, 'w', path='/r2d2')
#     hkl.dump(self.atlasDICT, fnam, 'r+', path='/atlas')
