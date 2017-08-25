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
        self.conn           = pull.thin_connect(self.number, server=server)
        self.shots          = pull.pull_shots(self.conn, self.letter + 'D')
        self.r_probe        = pull.pull_rprobe(self.conn, self.letter)

        # Relevant lists for each probe, upstream (U) and downstream (D).
        self.locations_U    = []
        self.w_areal_U      = []
        self.w_areal_err_U  = []
        self.locations_D    = []
        self.w_areal_D      = []
        self.w_areal_err_D  = []

        # this is needed b/c 167206 has no EFITs (as of 08/22/2017) -- EAU
        ind = np.argwhere(self.shots == 167206)
        self.shots = np.delete(self.shots, ind)

        # Data from R2D2.
        for run in range(1, 1000):
            # Upstream data.
            try:
                print "AU Run: " + str(run)

                loc = pull.pull_rbs_loc(self.conn, self.letter + 'U', run) / 10.0
                areal = pull.pull_rbs_areal(self.conn, self.letter + 'U', run)
                areal_err = pull.pull_rbs_areal_err(self.conn, self.letter + 'U', run)
                self.locations_U.append(loc)
                self.w_areal_U.append(areal)
                self.w_areal_err_U.append(areal_err)
            except:
                break

        for run in range(1, 1000):
            # Downstream data.
            try:
                print "AD Run: " + str(run)
                loc = pull.pull_rbs_loc(self.conn, self.letter + 'D', run) / 10.0
                areal = pull.pull_rbs_areal(self.conn, self.letter + 'D', run)
                areal_err = pull.pull_rbs_areal_err(self.conn, self.letter + 'D', run)
                self.locations_D.append(loc)
                self.w_areal_D.append(areal)
                self.w_areal_err_D.append(areal_err)
            except:
                break

    def atlas(self, server='atlas.gat.com', EFIT='EFIT01', startTime=2500, endTime=5000, step=500):

        self.rminrsep_U         = []
        self.rminrsep_err_U     = []
        self.rminrsep_omp_U     = []
        self.rminrsep_omp_err_U = []
        self.rminrsep_D         = []
        self.rminrsep_err_D     = []
        self.rminrsep_omp_D     = []
        self.rminrsep_omp_err_D = []

        print "Analyzing " + self.letter + "U" + str(self.number) + " data..."
        avg_dict_U = get.avg_Rsep_all(self.shots, self.r_probe, self.locations_U,
                                      server=server, Etree=EFIT, startTime=startTime,
                                      endTime=endTime, step=500)
        print "Analyzing " + self.letter + "D" + str(self.number) + " data..."
        avg_dict_D = get.avg_Rsep_all(self.shots, self.r_probe, self.locations_D,
                                      server=server, Etree=EFIT, startTime=startTime,
                                      endTime=endTime, step=500)

        probe_name = self.letter + 'U'
        for loc in self.locations_U:
            # print "U Loc: " + str(loc)
            self.rminrsep_U.append(avg_dict_U[probe_name.lower()][str(loc)])
            self.rminrsep_err_U.append(avg_dict_U[probe_name.lower() + '_err'][str(loc)])
            self.rminrsep_omp_U.append(avg_dict_U[probe_name.lower() + '_omp'][str(loc)])
            self.rminrsep_omp_err_U.append(avg_dict_U[probe_name.lower() + '_omp_err'][str(loc)])

        probe_name = self.letter + 'D'
        # for key in avg_dict[probe_name.lower()].keys():
        #    print key

        for loc in self.locations_D:
            # print "D Loc: " + str(loc)
            self.rminrsep_D.append(avg_dict_D[probe_name.lower()][str(loc)])
            self.rminrsep_err_D.append(avg_dict_D[probe_name.lower() + '_err'][str(loc)])
            self.rminrsep_omp_D.append(avg_dict_D[probe_name.lower() + '_omp'][str(loc)])
            self.rminrsep_omp_err_D.append(avg_dict_D[probe_name.lower() + '_omp_err'][str(loc)])

    def plot_norm(self, limit=6):

        plt.errorbar(x=self.rminrsep_U[limit:], y=self.w_areal_U[limit:],
                     xerr=self.rminrsep_err_U[limit:], yerr=self.w_areal_err_U[limit:],
                     label=self.letter + 'U' + str(self.number))
        plt.errorbar(x=self.rminrsep_D[limit:], y=self.w_areal_D[limit:],
                     xerr=self.rminrsep_err_D[limit:], yerr=self.w_areal_err_D[limit:],
                     label=self.letter + 'D' + str(self.number))
        plt.legend(loc='upper right')
        plt.xlabel("R - Rsep (cm)")
        plt.ylabel('W Areal Density (x10^16 cm^-2)')
        plt.title('W Areal Density for ' + self.letter + 'U/' + self.letter + 'D Probes')

        plt.show()

    def plot_omp(self, limit=6):
        plt.errorbar(x=self.rminrsep_omp_U[limit:], y=self.w_areal_U[limit:],
                     xerr=self.rminrsep_omp_err_U[limit:], yerr=self.w_areal_err_U[limit:],
                     label=self.letter + 'U' + str(self.number) + ' omp', fmt='o')
        plt.errorbar(x=self.rminrsep_omp_D[limit:], y=self.w_areal_D[limit:],
                     xerr=self.rminrsep_omp_err_D[limit:], yerr=self.w_areal_err_D[limit:],
                     label=self.letter + 'D' + str(self.number) + ' omp', fmt='o')
        plt.legend(loc='upper right')
        plt.xlabel("R - Rsep_omp (cm)")
        plt.ylabel("W Areal Density (x10^16 cm^-2)")
        plt.title('W Areal Density for ' + self.letter + 'U/' + self.letter + 'D Probes at OMP')

        plt.show()

    def to_matlab(self):
        tmp_dict={}
        arr_loc_U   = np.array(self.locations_U)
        arr_sep_U   = np.array(self.rminrsep_U)
        arr_omp_U   = np.array(self.rminrsep_omp_U)
        arr_areal_U = np.array(self.w_areal_U)
        tmp_dict['locations_U']       = arr_loc_U
        tmp_dict['rminrsep_U']        = arr_sep_U
        tmp_dict['rminrsep_omp_U']    = arr_omp_U
        tmp_dict['w_areal_density_U'] = arr_areal_U

        arr_loc_D   = np.array(self.locations_D)
        arr_sep_D   = np.array(self.rminrsep_D)
        arr_omp_D   = np.array(self.rminrsep_omp_D)
        arr_areal_D = np.array(self.w_areal_D)
        tmp_dict['locations_D']       = arr_loc_D
        tmp_dict['rminrsep_D']        = arr_sep_D
        tmp_dict['rminrsep_omp_D']    = arr_omp_D
        tmp_dict['w_areal_density_D'] = arr_areal_D

        filename = self.letter.upper() + str(self.number) + 'data.mat'
        sio.savemat(filename, tmp_dict)
