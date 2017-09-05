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
    """Contains information about probes. Must run r2d2, then atlas to
    fill out the variables. A list of available data after running each
    is as follows:

    letter         - Letter of probe: A, B or C
    number         - Number of probe.
    shots          - List of shots the probe was in for.
    r_probe        - Radial location of probe holder tip.
    locations_U    - The actual location along the probe in cm. 0 would be the
                     side closest to the plasma.
    w_areal_U      - Areal density of tungsten in W/cm^2. Corresponds to each location
                     or rminrsep. The lists are ordered to match up.
    w_areal_err_U  - Error of the above.
    rminrsep_U     - Distance between the average position of the sepatrix and the
                     corresponding location along the probe. The list is ordered and lines
                     up with w_areal and locations.
    rminrsep_err_U - Std. dev. of the above.
    rminrsep_omp_U - The average distance from the sepatrix after translating up
                     to the outboard midplane (omp). This is the magnetic axis of the plasma.
    rminrsep_err_U - Std. dev. of the above.

    The 'U' can be swapped out for a 'D' to get the D side of the probe data.

    """

    def __init__(self, letter, number):
        self.letter = letter
        self.number = number

    def r2d2(self, server='r2d2.gat.com'):
        self.conn           = pull.thin_connect(self.number, server=server)
        self.shots          = pull.pull_shots(self.conn, self.letter + 'D')
        self.r_probe        = pull.pull_rprobe(self.conn, self.letter)

        # Special cases of bad shots.
        # This corresponds to shot 167206.
        if self.number==2:
            self.shots = np.delete(self.shots, 10)

        # Relevant lists for each probe, upstream (U) and downstream (D).
        self.locations_U    = []
        self.w_areal_U      = []
        self.w_areal_err_U  = []
        self.locations_D    = []
        self.w_areal_D      = []
        self.w_areal_err_D  = []


        # Data from R2D2.
        for run in range(1,1000):
            # Upstream data.
            try:
                print self.letter + "U Run: " + str(run)

                loc = pull.pull_rbs_loc(self.conn, self.letter + 'U', run) / 10.0
                areal = pull.pull_rbs_areal(self.conn, self.letter + 'U', run)
                areal_err = pull.pull_rbs_areal_err(self.conn, self.letter + 'U', run)
                self.locations_U.append(loc)
                self.w_areal_U.append(areal)
                self.w_areal_err_U.append(areal_err)
            except:
                break

        for run in range(1,1000):
            # Downstream data.
            try:
                print self.letter + "D Run: " + str(run)
                loc = pull.pull_rbs_loc(self.conn, self.letter + 'D', run) / 10.0
                areal = pull.pull_rbs_areal(self.conn, self.letter + 'D', run)
                areal_err = pull.pull_rbs_areal_err(self.conn, self.letter + 'D', run)
                self.locations_D.append(loc)
                self.w_areal_D.append(areal)
                self.w_areal_err_D.append(areal_err)
            except:
                break

    def atlas(self, server='atlas.gat.com', startTime=2500, endTime=5000, step=500, efitTree='EFIT01'):

        # Create lists required to hold the data we want.
        self.rminrsep_U         = []
        self.rminrsep_err_U     = []
        self.rminrsep_omp_U     = []
        self.rminrsep_omp_err_U = []
        self.rminrsep_D         = []
        self.rminrsep_err_D     = []
        self.rminrsep_omp_D     = []
        self.rminrsep_omp_err_D = []

        # Get the dictionaries from the avg_Rsep_all function that give us all our
        # rminrsep relevant data.
        print "Analyzing " + self.letter + "U" + str(self.number) + " data..."
        avg_dict_U = get.avg_Rsep_all(self.shots, self.r_probe, self.locations_U,
            server=server, startTime=startTime, endTime=endTime, step=500, Etree=efitTree)
        print "Analyzing " + self.letter + "D" + str(self.number) + " data..."
        avg_dict_D = get.avg_Rsep_all(self.shots, self.r_probe, self.locations_D,
            server=server, startTime=startTime, endTime=endTime, step=500, Etree=efitTree)

        # Pull the data out of the dicitonaries and put into the lists in the Probe class.
        probe_name = self.letter + 'U'
        for loc in self.locations_U:
            self.rminrsep_U.append(avg_dict_U[probe_name.lower()][str(loc)])
            self.rminrsep_err_U.append(avg_dict_U[probe_name.lower() + '_err'][str(loc)])
            self.rminrsep_omp_U.append(avg_dict_U[probe_name.lower() + '_omp'][str(loc)])
            self.rminrsep_omp_err_U.append(avg_dict_U[probe_name.lower() + '_omp_err'][str(loc)])

        probe_name = self.letter + 'D'
        for loc in self.locations_D:
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
        plt.ylabel("W Areal Density (cm^-2)")
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
        plt.ylabel("W Areal Density (cm^-2)")
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

def get_multiple(aNumber=None, bNumber=None, cNumber=None, startTime=2500, endTime=5000, step=500, efitTree='EFIT01'):
    """ Allows filling out of multiple probe classes at once for probes that
        were inserted together. Returns a list of the up to three probes requested,
        each of class Probe so as to preserve each individual function in the
        Probe class (such as to_matlab)."""
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
    for p in pList:
        p.r2d2(server='localhost')

    # Get the Atlas data.
    raw_input("SSH into Atlas. Press enter when finished...")
    for p in pList:
        p.atlas(server='localhost', startTime=startTime, endTime=endTime, step=step, efitTree=efitTree)

    # Give a warning if the shots don't match up.
    # Case if only two probes are given.
    if len(pList==2):
        if (pList[0].shots.all() != pList[1].shots.all()):
            print("Error \n-------")
            print(pList[0].letter + " probe shots do not match " + pList[1].letter + " probe shots.")
            print(pList[0].letter + " shots: " + str(pList[0].shots))
            print(pList[1].letter + " shots: " + str(pList[1].shots) + "\n")

    # Case if three probes are given.
    if len(pList==3):
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

    # Return list. Could have up to three probes, but will still be in order or A, B then C.
    return pList
