# Description:
# Pulls data from the dp_probes tree on r2d2 and puts it into dictionaries
# in the class. Then uses this data to access EFIT on atlas and get R - Rsep and
# R - Rsep_omp. These can then be plotted using the plot functions or the data
# can be used in whatever way seen fit.

import pull_data_dp as pull
import get_Rsep as get

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import warnings
#from __future__ import print_function


class Probe():
    """
    Contains information about probes. Must run from collector probe MDS+ repository first
    (r2d2.gat.com), then from the EFIT MDS+ repo (atlas.gat.com) to fill out the variables.

    Contains information about probes. Must run r2d2, then atlas to
    fill out the variables from corresponding MDS+ repositories. A list of available
    data after running each is put into dictionaries inside the class object.:

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
    EFIT tree       - Which EFIT was used.
    EFIT start time - Time of shot where the range of the averages is started.
    EFIT end time   - Time of shot where the range of averages ends.
    EFIT time step  - Time increment between start and end times.

    The 'U' can be swapped out for a 'D' to get the D side of the probe data.
    """

    def __init__(self, letter, number):
        self.letter = letter
        self.number = number

    def r2d2(self, server='r2d2.gat.com'):
        """This functions pulls all the relevant data stored on the MDS+
           repository on the R2D2 server."""

        conn    = pull.thin_connect(self.number, server=server)
        shots   = pull.pull_shots(conn,  self.letter + 'D')
        r_probe = pull.pull_rprobe(conn, self.letter)

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
        print "\n"
        for run in range(1, 1000):
            # U-face data.
            try:
                #print "\033[F \033[F"
                loc = pull.pull_rbs_loc(conn, self.letter + 'U', run) / 10.0
                areal = pull.pull_rbs_areal(conn, self.letter + 'U', run)
                areal_err = pull.pull_rbs_areal_err(conn, self.letter + 'U', run)
                locations_U.append(loc)
                w_areal_U.append(areal)
                w_areal_err_U.append(areal_err)
                print "\033[F \033[F"
                print self.letter + "U Run: " + str(run)
            except:
                break

        print "\n"
        for run in range(1, 1000):
            # D-face data.
            try:
                #print "\033[F \033[F"
                loc = pull.pull_rbs_loc(conn, self.letter + 'D', run) / 10.0
                areal = pull.pull_rbs_areal(conn, self.letter + 'D', run)
                areal_err = pull.pull_rbs_areal_err(conn, self.letter + 'D', run)
                locations_D.append(loc)
                w_areal_D.append(areal)
                w_areal_err_D.append(areal_err)
                print "\033[F \033[F"
                print self.letter + "D Run: " + str(run)
            except:
                break

        # Put all the data into a dictionary inside the class.
        self.r2d2DICT = {'shots': shots, 'r_probe': r_probe,
                         'locations_U': np.array(locations_U),
                         'locations_D': np.array(locations_D),
                         'w_areal_U': np.array(w_areal_U),
                         'w_areal_D': np.array(w_areal_D),
                         'w_areal_err_U': np.array(w_areal_err_U),
                         'w_areal_err_D': np.array(w_areal_err_D)}

    def atlas(self, server='atlas.gat.com', EFIT='EFIT01', startTime=2500, endTime=5000, step=500,
              probe_tip_corr=1.5):
        """
        This function access the EFIT data on atlas to get the R - Rsep
        and R_omp - Rsep_omp values. Averages between the start and end times
        are returned. Run this after r2d2.
        Input varibles/keywords (default values if any):
        server ('atlas.gat.com')  -- the MDSplus connection information, can be localhost
                                     if tunneling.
        EFIT ('EFIT01')           -- name of the EFIT MDSplus tree used in calculations.
        startTime (2500)          -- time in msec for first EFIT.
        endTime (5000)            -- time in msec for last EFIT.
        step (500)                -- del_time between EFIT equilibrium.
        probe_tip_corr (1.5)      -- correction to PLC insertion value in
                                     cm (from Dmitry's calibration).
        """

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

        # Get the Rsep data from EFIT and perform necessary operations in the
        # get_Rsep file.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            print "Analyzing " + self.letter + "U" + str(self.number) + " data..."
            avg_dict_U = get.avg_Rsep_all(self.r2d2DICT['shots'],
                                          self.r2d2DICT['r_probe']+probe_tip_corr,
                                          self.r2d2DICT['locations_U'],
                                          server=server,
                                          Etree=EFIT,
                                          startTime=startTime,
                                          endTime=endTime,
                                          step=step)
            print "Analyzing " + self.letter + "D" + str(self.number) + " data..."
            avg_dict_D = get.avg_Rsep_all(self.r2d2DICT['shots'],
                                          self.r2d2DICT['r_probe']+probe_tip_corr,
                                          self.r2d2DICT['locations_D'],
                                          server=server,
                                          Etree=EFIT,
                                          startTime=startTime,
                                          endTime=endTime,
                                          step=step)

        # Pull the data from the returned dictionary from get_Rsep and put it
        # into lists for each the U and D probes.
        probe_name = self.letter + 'U'
        for loc in self.r2d2DICT['locations_U']:
            rminrsep_U.append(avg_dict_U[probe_name.lower()][str(loc)])
            rminrsep_err_U.append(avg_dict_U[probe_name.lower() + '_err'][str(loc)])
            rminrsep_omp_U.append(avg_dict_U[probe_name.lower() + '_omp'][str(loc)])
            rminrsep_omp_err_U.append(avg_dict_U[probe_name.lower() + '_omp_err'][str(loc)])

        probe_name = self.letter + 'D'
        for loc in self.r2d2DICT['locations_D']:
            # print "D Loc: " + str(loc)
            rminrsep_D.append(avg_dict_D[probe_name.lower()][str(loc)])
            rminrsep_err_D.append(avg_dict_D[probe_name.lower() + '_err'][str(loc)])
            rminrsep_omp_D.append(avg_dict_D[probe_name.lower() + '_omp'][str(loc)])
            rminrsep_omp_err_D.append(avg_dict_D[probe_name.lower() + '_omp_err'][str(loc)])

        # Create atlas dictionary and store it in the class.
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

    # Plot the R - Rsep data. The limit flag is to match up with the omp graph.
    # For whatever reason the last ~6 points are garbage-like. Something with
    # the interpolation maybe.
    def plot_norm(self, newFIG=True):
        if newFIG:
            plt.figure()
        plt.errorbar(x=self.atlasDICT['rminrsep_U'],
                     y=self.r2d2DICT['w_areal_U'],
                     xerr=self.atlasDICT['rminrsep_err_U'],
                     yerr=self.r2d2DICT['w_areal_err_U'],
                     label=self.letter + 'U' + str(self.number))
        plt.errorbar(x=self.atlasDICT['rminrsep_D'],
                     y=self.r2d2DICT['w_areal_D'],
                     xerr=self.atlasDICT['rminrsep_err_D'],
                     yerr=self.r2d2DICT['w_areal_err_D'],
                     label=self.letter + 'D' + str(self.number))
        plt.legend(loc='upper right')
        plt.xlabel("R - Rsep (cm)")
        plt.ylabel('W Areal Density (x10^16 cm^-2)')
        plt.title('W Areal Density for ' + self.letter + 'U/' + self.letter + 'D Probes from ' + str(self.atlasDICT['EFIT start time']) + ' to ' + str(self.atlasDICT['EFIT end time']))
        plt.show()

    # Plot the R - Rsep omp data.
    def plot_omp(self, newFIG=True):
        if newFIG:
            plt.figure()
        plt.errorbar(x=self.atlasDICT['rminrsep_omp_U'],
                     y=self.r2d2DICT['w_areal_U'],
                     xerr=self.atlasDICT['rminrsep_omp_err_U'],
                     yerr=self.r2d2DICT['w_areal_err_U'],
                     label=self.letter + 'U' + str(self.number) + ' omp', fmt='o')
        plt.errorbar(x=self.atlasDICT['rminrsep_omp_D'],
                     y=self.r2d2DICT['w_areal_D'],
                     xerr=self.atlasDICT['rminrsep_omp_err_D'],
                     yerr=self.r2d2DICT['w_areal_err_D'],
                     label=self.letter + 'D' + str(self.number) + ' omp', fmt='o')
        plt.legend(loc='upper right')
        plt.xlabel("R - Rsep_omp (cm)")
        plt.ylabel("W Areal Density (x10^16 cm^-2)")
        plt.title('W Areal Density for ' + self.letter + 'U/' + self.letter + 'D Probes at OMP from ' + str(self.atlasDICT['EFIT start time']) + ' to ' + str(self.atlasDICT['EFIT end time']))
        plt.show()

    # Output to basic matlab file for curve fitting or whatever.
    def to_matlab(self):
        tmp_dict = {}
        arr_loc_U   = np.array(self.r2d2DICT['locations_U'])
        arr_sep_U   = np.array(self.atlasDICT['rminrsep_U'])
        arr_sep_err_U   = np.array(self.atlasDICT['rminrsep_err_U'])
        arr_omp_U   = np.array(self.atlasDICT['rminrsep_omp_U'])
        arr_omp_err_U   = np.array(self.atlasDICT['rminrsep_omp_err_U'])
        arr_areal_U = np.array(self.r2d2DICT['w_areal_U'])
        arr_areal_err_U = np.array(self.r2d2DICT['w_areal_err_U'])
        tmp_dict['locations_U']       = arr_loc_U
        tmp_dict['rminrsep_U']        = arr_sep_U
        tmp_dict['rminrsep_err_U']        = arr_sep_err_U
        tmp_dict['rminrsep_omp_U']    = arr_omp_U
        tmp_dict['rminrsep_omp_err_U']    = arr_omp_err_U
        tmp_dict['w_areal_density_U'] = arr_areal_U
        tmp_dict['w_areal_density_err_U'] = arr_areal_err_U

        arr_loc_D   = np.array(self.r2d2DICT['locations_D'])
        arr_sep_D   = np.array(self.atlasDICT['rminrsep_D'])
        arr_sep_err_D   = np.array(self.atlasDICT['rminrsep_err_D'])
        arr_omp_D   = np.array(self.atlasDICT['rminrsep_omp_D'])
        arr_omp_err_D   = np.array(self.atlasDICT['rminrsep_omp_err_D'])
        arr_areal_D = self.r2d2DICT['w_areal_D']
        arr_areal_err_D = self.r2d2DICT['w_areal_err_D']
        tmp_dict['locations_D']       = arr_loc_D
        tmp_dict['rminrsep_D']        = arr_sep_D
        tmp_dict['rminrsep_err_D']        = arr_sep_err_D
        tmp_dict['rminrsep_omp_D']    = arr_omp_D
        tmp_dict['rminrsep_omp_err_D']    = arr_omp_err_D
        tmp_dict['w_areal_density_D'] = arr_areal_D
        tmp_dict['w_areal_density_err_D'] = arr_areal_err_D

        filename = self.letter.upper() + str(self.number) + 'data.mat'
        sio.savemat(filename, tmp_dict)


def get_multiple(aNumber=None, bNumber=None, cNumber=None, MDStunnel=False, startTime=2500,
                 endTime=5000, step=500, efitTree='EFIT01', toHDF5=False):
    """
    Allows filling out of multiple probe classes at once for probes that
    were inserted together. Returns a list of the up to three probes requested,
    each of class Probe and thus with its own corresponding data in it. try
    using 'vars(ProbeObjectHere)' to get a dictionary of the data returned.
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
        if (pList[0].r2d2DICT['shots'].all() != pList[1].r2d2DICT['shots'].all()):
            print("Error \n-------")
            print(pList[0].letter + " probe shots do not match " + pList[1].letter +
                  " probe shots.")
            print(pList[0].letter + " shots: " + str(pList[0].r2d2DICT['shots']))
            print(pList[1].letter + " shots: " + str(pList[1].r2d2DICT['shots']) + "\n")

    # Case if three probes are given.
    if len(pList) == 3:
        if (pList[0].r2d2DICT['shots'].all() != pList[1].r2d2DICT['shots'].all()):
            print("Error \n-------")
            print(pList[0].letter + " probe shots do not match " + pList[1].letter +
                  " probe shots.")
            print(pList[0].letter + " shots: " + str(pList[0].r2d2DICT['shots']))
            print(pList[1].letter + " shots: " + str(pList[1].r2d2DICT['shots']) + "\n")
        if (pList[0].r2d2DICT['shots'].all() != pList[2].r2d2DICT['shots'].all()):
            print("Error \n-------")
            print(pList[0].letter + " probe shots do not match " + pList[2].letter +
                  " probe shots.")
            print(pList[0].letter + " shots: " + str(pList[0].r2d2DICT['shots']))
            print(pList[2].letter + " shots: " + str(pList[2].shots) + "\n")
        if (pList[1].r2d2DICT['shots'].all() != pList[2].r2d2DICT['shots'].all()):
            print("Error \n-------")
            print(pList[1].letter + " probe shots do not match " + pList[2].letter +
                  " probe shots.")
            print(pList[1].letter + " shots: " + str(pList[1].r2d2DICT['shots']))
            print(pList[2].letter + " shots: " + str(pList[2].r2d2DICT['shots']) + "\n")

    # Return list. Could have up to three probes, but will still be in order or A, B then C.
    return pList


# Output to HDF5 file.
def dump2HDF5(pList):
    """
    Warning: make sure you have the hickle package installed via:
             pip install hickle
    """

    import Misc.hickle.hickle as hkl
    suffix = ''
    for i, j in enumerate(pList):
        suffix += str(j.letter) + str(j.number)

    fnam = 'CPdata_mdsplus_'+suffix+'.h5'
    print "HDF5 dump: to "+fnam

    hkl.dump(pList[0].r2d2DICT, fnam, 'w', path='/'+str(pList[0].letter) + str(pList[0].number) + "_r2d2")
    hkl.dump(pList[0].atlasDICT, fnam, 'r+', path='/'+str(pList[0].letter) + str(pList[0].number) + "_atlas")

    for i, j in enumerate(pList):
        if i != 0:
            hkl.dump(pList[i].r2d2DICT, fnam, 'r+', path='/'+str(pList[i].letter) + str(pList[i].number) + "_r2d2")
            hkl.dump(pList[i].atlasDICT, fnam, 'r+', path='/'+str(pList[i].letter) + str(pList[i].number) + "_atlas")
