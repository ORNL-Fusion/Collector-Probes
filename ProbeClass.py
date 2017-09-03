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

    def __init__(self, Anumber=None, Bnumber=None, Cnumber=None):
            self.Anumber = Anumber
            self.Bnumber = Bnumber
            self.Cnumber = Cnumber
            if Anumber is not None:
                self.Aletter = 'A'
            if Bnumber is not None:
                self.Bletter = 'B'
            if Cnumber is not None:
                self.Cletter = 'C'

    def r2d2(self, server='r2d2.gat.com'):

        # Setup each probe if used.
        if self.Anumber is not None:
            conn           = pull.thin_connect(self.Anumber, server=server)
            A_shots        = pull.pull_shots(conn, 'AD')
            self.A_r_probe = pull.pull_rprobe(conn, 'A')
            # Relevant lists for each probe, 'U' and 'D' faces.
            self.AU_locations    = []
            self.AU_w_areal      = []
            self.AU_w_areal_err  = []
            self.AD_locations    = []
            self.AD_w_areal      = []
            self.AD_w_areal_err  = []
        if self.Bnumber is not None:
            Bconn      = pull.thin_connect(self.Bnumber, server=server)
            B_shots   = pull.pull_shots(Bconn, 'BD')
            self.B_r_probe = pull.pull_rprobe(Bconn, 'B')
            # Relevant lists for each probe, 'U' and 'D' faces.
            self.BU_locations    = []
            self.BU_w_areal      = []
            self.BU_w_areal_err  = []
            self.BD_locations    = []
            self.BD_w_areal      = []
            self.BD_w_areal_err  = []
        if self.Cnumber is not None:
            self.conn      = pull.thin_connect(self.Cnumber, server=server)
            C_shots   = pull.pull_shots(conn, self.Cletter + 'D')
            self.C_r_probe = pull.pull_rprobe(conn, self.Cletter)
            ## Relevant lists for each probe, 'U' and 'D' faces.
            self.CU_locations    = []
            self.CU_w_areal      = []
            self.CU_w_areal_err  = []
            self.CD_locations    = []
            self.CD_w_areal      = []
            self.CD_w_areal_err  = []

        # Now check that the shots are consistent. Implies that the probe set is equal.
        # In most cases this is true.
        # Can default to single probe pulls.
        if self.Anumber is not None and self.Bnumber is not None and self.Cnumber is not None:
            if all(A_shots) == all(B_shots) == all(C_shots):
                shots2DICT = A_shots
            else:
                return 'Shot lists are not the same for each probe, check documentation for correct probes or get one probe at a time. Cannot continue.'
        if self.Anumber is not None and self.Bnumber is not None and self.Cnumber is None:
            if all(A_shots) == all(B_shots):
                shots2DICT = A_shots
            else:
                return 'Shot lists are not the same for each probe, check documentation for correct probes or get one probe at a time. Cannot continue.'
        elif self.Anumber is not None:
            shots2DICT = A_shots
        elif self.Bnumber is not None:
            shots2DICT = B_shots
        elif self.Cnumber is not None:
            shots2DICT = C_shots

        # this is needed b/c 167206 has no EFITs (as of 08/22/2017) -- EAU
        ind = np.argwhere(shots2DICT == 167206)
        shots2DICT = np.delete(shots2DICT, ind)
        # temporary until we get the EFIT04s in the database.
        ind = np.argwhere(shots2DICT == 167220)
        shots2DICT = np.delete(shots2DICT, ind)

        # Collect Data from R2D2 for each probe if approriate.
        if self.Anumber is not None:
            print "AU Runs"
            for run in range(1, 1000):
                # U-face data.
                try:
                    loc = pull.pull_rbs_loc(Aconn, 'AU', run) / 10.0
                    areal = pull.pull_rbs_areal(Aconn, 'AU', run)
                    areal_err = pull.pull_rbs_areal_err(Aconn, 'AU', run)
                    self.AU_locations.append(loc)
                    self.AU_w_areal.append(areal)
                    self.AU_w_areal_err.append(areal_err)
                except:
                    print "AU Stop at Run: " + str(run)
                    break
            print "AD Runs"
            for run in range(1, 1000):
                # D-face data.
                try:
                    loc = pull.pull_rbs_loc(self.conn, 'AD', run) / 10.0
                    areal = pull.pull_rbs_areal(self.conn, 'AD', run)
                    areal_err = pull.pull_rbs_areal_err(self.conn, 'AD', run)
                    self.AD_locations.append(loc)
                    self.AD_w_areal.append(areal)
                    self.AD_w_areal_err.append(areal_err)
                except:
                    print "AD Stop at Run: " + str(run)
                    break

        if self.Bnumber is not None:
            print "BU Runs"
            for run in range(1, 1000):
                try:
                    loc = pull.pull_rbs_loc(self.conn, 'BU', run) / 10.0
                    areal = pull.pull_rbs_areal(self.conn, 'BU', run)
                    areal_err = pull.pull_rbs_areal_err(self.conn, 'BU', run)
                    self.BU_locations.append(loc)
                    self.BU_w_areal.append(areal)
                    self.BU_w_areal_err.append(areal_err)
                except:
                    print "BU Stop at Run: " + str(run)
                    break
            print "BD Runs"
            for run in range(1, 1000):
                # D-face data.
                    try:
                        loc = pull.pull_rbs_loc(self.conn, self.Bletter + 'D', run) / 10.0
                        areal = pull.pull_rbs_areal(self.conn, self.Bletter + 'D', run)
                        areal_err = pull.pull_rbs_areal_err(self.conn, self.Bletter + 'D', run)
                        self.BD_locations.append(loc)
                        self.BD_w_areal.append(areal)
                        self.BD_w_areal_err.append(areal_err)
                    except:
                        print "BD Stop at Run: " + str(run)
                        break

        if self.Cnumber is not None:
            print "CU Runs"
            for run in range(1, 1000):
                try:
                    loc = pull.pull_rbs_loc(self.conn, 'CU', run) / 10.0
                    areal = pull.pull_rbs_areal(self.conn, 'CU', run)
                    areal_err = pull.pull_rbs_areal_err(self.conn, 'CU', run)
                    self.CU_locations.append(loc)
                    self.CU_w_areal.append(areal)
                    self.CU_w_areal_err.append(areal_err)
                except:
                    print "CU Stop at Run: " + str(run)
                    break
            print "CD Runs"
            for run in range(1, 1000):
                # D-face data.
                    try:
                        loc = pull.pull_rbs_loc(self.conn, 'CD', run) / 10.0
                        areal = pull.pull_rbs_areal(self.conn, 'CD', run)
                        areal_err = pull.pull_rbs_areal_err(self.conn, 'CD', run)
                        self.CD_locations.append(loc)
                        self.CD_w_areal.append(areal)
                        self.CD_w_areal_err.append(areal_err)
                    except:
                        print "CD Stop at Run: " + str(run)
                        break

        # Finally, do a check to see that r_probe is the same for all probes.
        # It should be, if not should kill script and think about what's what.
        # If it's OK, return a dictionary with all the data.
        if self.Anumber is not None:
            r_probe2DICT = self.A_r_probe
            self.r2d2DICT = {'shots': shots2DICT, 'r_probe': r_probe2DICT,
                             'AU_locations': np.array(self.AU_locations),
                             'AD_locations': np.array(self.AD_locations),
                             'AU_w_areal': np.array(self.AU_w_areal),
                             'AD_w_areal': np.array(self.AD_w_areal),
                             'AU_w_areal_err': np.array(self.AU_w_areal_err),
                             'AD_w_areal_err': np.array(self.AD_w_areal_err)}
        if self.Bnumber is not None:
            r_probe2DICT = self.B_r_probe
            self.r2d2DICT = {'shots': shots2DICT, 'r_probe': r_probe2DICT,
                             'BU_locations': np.array(self.BU_locations),
                             'BD_locations': np.array(self.BD_locations),
                             'BU_w_areal': np.array(self.BU_w_areal),
                             'BD_w_areal': np.array(self.BD_w_areal),
                             'BU_w_areal_err': np.array(self.BU_w_areal_err),
                             'BD_w_areal_err': np.array(self.BD_w_areal_err)}
        if self.Cnumber is not None:
            r_probe2DICT = self.C_r_probe
            self.r2d2DICT = {'shots': shots2DICT, 'r_probe': r_probe2DICT,
                             'CU_locations': np.array(self.CU_locations),
                             'CD_locations': np.array(self.CD_locations),
                             'CU_w_areal': np.array(self.CU_w_areal),
                             'CD_w_areal': np.array(self.CD_w_areal),
                             'CU_w_areal_err': np.array(self.CU_w_areal_err),
                             'CD_w_areal_err': np.array(self.CD_w_areal_err)}
        if self.Anumber is not None and self.Bnumber is not None:
            if self.A_r_probe == self.B_r_probe:
                r_probe2DICT = self.A_r_probe
                self.r2d2DICT = {'shots': shots2DICT, 'r_probe': r_probe2DICT,
                                 'AU_locations': np.array(self.AU_locations),
                                 'AD_locations': np.array(self.AD_locations),
                                 'AU_w_areal': np.array(self.AU_w_areal),
                                 'AD_w_areal': np.array(self.AD_w_areal),
                                 'AU_w_areal_err': np.array(self.AU_w_areal_err),
                                 'AD_w_areal_err': np.array(self.AD_w_areal_err),
                                 'BU_locations': np.array(self.BU_locations),
                                 'BD_locations': np.array(self.BD_locations),
                                 'BU_w_areal': np.array(self.BU_w_areal),
                                 'BD_w_areal': np.array(self.BD_w_areal),
                                 'BU_w_areal_err': np.array(self.BU_w_areal_err),
                                 'BD_w_areal_err': np.array(self.BD_w_areal_err)}
            else:
                print 'r_probe values not consistent between probes. Check data.'
                self.r2d2DICT ={}
        if self.Anumber is not None and self.Bnumber is not None and self.Cnumber is not None:
            if self.A_r_probe == self.B_r_probe == self.C_r_probe:
                r_probe2DICT = self.A_r_probe
                self.r2d2DICT = {'shots': shots2DICT, 'r_probe': r_probe2DICT,
                                 'AU_locations': np.array(self.AU_locations),
                                 'AD_locations': np.array(self.AD_locations),
                                 'AU_w_areal': np.array(self.AU_w_areal),
                                 'AD_w_areal': np.array(self.AD_w_areal),
                                 'AU_w_areal_err': np.array(self.AU_w_areal_err),
                                 'AD_w_areal_err': np.array(self.AD_w_areal_err),
                                 'BU_locations': np.array(self.BU_locations),
                                 'BD_locations': np.array(self.BD_locations),
                                 'BU_w_areal': np.array(self.BU_w_areal),
                                 'BD_w_areal': np.array(self.BD_w_areal),
                                 'BU_w_areal_err': np.array(self.BU_w_areal_err),
                                 'BD_w_areal_err': np.array(self.BD_w_areal_err),
                                 'CU_locations': np.array(self.CU_locations),
                                 'CD_locations': np.array(self.CD_locations),
                                 'CU_w_areal': np.array(self.CU_w_areal),
                                 'CD_w_areal': np.array(self.CD_w_areal),
                                 'CU_w_areal_err': np.array(self.CU_w_areal_err),
                                 'CD_w_areal_err': np.array(self.CD_w_areal_err)}
            else:
                print 'r_probe values not consistent between probes. Check data.'
                self.r2d2DICT ={}

        return self.r2d2DICT

    def atlas(self, server='atlas.gat.com', EFIT='EFIT01', startTime=2500, endTime=5000, step=500):
        # inputs for final dictionary.
        self.EFIT_tstart = startTime
        self.EFIT_tend = endTime
        self.EFIT_step = step

        # for A probes
        if self.Anumber is not None:
            AU_rminrsep         = []
            AU_rminrsep_err     = []
            AU_rminrsep_omp     = []
            AU_rminrsep_omp_err = []
            AD_rminrsep         = []
            AD_rminrsep_err     = []
            AD_rminrsep_omp     = []
            AD_rminrsep_omp_err = []

            print "Analyzing AU data..."
            avg_dict_U = get.avg_Rsep_all(self.r2d2DICT['shots'],
                                          self.r2d2DICT['r_probe'],
                                          self.r2d2DICT['AU_locations'],
                                          server=server,
                                          Etree=EFIT,
                                          startTime=startTime,
                                          endTime=endTime,
                                          step=step)
            print "Analyzing AD data..."
            avg_dict_D = get.avg_Rsep_all(self.r2d2DICT['shots'],
                                          self.r2d2DICT['r_probe'],
                                          self.r2d2DICT['AD_locations'],
                                          server=server,
                                          Etree=EFIT,
                                          startTime=startTime,
                                          endTime=endTime,
                                          step=step)

            probe_name = 'AU'
            for loc in self.r2d2DICT['AU_locations']:
                # print "U Loc: " + str(loc)
                AU_rminrsep.append(avg_dict_U[probe_name.lower()][str(loc)])
                AU_rminrsep_err.append(avg_dict_U[probe_name.lower() + '_err'][str(loc)])
                AU_rminrsep_omp.append(avg_dict_U[probe_name.lower() + '_omp'][str(loc)])
                AU_rminrsep_omp_err.append(avg_dict_U[probe_name.lower() + '_omp_err'][str(loc)])

            probe_name = 'AD'
            # for key in avg_dict[probe_name.lower()].keys():
            #    print key
            for loc in self.r2d2DICT['AD_locations']:
                # print "D Loc: " + str(loc)
                AD_rminrsep.append(avg_dict_D[probe_name.lower()][str(loc)])
                AD_rminrsep_err.append(avg_dict_D[probe_name.lower() + '_err'][str(loc)])
                AD_rminrsep_omp.append(avg_dict_D[probe_name.lower() + '_omp'][str(loc)])
                AD_rminrsep_omp_err.append(avg_dict_D[probe_name.lower() + '_omp_err'][str(loc)])

        # for B probes
        if self.Bnumber is not None:
            BU_rminrsep         = []
            BU_rminrsep_err     = []
            BU_rminrsep_omp     = []
            BU_rminrsep_omp_err = []
            BD_rminrsep         = []
            BD_rminrsep_err     = []
            BD_rminrsep_omp     = []
            BD_rminrsep_omp_err = []

            print "Analyzing BU data..."
            avg_dict_U = get.avg_Rsep_all(self.r2d2DICT['shots'],
                                          self.r2d2DICT['r_probe'],
                                          self.r2d2DICT['BU_locations'],
                                          server=server,
                                          Etree=EFIT,
                                          startTime=startTime,
                                          endTime=endTime,
                                          step=step)
            print "Analyzing BD data..."
            avg_dict_D = get.avg_Rsep_all(self.r2d2DICT['shots'],
                                          self.r2d2DICT['r_probe'],
                                          self.r2d2DICT['BD_locations'],
                                          server=server,
                                          Etree=EFIT,
                                          startTime=startTime,
                                          endTime=endTime,
                                          step=step)

            probe_name = 'BU'
            for loc in self.r2d2DICT['BU_locations']:
                # print "U Loc: " + str(loc)
                BU_rminrsep.append(avg_dict_U[probe_name.lower()][str(loc)])
                BU_rminrsep_err.append(avg_dict_U[probe_name.lower() + '_err'][str(loc)])
                BU_rminrsep_omp.append(avg_dict_U[probe_name.lower() + '_omp'][str(loc)])
                BU_rminrsep_omp_err.append(avg_dict_U[probe_name.lower() + '_omp_err'][str(loc)])

            probe_name = 'BD'
            # for key in avg_dict[probe_name.lower()].keys():
            #    print key
            for loc in self.r2d2DICT['BD_locations']:
                # print "D Loc: " + str(loc)
                BD_rminrsep.append(avg_dict_D[probe_name.lower()][str(loc)])
                BD_rminrsep_err.append(avg_dict_D[probe_name.lower() + '_err'][str(loc)])
                BD_rminrsep_omp.append(avg_dict_D[probe_name.lower() + '_omp'][str(loc)])
                BD_rminrsep_omp_err.append(avg_dict_D[probe_name.lower() + '_omp_err'][str(loc)])

        # for C probes
        if self.Cnumber is not None:
            CU_rminrsep         = []
            CU_rminrsep_err     = []
            CU_rminrsep_omp     = []
            CU_rminrsep_omp_err = []
            CD_rminrsep         = []
            CD_rminrsep_err     = []
            CD_rminrsep_omp     = []
            CD_rminrsep_omp_err = []

            print "Analyzing CU data..."
            avg_dict_U = get.avg_Rsep_all(self.r2d2DICT['shots'],
                                          self.r2d2DICT['r_probe'],
                                          self.r2d2DICT['CU_locations'],
                                          server=server,
                                          Etree=EFIT,
                                          startTime=startTime,
                                          endTime=endTime,
                                          step=step)
            print "Analyzing CD data..."
            avg_dict_D = get.avg_Rsep_all(self.r2d2DICT['shots'],
                                          self.r2d2DICT['r_probe'],
                                          self.r2d2DICT['CD_locations'],
                                          server=server,
                                          Etree=EFIT,
                                          startTime=startTime,
                                          endTime=endTime,
                                          step=step)

            probe_name = 'CU'
            for loc in self.r2d2DICT['CU_locations']:
                # print "U Loc: " + str(loc)
                CU_rminrsep.append(avg_dict_U[probe_name.lower()][str(loc)])
                CU_rminrsep_err.append(avg_dict_U[probe_name.lower() + '_err'][str(loc)])
                CU_rminrsep_omp.append(avg_dict_U[probe_name.lower() + '_omp'][str(loc)])
                CU_rminrsep_omp_err.append(avg_dict_U[probe_name.lower() + '_omp_err'][str(loc)])

            probe_name = 'CD'
            # for key in avg_dict[probe_name.lower()].keys():
            #    print key
            for loc in self.r2d2DICT['CD_locations']:
                # print "D Loc: " + str(loc)
                CD_rminrsep.append(avg_dict_D[probe_name.lower()][str(loc)])
                CD_rminrsep_err.append(avg_dict_D[probe_name.lower() + '_err'][str(loc)])
                CD_rminrsep_omp.append(avg_dict_D[probe_name.lower() + '_omp'][str(loc)])
                CD_rminrsep_omp_err.append(avg_dict_D[probe_name.lower() + '_omp_err'][str(loc)])

        if self.Anumber is not None:
            self.atlasDICT = {'EFIT tree': EFIT, "EFIT start time": startTime,
                              'EFIT end time': endTime, 'EFIT time step': step,
                              'AU_rminrsep': np.array(AU_rminrsep),
                              'AUrminrsep_err': np.array(AU_rminrsep_err),
                              'AU_rminrsep_omp': np.array(AU_rminrsep_omp),
                              'AU_rminrsep_omp_err': np.array(AU_rminrsep_omp_err),
                              'AD_rminrsep': np.array(AD_rminrsep),
                              'AD_rminrsep_err': np.array(AD_rminrsep_err),
                              'AD_rminrsep_omp': np.array(AD_rminrsep_omp),
                              'AD_rminrsep_omp_err': np.array(AD_rminrsep_omp_err)}
        if self.Bnumber is not None:
            self.atlasDICT = {'EFIT tree': EFIT, "EFIT start time": startTime,
                              'EFIT end time': endTime, 'EFIT time step': step,
                              'BU_rminrsep': np.array(BU_rminrsep),
                              'BUrminrsep_err': np.array(BU_rminrsep_err),
                              'BU_rminrsep_omp': np.array(BU_rminrsep_omp),
                              'BU_rminrsep_omp_err': np.array(BU_rminrsep_omp_err),
                              'BD_rminrsep': np.array(BD_rminrsep),
                              'BD_rminrsep_err': np.array(BD_rminrsep_err),
                              'BD_rminrsep_omp': np.array(BD_rminrsep_omp),
                              'BD_rminrsep_omp_err': np.array(BD_rminrsep_omp_err)}
        if self.Cnumber is not None:
            self.atlasDICT = {'EFIT tree': EFIT, "EFIT start time": startTime,
                              'EFIT end time': endTime, 'EFIT time step': step,
                              'CU_rminrsep': np.array(CU_rminrsep),
                              'CUrminrsep_err': np.array(CU_rminrsep_err),
                              'CU_rminrsep_omp': np.array(CU_rminrsep_omp),
                              'CU_rminrsep_omp_err': np.array(CU_rminrsep_omp_err),
                              'CD_rminrsep': np.array(CD_rminrsep),
                              'CD_rminrsep_err': np.array(CD_rminrsep_err),
                              'CD_rminrsep_omp': np.array(CD_rminrsep_omp),
                              'CD_rminrsep_omp_err': np.array(CD_rminrsep_omp_err)}
        if self.Anumber is not None and self.Bnumber is not None:
            self.atlasDICT = {'EFIT tree': EFIT, "EFIT start time": startTime,
                              'EFIT end time': endTime, 'EFIT time step': step,
                              'AU_rminrsep': np.array(AU_rminrsep),
                              'AUrminrsep_err': np.array(AU_rminrsep_err),
                              'AU_rminrsep_omp': np.array(AU_rminrsep_omp),
                              'AU_rminrsep_omp_err': np.array(AU_rminrsep_omp_err),
                              'AD_rminrsep': np.array(AD_rminrsep),
                              'AD_rminrsep_err': np.array(AD_rminrsep_err),
                              'AD_rminrsep_omp': np.array(AD_rminrsep_omp),
                              'AD_rminrsep_omp_err': np.array(AD_rminrsep_omp_err),
                              'BU_rminrsep': np.array(BU_rminrsep),
                              'BUrminrsep_err': np.array(BU_rminrsep_err),
                              'BU_rminrsep_omp': np.array(BU_rminrsep_omp),
                              'BU_rminrsep_omp_err': np.array(BU_rminrsep_omp_err),
                              'BD_rminrsep': np.array(BD_rminrsep),
                              'BD_rminrsep_err': np.array(BD_rminrsep_err),
                              'BD_rminrsep_omp': np.array(BD_rminrsep_omp),
                              'BD_rminrsep_omp_err': np.array(BD_rminrsep_omp_err)}
        if self.Anumber is not None and self.Bnumber is not None and self.Cnumber is not None:
            self.atlasDICT = {'EFIT tree': EFIT, "EFIT start time": startTime,
                              'EFIT end time': endTime, 'EFIT time step': step,
                              'AU_rminrsep': np.array(AU_rminrsep),
                              'AUrminrsep_err': np.array(AU_rminrsep_err),
                              'AU_rminrsep_omp': np.array(AU_rminrsep_omp),
                              'AU_rminrsep_omp_err': np.array(AU_rminrsep_omp_err),
                              'AD_rminrsep': np.array(AD_rminrsep),
                              'AD_rminrsep_err': np.array(AD_rminrsep_err),
                              'AD_rminrsep_omp': np.array(AD_rminrsep_omp),
                              'AD_rminrsep_omp_err': np.array(AD_rminrsep_omp_err),
                              'BU_rminrsep': np.array(BU_rminrsep),
                              'BUrminrsep_err': np.array(BU_rminrsep_err),
                              'BU_rminrsep_omp': np.array(BU_rminrsep_omp),
                              'BU_rminrsep_omp_err': np.array(BU_rminrsep_omp_err),
                              'BD_rminrsep': np.array(BD_rminrsep),
                              'BD_rminrsep_err': np.array(BD_rminrsep_err),
                              'BD_rminrsep_omp': np.array(BD_rminrsep_omp),
                              'BD_rminrsep_omp_err': np.array(BD_rminrsep_omp_err),
                              'CU_rminrsep': np.array(CU_rminrsep),
                              'CUrminrsep_err': np.array(CU_rminrsep_err),
                              'CU_rminrsep_omp': np.array(CU_rminrsep_omp),
                              'CU_rminrsep_omp_err': np.array(CU_rminrsep_omp_err),
                              'CD_rminrsep': np.array(CD_rminrsep),
                              'CD_rminrsep_err': np.array(CD_rminrsep_err),
                              'CD_rminrsep_omp': np.array(CD_rminrsep_omp),
                              'CD_rminrsep_omp_err': np.array(CD_rminrsep_omp_err)}

        return self.atlasDICT

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
        plt.title('W Areal Density for ' + self.letter + 'U/' + self.letter + 'D Probes from ' + str(self.EFIT_tstart) + ' to ' + str(self.EFIT_tend))
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
        plt.title('W Areal Density for ' + self.letter + 'U/' + self.letter + 'D Probes at OMP from ' + str(self.EFIT_tstart) + ' to ' + str(self.EFIT_tend))
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

    def to_hickle(self):
        import Misc.hickle.hickle as hkl

        if self.Anumber is not None:
            suffix = 'A'+str(self.Anumber)
        if self.Bnumber is not None:
            suffix = 'B'+str(self.Bnumber)
        if self.Cnumber is not None:
            suffix = 'C'+str(self.Cnumber)
        if self.Anumber is not None:
            suffix = 'A'+str(self.Anumber)+'_B'+str(self.Bnumber)
        if self.Anumber is not None and self.Bnumber is not None and self.Cnumber is not None:
            suffix = 'A'+str(self.Anumber)+'_B'+str(self.Bnumber)+'_C'+str(self.Cnumber)

        fnam = 'CPdata_mdsplus_'+suffix+'.h5'
        print fnam
        hkl.dump(self.r2d2DICT, fnam, 'w')
        hkl.dump(self.atlasDICT, fnam, 'w')
