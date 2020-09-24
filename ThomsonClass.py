# Imports for python2 implementation.
from __future__ import print_function, unicode_literals
from __future__ import absolute_import, division

import sys, os
import numpy             as np
import pandas            as pd
import MDSplus           as mds
import matplotlib        as mpl
import matplotlib.pyplot as plt
from matplotlib        import ticker, cm, colors
from scipy.interpolate import Rbf, interp1d
from gadata            import gadata


class ThomsonClass:
    """
    Examples of how to use ThomsonClass object.

    Example #1:
      ts = ThomsonClass(176343, 'divertor')
      ts.load_ts()
      ts.map_to_efit(times=np.linspace(2500, 5000, 10), ref_time=2500)
      ts.heatmap(detach_front_te=5, offset=0.1)

    Example #2:
      ts = ThomsonClass(176343, 'core')
      ts.load_ts()
      ts.map_to_efit(np.linspace(2500, 5000, 10))
      ts.avg_omp  # Dictionary of average omp values.

    """

    def __init__(self, shot, system):
        self.shot   = shot
        self.system = system
        self.conn = None
        self.ts_dict = {}

    def __repr__(self):
        return "ThomsonClass Object\n" \
               + "  Shot:   " + str(self.shot) + "\n" \
               + "  System: " + str(self.system)

    def load_ts(self, verbal=True, times=None, tunnel=True, filter=False,
                fs='FS04', avg_thresh=2, method='simple', window_len=11):
        """
        Function to get all the data from the Thomson Scattering "BLESSED" tree
        on atlas. The data most people probably care about is the temperature (eV)
        and density (m-3) data in ts_dict. In the 'temp' and 'density' entries,
        'X' is the time, and 'Y' is a 2D array where each row is the data for a
        single chord at each of those times. Mapping to psin coordinates and such
        are done in a different function.

        tunnel: If using locally, set to True. This mean you need an ssh tunnel connected
                  to local host. The command:

                    ssh -Y -p 2039 -L 8000:atlas.gat.com:8000 username@cybele.gat.com

                  should work in a separate terminal. Set to False if on DIII-D network.
        times: Provide times that a polynomial fit will be applied to, and then
               averaged over. This is my attempt at smoothing the TS data when
               it's noisy (maybe to help with ELMs).
        filter: Run the filter function to do a simple filter of the data.
        """

        # Create thin connection to MDSplus on atlas. Tunnel if connection locally.
        if tunnel:
            conn = mds.Connection("localhost")
        else:
            conn = mds.Connection("atlas.gat.com")

        # Store the connection object in the class for later use.
        self.conn = conn

        # Open the tree of the shot we want.
        tree = conn.openTree("d3d", self.shot)

        # Will need to probably update this as I go, since I'm not sure what
        # the earliest BLESSED shot was (or which is REVISION01 for that matter).
        if self.shot <= 94741:
            base = "\\D3D::TOP.ELECTRONS.TS.REVISIONS.REVISION00"
        else:
            base = "\\D3D::TOP.ELECTRONS.TS.BLESSED"

        # Specify which system we want.
        print("Thomson system: " + self.system)
        if self.system is "core":
            base = base + ".CORE"
        elif self.system is "divertor":
            base = base + ".DIVERTOR"
        elif self.system is "tangential":
            base = base + ".TANGENTIAL"

        # List containing all the names of the nodes under BLESSED.
        nodes = ["CALCMASK", "CDPOLYBOX", "CDPULSE", "CHANNEL", "CHI_MAX",
                 "CHI_MIN", "DCDATA", "DENSITY", "DENSITY_E", "DETOPT",
                 "DCPEDESTAL", "DETOPT", "FITDATA", "FITTHRESHOLD", "FRACCHI",
                 "INIT_NE", "INIT_TE", "ITMICRO", "LFORDER", "LPROF", "MAXFITS",
                 "PHI", "PLDATA", "PLERROR", "PLPEDESTAL", "PLPEDVAR", "R",
                 "REDCHISQ", "SUPOPT", "TEMP", "TEMP_E", "THETA", "TIME", "Z"]

        if not verbal:
            print("Loading Thomson data...")

        # Get the data for each node, and put it in a dictionary of X and Y values.
        # For 1D data, like "Z", the X will be the channel number, and the Y will
        # be the Z coordinate.
        for node in nodes:
            try:
                # Make the path to the node.
                path = base + "." + node
                if verbal:
                    print("Getting data from node: " + path)

                # Get the data(Y) and dimension data(X) from the node.
                data = conn.get(path).data()
                data_dim = conn.get("DIM_OF(" + path + ")").data()

                # Put into a dictionary then put it into the master dictionary.
                data_dict = {"X":data_dim, "Y":data}
                self.ts_dict[node.lower()] = data_dict

                # If the node is TEMP or DENSITY, there are nodes beneath it. I haven't
                # seen data in these nodes ever (except R and Z, but they're redundant
                # and the same as the R and Z nodes as above), but check anyways.
                if node in ["TEMP", "DENSITY"]:
                    for subnode in ["PHI", "PSI01", "PSI02", "R", "RHO01", "RHO02", "Z"]:
                        path = base + "." + node + "." + subnode
                        if verbal:
                            print("Getting data from node: " + path)
                        try:
                            data = conn.get(path).data()
                            data_dim = conn.get("dim_of(" + path + ")").data()
                            data_dict = {"Time":data_dim, subnode.lower():data}
                            self.ts_dict[node.lower() + "." + subnode.lower()] = data_dict
                        except (mds.MdsIpException, mds.TreeNODATA):
                            if verbal:
                                print("  Node has no data.")

            # This error is returned if the node is empty. Catch it.
            #except (mds.MdsIpException, mds.TreeNODATA, mds.MdsException):
            #except (mds.MdsException):
            #    if verbal:
            #        print("  Node has no data.")

            # For compatablity with some of the older shots.
            #except mds.TreeNNF:
            #    if verbal:
            #        print("  Node not found.")

            except:
                if verbal:
                    print("  Node not found/no data.")

        # Pull these into DataFrames, a logical and easy way to represent the
        # data. Initially the rows are each a TS chord, and the columns are at
        # each time.
        self.temp_df = pd.DataFrame(columns=self.ts_dict['temp']['X'],    data=self.ts_dict['temp']['Y'])
        self.dens_df = pd.DataFrame(columns=self.ts_dict['density']['X'], data=self.ts_dict['density']['Y'])
        self.temp_df.index.name   = 'Chord'
        self.temp_df.columns.name = 'Time (ms)'
        self.dens_df.index.name   = 'Chord'
        self.dens_df.columns.name = 'Time (ms)'

        # Transpose the data so each row is at a specific time, and the columns
        # are the chords.
        self.temp_df = self.temp_df.transpose()
        self.dens_df = self.dens_df.transpose()

        # Filter the data from ELMs and just replace with the filtered dataframes.
        if filter:
            print("Filtering data...")
            self.temp_df_unfiltered = self.temp_df
            self.dens_df_unfiltered = self.dens_df
            self.filter_elms(fs=fs, avg_thresh=avg_thresh, method=method, window_len=window_len)
            self.temp_df = self.temp_df_filt
            self.dens_df = self.dens_df_filt

        # Do a polynomial fit to the data. Fifth-order.
        if times is not None:
            self.temp_df_poly = pd.DataFrame()
            self.dens_df_poly = pd.DataFrame()
            xp = np.linspace(times.min(), times.max(), 1000)
            for chord in self.temp_df:

                # Limit time between desired times.
                x = self.temp_df.index.values
                idxs = np.where(np.logical_and(x >= times.min(), x <= times.max()))[0]
                x = x[idxs]
                y_te = self.temp_df[chord].values[idxs]
                y_ne = self.dens_df[chord].values[idxs]
                z_te = np.polyfit(x, y_te, 5) # Returns five exponents.
                z_ne = np.polyfit(x, y_ne, 5)
                p_te = np.poly1d(z_te) # Creates the fit with the exponents.
                p_ne = np.poly1d(z_ne)
                yp_te = p_te(xp) # Use the fit to create 1,000 points.
                yp_ne = p_ne(xp)

                # Put into the Dataframe.
                self.temp_df_poly[chord] = yp_te
                self.dens_df_poly[chord] = yp_ne

            # Give the index the times.
            self.temp_df_poly.index = xp
            self.dens_df_poly.index = xp

            # Can we just swap temp_df out with temp_df_poly?
            #self.temp_df = self.temp_df_poly
            #self.dens_df = self.dens_df_poly


    def load_gfile_mds(self, shot, time, tree="EFIT04", exact=False,
                       connection=None, tunnel=True, verbal=True):
        """
        This is scavenged from the load_gfile_d3d script on the EFIT repository,
        except updated to run on python3.

        shot:       Shot to get gfile for.
        time:       Time of the shot to load gfile for, in ms.
        tree:       One of the EFIT trees to get the data from.
        exact:      If True will raise error if time does not exactly match any gfile
                    times. False will grab the closest time.
        connection: An MDSplus connection to atlas.
        tunnel:     Set to True if accessing outside DIII-D network.

        returns:    The requested gfile as a dictionary.
        """

        # Connect to server, open tree and go to g-file
        if connection is None:
            if tunnel is True:
                connection = mds.Connection("localhost")
            else:
                connection = mds.Connection('atlas.gat.com')
        connection.openTree(tree, shot)

        base = 'RESULTS:GEQDSK:'

        # get time slice
        if verbal:
            print("Loading gfile:")
            print("  Shot: " + str(shot))
            print("  Tree: " + tree)
            print("  Time: " + str(time))
        signal = 'GTIME'
        k = np.argmin(np.abs(connection.get(base + signal).data() - time))
        time0 = int(connection.get(base + signal).data()[k])

        if (time != time0):
            if exact:
                raise RuntimeError(tree + ' does not exactly contain time %.2f'\
                                   %time + '  ->  Abort')
            else:
                if verbal:
                    print('Warning: Closest time is ' + str(time0) +'.')
                    #print('Fetching time slice ' + str(time0))
                time = time0

        # store data in dictionary
        g = {'shot': shot, 'time': time}

        # get header line
        try:
            header = connection.get(base + 'ECASE').data()[k]
        except:
            print("  No header line.")

        # get all signals, use same names as in read_g_file
        translate = {'MW': 'NR', 'MH': 'NZ', 'XDIM': 'Xdim', 'ZDIM': 'Zdim', 'RZERO': 'R0',
                     'RMAXIS': 'RmAxis', 'ZMAXIS': 'ZmAxis', 'SSIMAG': 'psiAxis', 'SSIBRY': 'psiSep',
                     'BCENTR': 'Bt0', 'CPASMA': 'Ip', 'FPOL': 'Fpol', 'PRES': 'Pres',
                     'FFPRIM': 'FFprime', 'PPRIME': 'Pprime', 'PSIRZ': 'psiRZ', 'QPSI': 'qpsi',
                     'NBBBS': 'Nlcfs', 'LIMITR': 'Nwall'}
        for signal in translate:
            try:
                g[translate[signal]] = connection.get(base + signal).data()[k]
            except:
                print("  Node not found: " + base + signal)

        g['R1'] = connection.get(base + 'RGRID').data()[0]
        g['Zmid'] = 0.0

        RLIM = connection.get(base + 'LIM').data()[:, 0]
        ZLIM = connection.get(base + 'LIM').data()[:, 1]
        g['wall'] = np.vstack((RLIM, ZLIM)).T

        RBBBS = connection.get(base + 'RBBBS').data()[k][:int(g['Nlcfs'])]
        ZBBBS = connection.get(base + 'ZBBBS').data()[k][:int(g['Nlcfs'])]
        g['lcfs'] = np.vstack((RBBBS, ZBBBS)).T

        KVTOR = 0
        RVTOR = 1.7
        NMASS = 0
        RHOVN = connection.get(base + 'RHOVN').data()[k]

        # convert floats to integers
        for item in ['NR', 'NZ', 'Nlcfs', 'Nwall']:
            g[item] = int(g[item])

        # convert single (float32) to double (float64) and round
        for item in ['Xdim', 'Zdim', 'R0', 'R1', 'RmAxis', 'ZmAxis', 'psiAxis',
                     'psiSep', 'Bt0', 'Ip']:
            g[item] = np.round(np.float64(g[item]), 7)

        # convert single arrays (float32) to double arrays (float64)
        for item in ['Fpol', 'Pres', 'FFprime', 'Pprime', 'psiRZ', 'qpsi',
                     'lcfs', 'wall']:
            g[item] = np.array(g[item], dtype=np.float64)

        # Construct (R,Z) grid for psiRZ
        g['dR'] = g['Xdim']/(g['NR'] - 1)
        g['R'] = g['R1'] + np.arange(g['NR'])*g['dR']

        g['dZ'] = g['Zdim']/(g['NZ'] - 1)
        NZ2 = int(np.floor(0.5*g['NZ']))
        g['Z'] = g['Zmid'] + np.arange(-NZ2, NZ2+1)*g['dZ']

        # normalize psiRZ
        g['psiRZn'] = (g['psiRZ'] - g['psiAxis']) / (g['psiSep'] - g['psiAxis'])

        return g

    def map_to_efit(self, times=None, ref_time=None, average_ts=5, debug=False,
                    tree='EFIT04', choose_interp_region=False, trunc_div=False):
        """
        This function uses the gfiles of each time in times to find the location
        of each TS chord relative to the X-point in polar coordinates, (d, theta).
        By doing this for a swept strike point and mapping each (d, theta) back
        to a reference frame, 2D profiles of Te and ne can be obtained.

        BUG: Mapping the core to the X-point during a sweep doesn't work. The
          X-point can move while the core plasma doesn't, thus a moving X-point
          could make it seem like the core sweeps a range too when it doesn't.
          A fix here could be mapping the core to the plasma center, while the
          divertor TS stays with the X-point.

        Note: Currently only written for LSN. Needs to be updated for USN.

        times:      A list of times to average over and map to. Can be a single float
                      or list. The gfile for each time will be loaded.
        ref_time:   The time for the reference frame, which the 2D plots will be
                      plotted over. For a left to right sweep, choose the time
                      where the strike point is furthest to the left.
        average_ts: A parameter to help smooth the data a bit. Instead of just
                      taking the TS data at the time that matches the current
                      gile, it will take the previous 5 and next 5 TS times and
                      average them.
        """


        # Just use the first time as the ref_time (this may be run without
        # needing it so catch it).
        if ref_time is None:
            ref_time = times[0]

        times = np.array(times)

        # Store times for later use in plotting function.
        self.times = np.array(times)
        self.ref_time = ref_time
        self.temp_df_omp = pd.DataFrame()
        self.dens_df_omp = pd.DataFrame()

        # Load gfile(s).
        self.ref_df = pd.DataFrame()
        self.ref_df.index.name =   'Chord'
        self.ref_df.columns.name = 'Time (R, Z)'
        count = 1

        # Current implementation makes the first time the reference frame. So make
        # times have the ref_time in front.
        ref_idx = np.where(times == ref_time)[0][0]
        times = np.append(ref_time, np.append(times[:ref_idx], times[ref_idx+1:]))
        ref_flag = True

        # Defaults for fitting the LCFS interpolation functions.
        start = 13; end = 71
        for time in times:
            print('\nLoading gfile (' + str(count) + '/' + str(len(times)) + ')...')
            try:
                gfile = self.load_gfile_mds(shot=self.shot, time=time, connection=self.conn, verbal=True, tree=tree)

                # Store for plotting function.
                if ref_flag:
                    self.ref_gfile = gfile

                # Z's and R's of the separatrix, in m.
                Zes = np.copy(gfile['lcfs'][:, 1])
                Res = np.copy(gfile['lcfs'][:, 0])

                # The location of the X-point is here where the lowest Z is (LSN).
                xpoint_idx = np.where(Zes == Zes.min())[0][0]
                Rx = Res[xpoint_idx]
                Zx = Zes[xpoint_idx]
                if debug:
                    print('X-point (R, Z): ({:.2f}, {:.2f})'.format(Rx, Zx))

                # Use polar coordinates with the X-point as the origin.
                # The distance from the X-point is just the distance formula.
                rs = self.ts_dict['r']['Y']; zs = self.ts_dict['z']['Y']

                # If this is the first time(the reference frame), save these values.
                # These ref rs, zs are already (R, Z) in the reference frame (obviously),
                # so we don't need to convert to polar and then back to (R, Z) in the
                # reference frame.
                if ref_flag:
                    Rx_ref = Rx; Zx_ref = Zx
                    rs_ref = rs; zs_ref = zs
                    self.ref_df[str(time)] = list(zip(rs_ref, zs_ref))
                    ref_flag = False
                else:
                    # The angle is computed using arctan2 to get the correct quadrant (radians).
                    theta = np.arctan2(zs - Zx, rs - Rx)
                    #print("Theta: ", end=''); print(*theta, sep=', ')
                    d = np.sqrt((rs - Rx)**2 + (zs - Zx)**2)
                    #print("Distance from X-point: ", end=''); print(*d, sep=', ')

                    # Now convert back to (R, Z), but in the reference frame.
                    rs_into_ref = d * np.cos(theta) + Rx_ref
                    zs_into_ref = d * np.sin(theta) + Zx_ref
                    self.ref_df[str(time)] = list(zip(rs_into_ref, zs_into_ref))

                # Find the Te and ne data to put in for these times as well. Average
                # the neighboring +/-"average_ts" points. Ex. If average_ts = 5, and
                # the time is 2000, it will take the last 5 and next 5 TS points and
                # average the resulting 11 profiles together into one for 2000 ms.

                # First find closest time in the index.
                idx = np.abs(self.temp_df.index.values - time).argmin()

                # Average the desired range of values.
                #self.avg_temp_df = self.temp_df.iloc[idx - average_ts : idx + average_ts].mean()
                #self.avg_dens_df = self.dens_df.iloc[idx - average_ts : idx + average_ts].mean()
                self.avg_temp_df = self.temp_df.iloc[idx]
                self.avg_dens_df = self.dens_df.iloc[idx]

                # Append this to our ref_df.
                self.ref_df['Te at ' + str(time)] = self.avg_temp_df
                self.ref_df['Ne at ' + str(time)] = self.avg_dens_df

                # Let's also map each chord to R-Rsep omp and store in a DataFrame.
                # Get the additional information needed for this.
                Rs, Zs = np.meshgrid(gfile['R'], gfile['Z'])
                Z_axis = gfile['ZmAxis']
                R_axis = gfile['RmAxis']
                if trunc_div:
                    Rs_trunc = Rs > self.ts_dict['r']['Y'][0] * 0.95
                else:
                    Rs_trunc = Rs > R_axis

                # Only want the outboard half since thats where we're mapping R-Rsep OMP to.
                if choose_interp_region:
                    lcfs_rs = gfile['lcfs'][:, 0]
                    lcfs_zs = gfile['lcfs'][:, 1]
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    ax.plot(lcfs_rs, lcfs_zs, 'k.')
                    for i in range(0, len(lcfs_rs)):
                        ax.annotate(i, (lcfs_rs[i], lcfs_zs[i]))
                    ax.plot(rs, zs, 'r.')
                    fig.show()
                    start = int(input('Enter start index for interpolation region: '))
                    end   = int(input('Enter end index for interpolation region: '))
                    choose_interp_region = False

                Zes_outboard = np.copy(gfile['lcfs'][:, 1][start:end])
                Res_outboard = np.copy(gfile['lcfs'][:, 0][start:end])

                #else:
                #    Zes_outboard = np.copy(gfile['lcfs'][:, 1][13:-17])
                #    Res_outboard = np.copy(gfile['lcfs'][:, 0][13:-17])

                try:
                    # Interpolation functions of psin(R, Z) and R(psin, Z). Rs_trunc
                    # helps with not interpolating the entire plasma, and just that
                    # to the right of the magnetic axis, which is normally good enough.
                    f_psiN = Rbf(Rs[Rs_trunc], Zs[Rs_trunc], gfile['psiRZn'][Rs_trunc])
                    f_Romp = Rbf(gfile['psiRZn'][Rs_trunc], Zs[Rs_trunc], Rs[Rs_trunc], epsilon=0.00001)
                    f_Rs = interp1d(Zes_outboard, Res_outboard, assume_sorted=False)

                    # The process is to get the (R, Z) of each chord...
                    chord_rs = self.ts_dict['r']['Y']
                    chord_zs = self.ts_dict['z']['Y']

                    # ...then find the corresponding psin of each...
                    chord_psins = f_psiN(chord_rs, chord_zs)

                    # ...then use this psin to find the value at the omp (Z_axis)...
                    chord_omps = f_Romp(chord_psins, np.full(len(chord_psins), Z_axis))

                    # ... get the value of the separatrix at the omp the calculate
                    # R-Rsep omp...
                    rsep_omp = f_Rs(Z_axis)
                    chord_rminrsep_omps = chord_omps - rsep_omp

                    if debug:
                        print("Rsep OMP: {:.3f}".format(rsep_omp))
                        print("Chord OMPs: ", end=""); print(chord_omps)

                    # ...and then wrap each omp value with the corresponding temperature,
                    # and put each point into a dataframe.
                    te_idx = np.where(np.abs(self.temp_df.index.values-time) ==
                                      np.abs(self.temp_df.index.values-time).min())[0][0]
                    chord_tes = self.temp_df.iloc[te_idx]
                    chord_nes = self.dens_df.iloc[te_idx]
                    self.temp_df_omp[time] = list(zip(chord_rminrsep_omps, chord_tes))
                    self.dens_df_omp[time] = list(zip(chord_rminrsep_omps, chord_nes))

                    # Put the psin that the chord is on too.
                    self.temp_df_omp[str(time) + ' psin'] = list(zip(chord_psins, chord_tes))
                    self.dens_df_omp[str(time) + ' psin'] = list(zip(chord_psins, chord_nes))

                except Exception as e:
                    print("Error in the OMP steps: \n  " + str(e))
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    print(exc_type, fname, exc_tb.tb_lineno)

            except Exception as e:
                print("Error loading gfile: \n " + str(e))
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                print(exc_type, fname, exc_tb.tb_lineno)

            count += 1

        # End "for time in times" loop.

        # Create average Te and ne values mapped to the OMP.
        self.avg_omp = {}
        try:
            avg_psins = np.array([])
            avg_omps = np.array([])
            avg_tes  = np.array([])
            avg_nes  = np.array([])
            avg_omps_err = np.array([])
            avg_tes_err  = np.array([])
            avg_nes_err  = np.array([])
            for chord in range(0, len(self.temp_df_omp.index)):
                tmp_psins = np.array([])
                tmp_omps = np.array([])
                tmp_tes  = np.array([])
                tmp_nes  = np.array([])

                #for time in range(0, len(times)):
                for time in times:

                    # Get the tuple data point for this chord at this time (r-rsep_omp, Te).
                    tmp_p = self.temp_df_omp[str(time) + ' psin'].values[chord][0]
                    tmp_o = self.temp_df_omp[time].values[chord][0]
                    tmp_t = self.temp_df_omp[time].values[chord][1]
                    tmp_n = self.dens_df_omp[time].values[chord][1]
                    tmp_psins = np.append(tmp_psins, tmp_p)
                    tmp_omps  = np.append(tmp_omps, tmp_o)
                    tmp_tes   = np.append(tmp_tes,  tmp_t)
                    tmp_nes   = np.append(tmp_nes,  tmp_n)

                # Get the average for this chord. Append it to avg_omps/avg_tes.
                avg_psins = np.append(avg_psins, tmp_psins.mean())
                avg_omps  = np.append(avg_omps,  tmp_omps.mean())
                avg_tes   = np.append(avg_tes,   tmp_tes.mean())
                avg_nes   = np.append(avg_nes,   tmp_nes.mean())
                avg_omps_err = np.append(avg_omps_err, tmp_omps.std())
                avg_tes_err  = np.append(avg_tes_err, tmp_tes.std())
                avg_nes_err  = np.append(avg_nes_err, tmp_nes.std())

            # Store in dictionary in class.
            self.avg_omp['Psin']             = avg_psins
            self.avg_omp['RminRsep_omp']     = avg_omps
            self.avg_omp['Te_omp']           = avg_tes
            self.avg_omp['ne_omp']           = avg_nes
            self.avg_omp['RminRsep_omp_err'] = avg_omps_err
            self.avg_omp['Te_omp_err']       = avg_tes_err
            self.avg_omp['ne_omp_err']       = avg_nes_err

        except Exception as e:
            print("Error in calculating the average OMP values: \n  " + str(e))
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)



    def heatmap(self, te_clip=50, ne_clip=2e25, rlim_min=1.3, rlim_max=1.6,
                zlim_min=-1.3, zlim_max=-1.0, detach_front_te=5, offset=0.01,
                sp_hit_z=-1.25):
        """
        Function to produce the 2D maps of Te and ne on the poloidal cross
        section of ref_time. load_ts and map_to_efit must be run first to store
        the needed data into the class!

        te_clip:  Clips any Te data above this value. Use if you have problems
                    with high values messing up the countour plot.
        ne_clip:  Same, but for ne.
        rlim_min: Limits on the plotting windows.
        rlim_max: See above.
        zlim_min: See above.
        zlim_max: See above.
        te_lims:  Lower and upper limits for the contour plot levels. Won't plot
                    anything higher. Must be a tuple or None!!! Ex. (0, 100) will
                    only plot Te values in the range of 0-100 eV.
        ne_lims:  Same, but for ne.
        sp_hit_z: The shelf is at Z = -1.25, so the SP hits at this location. Can
                  change to the floor location if that comes around as needed.
        """

        # First make sure ref_time is in self.times.
        if self.ref_time in self.times:

            # Create empty arrays to hold Rs, Zs, Tes and Nes.
            rs  = np.zeros((len(self.ref_df.index), len(self.times)))
            zs  = np.zeros((len(self.ref_df.index), len(self.times)))
            tes = np.zeros((len(self.ref_df.index), len(self.times)))
            nes = np.zeros((len(self.ref_df.index), len(self.times)))

            idx = 0
            for time in self.times:
                rs[:, idx]  = [p[0] for p in self.ref_df[str(time)].values]
                zs[:, idx]  = [p[1] for p in self.ref_df[str(time)].values]
                tes[:, idx] = [p for p in self.ref_df['Te at ' + str(time)].values]
                nes[:, idx] = [p for p in self.ref_df['Ne at ' + str(time)].values]
                idx += 1

            # Perform clipping so it isn't skewed toward higher values.
            tes = np.clip(tes, 0, te_clip)
            nes = np.clip(nes, 0, ne_clip)

            # Store rs, zs, tes, nes in class for debugging.
            self.rs = rs; self.zs = zs; self.tes = tes; self.nes = nes

            # For reference in messing with the plotting limits.
            print('Te (min/max): ({:.2f}/{:.2f})'.format(tes.min(), tes.max()))
            print('Ne (min/max): ({:.2e}/{:.2e})'.format(nes.min(), nes.max()))

            # Function to plot wall with lcfs and strike point.
            def plot_shot(fig, ax):
                ax.plot(self.ref_gfile['wall'][:, 0], self.ref_gfile['wall'][:, 1], 'k')
                ax.plot(self.ref_gfile['lcfs'][:, 0], self.ref_gfile['lcfs'][:, 1], 'k-')
                ax.set_xlim([rlim_min, rlim_max])
                ax.set_ylim([zlim_min, zlim_max])
                ax.set_xlabel('R (m)')
                ax.set_ylabel('Z (m)')
                xpoint = np.where(self.ref_gfile['lcfs'][:,1] == self.ref_gfile['lcfs'][:,1].min())[0][0]

                # Left leg.
                sp_r1 = self.ref_gfile['lcfs'][:,0][xpoint-1]
                sp_z1 = self.ref_gfile['lcfs'][:,1][xpoint-1]
                sp_r0 = self.ref_gfile['lcfs'][:,0][xpoint]
                sp_z0 = self.ref_gfile['lcfs'][:,1][xpoint]
                m = (sp_z1 - sp_z0) / (sp_r1 - sp_r0)
                sp_rs = np.linspace(sp_r0, 1.2, 10)
                b = sp_z1 - m * sp_r1
                sp_zs = m * sp_rs + b
                ax.plot(sp_rs, sp_zs, 'k')

                # Right leg.
                sp_r1 = self.ref_gfile['lcfs'][:,0][xpoint+1]
                sp_z1 = self.ref_gfile['lcfs'][:,1][xpoint+1]
                sp_r0 = self.ref_gfile['lcfs'][:,0][xpoint]
                sp_z0 = self.ref_gfile['lcfs'][:,1][xpoint]
                m = (sp_z1 - sp_z0) / (sp_r1 - sp_r0)
                sp_rs = np.linspace(sp_r0, 1.5, 1000)
                b = sp_z1 - m * sp_r1
                sp_zs = m * sp_rs + b
                ax.plot(sp_rs, sp_zs, 'k')

                # Store these values for use later in detach_length part.
                self.sp_r0 = sp_r0; self.sp_z0 = sp_z0
                self.sp_r1 = sp_r1; self.sp_z1 = sp_z1
                self.sp_rs = sp_rs; self.sp_zs = sp_zs

            # Te plot.
            fig = plt.figure(figsize=(13,6))
            ax1 = fig.add_subplot(121)
            plot_shot(fig, ax1)
            cont1 = ax1.contourf(rs, zs, tes,
                                 levels=[0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25],
                                 cmap='inferno', extend='max')
            cbar1 = fig.colorbar(cont1)
            ax1.set_title('Te (eV)')
            cbar1.ax.set_ylabel('Te (eV)', size=24)
            ax1.scatter(rs, zs, c=tes, edgecolor='k', cmap='inferno')

            # Ne plot.
            ax2 = fig.add_subplot(122)
            plot_shot(fig, ax2)
            cont2 = ax2.contourf(rs, zs, nes, levels=10, cmap='viridis', extend='max')
            ax2.set_title('ne (m-3)')
            cbar2 = fig.colorbar(cont2)
            cbar2.ax.set_ylabel('ne (m-3)', size=24)

            # Title it and plot.
            fig.suptitle(str(self.shot) + " (raw)")
            fig.tight_layout()
            fig.show()

            # Same as above, but with the interpolated data.
            f_te = Rbf(self.rs, self.zs, self.tes, epsilon=1e-9)
            f_ne = Rbf(self.rs, self.zs, self.nes)
            RS, ZS = np.meshgrid(np.linspace(self.rs.min(), self.rs.max(), 100),
                                 np.linspace(self.zs.min(), self.zs.max(), 100))

            """

            fig = plt.figure()
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)

            plot_shot(fig, ax1)
            cont1 = ax1.contourf(RS, ZS, f_te(RS, ZS),
                                 levels=[0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25],
                                 cmap='inferno')
            cbar1 = fig.colorbar(cont1, extend='max')
            ax1.set_title('Te (eV)')
            cbar1.ax.set_ylabel('Te (eV)', size=24)

            plot_shot(fig, ax2)
            cont2 = ax2.contourf(RS, ZS, f_ne(RS, ZS), levels=10, cmap='viridis')
            cbar2 = fig.colorbar(cont2)
            ax2.set_title('ne (m-3)')
            cbar2.ax.set_ylabel('ne (m-3)', size=24)

            fig.suptitle(str(self.shot) + " (interpolated)")
            fig.tight_layout()
            fig.show()

            """

            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            plot_shot(fig, ax1)
            cont1 = ax1.contourf(RS, ZS, np.clip(f_te(RS, ZS), 0, te_clip),
                                 levels=[0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25],
                                 cmap='inferno', extend='max')
            cbar1 = fig.colorbar(cont1)
            cbar1.ax.set_ylabel('Te (eV)', size=24)
            ax1.set_xlabel('R (m)', fontsize=24)
            ax1.set_ylabel('Z (m)', fontsize=24)
            ax1.set_title(str(self.shot) + ' Te (eV)', fontsize=24)

            # Routine to find how far down the leg the detachment goes.
            # Detachment front defined here as detach_front_te.
            if detach_front_te:
                TES = f_te(RS, ZS)
                def find_detach_front(offset=0):

                    dist_along_leg = 0
                    leg_rs = np.array([]); leg_zs = np.array([])

                    # Here's the part that measures it from the floor up.
                    # Slope of SP.
                    m = (self.sp_z1 - self.sp_z0) / (self.sp_r1 - self.sp_r0)
                    b = self.sp_z1 - m * self.sp_r1

                    # Find find the R of this where it hits the divertor (this is
                    # just point-slope formula).
                    sp_hit_r = 1/m * (sp_hit_z - self.sp_z0) + self.sp_r0
                    self.sp_hit_r = sp_hit_r

                    # Create a line between the X-point and this point where the SP hits.
                    length_rs = np.linspace(sp_hit_r, self.sp_r0, 1000)
                    length_zs = m * length_rs + b

                    self.length_rs = length_rs
                    self.length_zs = length_zs

                    # Go up this leg one point at a time until Te is > 5 eV.
                    for i in range(1, len(length_rs)):
                        tmp_r = length_rs[i] + offset
                        tmp_z = length_zs[i]

                        # Find the nearest point in the (R, Z) grid. To do this, first
                        # find the distance between our current leg point and each
                        # grid point.
                        grid_dist = np.sqrt(np.abs(RS-tmp_r)**2 + np.abs(ZS-tmp_z)**2)

                        # Then find index of the minimum distance.
                        closest_grid_point = np.where(grid_dist == grid_dist.min())

                        # Add the length travelled to our running distance from X-point.
                        dist_along_leg += np.sqrt((length_rs[i] - length_rs[i-1])**2
                                                + (length_zs[i] - length_zs[i-1])**2)

                        # If we hit detach_front_te, break and our answer is dist_along_leg.
                        #print("{:5.2f}  {:5.2f}  {5.2f}".format(tmp_r, tmp_z, TES[closest_grid_point]))
                        if TES[closest_grid_point] > detach_front_te:
                            print("Offset: {:.2f} cm".format(offset*100))
                            print("Detachment front is {:.2f} cm from the floor.".format(dist_along_leg*100))
                            self.detach_length = dist_along_leg
                            break

                        # Add the leg point to the list for plotting after.
                        leg_rs = np.append(leg_rs, tmp_r)
                        leg_zs = np.append(leg_zs, tmp_z)


                        """
                        dist_along_leg = 0
                        leg_rs = np.array([]); leg_zs = np.array([])
                        for i in range(1, len(self.sp_rs)):

                            # This finds the distance starting from the X-point. May
                            # be better to instead look at it from the floor up instead.

                            # Get the point along the leg, starting at the first point
                            # past the X-point (i starts at 1).
                            tmp_r = self.sp_rs[i] + offset
                            tmp_z = self.sp_zs[i]

                            # Find the nearest point in the (R, Z) grid. To do this, first
                            # find the distance between our current leg point and each
                            # grid point.
                            grid_dist = np.sqrt(np.abs(RS-tmp_r)**2 + np.abs(ZS-tmp_z)**2)

                            # Then find index of the minimum distance.
                            closest_grid_point = np.where(grid_dist == grid_dist.min())

                            # Add the length travelled to our running distance from X-point.
                            dist_along_leg += np.sqrt((self.sp_rs[i] - self.sp_rs[i-1])**2
                                                    + (self.sp_zs[i] - self.sp_zs[i-1])**2)

                            # If we hit detach_front_te, break and our answer is dist_along_leg.
                            if TES[closest_grid_point] < detach_front_te:
                                print("Offset: {:.2f} cm".format(offset*100))
                                print("Detachment front is {:.2f} cm from the X-point.".format(dist_along_leg*100))
                                self.detach_length = dist_along_leg
                                break

                            # Add the leg point to the list for plotting after.
                            leg_rs = np.append(leg_rs, tmp_r)
                            leg_zs = np.append(leg_zs, tmp_z)

                            """
                    return leg_rs, leg_zs


                self.leg_rs1, self.leg_zs1 = find_detach_front(offset=0)
                self.leg_rs2, self.leg_zs2 = find_detach_front(offset=offset)

                ax1.plot(self.leg_rs1, self.leg_zs1, 'r-')
                ax1.plot(self.leg_rs2, self.leg_zs2, 'r-')
                fig.tight_layout()
                fig.show()

        else:
            print("Error: ref_time not one of the times in map_to_efit.")

    def filter_elms(self, method='simple', fs='FS04', avg_thresh=2, plot_it=True,
                    window_len=11):
        """
        NOTE: You probably should just use create_omfit_excel.py instead of
               this, since OMFITprofiles is superior to filtering ELMs and such.
        Method to filter ELM data. Replace Te, ne data that was taken during an
        ELM with either exclude or with a linear fit between the value before the ELM
        and the value at the end.

        method     : One of 'simple' or 'median'.
        fs         : Simple method. Which filterscope to get data from/
        avg_thresh : Simple method. Anything above the average filterscope
                      value * avg_thresh will be considered an ELM and data in
                      that time range will be filtered out of the Thomson data.
        plot_it    : Plot the filtered data.
        window_len : Size of window for filtering/smoothing method. Must be odd.
        """
        print("Warning: Are you sure you want to use this? Consider using " + \
              "the workflow via create_omfit_excel.py. It is way better at " + \
              "getting ELM filtered data. Check there or the GitHub for more info.")

        if window_len % 2 == 0:
            raise ValueError("window_len must be an odd number.")

        if method == 'simple':

            # First pull in the filterscope data to detect ELMs.
            fs_obj = gadata(fs, shot=self.shot, connection=self.conn)

            # Indices of ELMs.
            abv_avg = fs_obj.zdata > fs_obj.zdata.mean() * avg_thresh

            # Plot it just to show it makes sense.
            if plot_it:
                fig, ax = plt.subplots()
                ax.plot(fs_obj.xdata[~abv_avg], fs_obj.zdata[~abv_avg], 'k')
                ax.plot(fs_obj.xdata[abv_avg],  fs_obj.zdata[abv_avg],  'r')
                ax.set_xlabel('Time (ms)')
                ax.set_ylabel(fs)
                fig.tight_layout()
                fig.show()

            # Create list of pairs of times, where anything between each pair is
            # data to be filtered.
            bad_times = []
            in_bad = False
            for i in range(0, len(abv_avg)):
                if not in_bad:
                    if abv_avg[i] == True:
                        bad_start = fs_obj.xdata[i]
                        in_bad = True
                else:
                    if abv_avg[i] == False:
                        bad_end = fs_obj.xdata[i-1]
                        in_bad = False
                        bad_range = (bad_start, bad_end)
                        bad_times.append(bad_range)

            # Get the indices (or just True/False array) of rows in the
            # temp_df/dens_df to filter out.
            filter_times = np.full(len(self.temp_df.index), False)
            filter_times_ref = []
            for bt in bad_times:

                # For the temp_df/dens_df.
                elm_times = np.logical_and(self.temp_df.index > bt[0], self.temp_df.index < bt[1])
                filter_times = np.logical_or(filter_times, elm_times)

            # Index only the data that didn't fall in an ELM time range.
            self.temp_df_filt = self.temp_df.iloc[~filter_times]
            self.dens_df_filt = self.dens_df.iloc[~filter_times]

            return None

        elif method == 'median':

            # Get the median filter from scipy.
            from scipy.signal import medfilt

            # Identify all the nonzero points since we don't care about zeros.
            #te_nonzero = self.temp_df != 0
            #ne_nonzero = self.dens_df != 0

            # Perform a median filter on the data.
            #te_filt = medfilt(self.temp_df[te_nonzero], window_len)
            #ne_filt = medfilt(self.dens_df[ne_nonzero], window_len)
            te_filt = medfilt(self.temp_df, window_len)
            ne_filt = medfilt(self.dens_df, window_len)

            # Plot all the chords if you want.
            if plot_it:
                num_col = 5
                num_rows = int(np.ceil(len(self.temp_df.columns) / num_col))
                for df, filt in [(self.temp_df, te_filt), (self.dens_df, ne_filt)]:
                    fig, axs = plt.subplots(num_rows, num_col, sharex=True)
                    axs = axs.flatten()
                    for chord in df.columns:
                        axs[chord].plot(df.index, df[chord], '-k.')
                        axs[chord].plot(df.index, filt[:, chord], '-r.')
                        if chord % num_col == 0:
                            if df.equals(self.temp_df):
                                axs[chord].set_ylabel('Te (eV)')
                            elif df.equals(self.dens_df):
                                axs[chord].set_ylabel('ne (m-3)')
                        if (num_rows * num_col - chord) < num_col:
                            axs[chord].set_xlabel('Time (ms)')
                    fig.tight_layout()
                    fig.show()

            # Store the filtered data.
            #self.temp_df_filt = self.temp_df.replace(te_filt)
            #self.dens_df_filt = self.dens_df.replace(ne_filt)

        elif method in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:

            def smooth(x, window_len=11, window='hanning'):
                """
                This smoothing algorithm is taken form the scipy cookbook:
                scipy-cookbook.readthedocs.io/items/SignalSmooth.html

                See there for an explanation, but essentially it just smooths
                the data.
                """

                if x.ndim != 1:
                    raise ValueError("smooth only accepts 1 dimension arrays.")

                if x.size < window_len:
                    raise ValueError("Input vector needs to be bigger than window size.")

                if window_len<3:
                    return x

                if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                    raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


                s = np.r_[x[window_len-1:0:-1], x, x[-2:-window_len-1:-1]]
                #print(len(s))
                if window == 'flat': #moving average
                    w = np.ones(window_len, 'd')
                else:
                    w = eval('np.' + window + '(window_len)')

                y = np.convolve(w / w.sum(), s, mode='valid')
                #return y
                return y[int(window_len / 2):-int(window_len / 2)]

            te_filt = np.zeros(self.temp_df.values.shape)
            ne_filt = np.zeros(self.temp_df.values.shape)

            for chord in range(0, len(self.temp_df.columns)):
                chord_te_filt = smooth(self.temp_df[chord].values, window_len, method)
                chord_ne_filt = smooth(self.dens_df[chord].values, window_len, method)
                te_filt[:, chord] = chord_te_filt
                ne_filt[:, chord] = chord_ne_filt

        # This method applies a rolling average and then a rolling median to
        # filter out ELMs (or just really smooth it).
        elif method == 'average_median':

            # Transpose just so each row is the time series data. Will transpose
            # back at the end.
            te_filt = np.zeros(self.temp_df.values.shape).T
            ne_filt = np.zeros(self.temp_df.values.shape).T
            times = self.temp_df.index.values

            for chord in range(0, te_filt.shape[0]):

                # First calculate the rolling average.
                roll_avg_te = np.zeros(len(times))
                roll_avg_ne = np.zeros(len(times))
                for i in range(0, len(times)):
                    if (i < window_len) or len(times) - i < window_len:
                        te_point = self.temp_df[chord].values[i]
                        ne_point = self.dens_df[chord].values[i]
                    else:
                        te_point = np.average(self.temp_df[chord].values[i-int(window_len/2):i+int(window_len/2)])
                        ne_point = np.average(self.dens_df[chord].values[i-int(window_len/2):i+int(window_len/2)])

                    # Put into array.
                    roll_avg_te[i] = te_point
                    roll_avg_ne[i] = ne_point

                # Then same thing, just do median of the rolling average array.
                roll_med_te = np.zeros(len(times))
                roll_med_ne = np.zeros(len(times))
                for i in range(0, len(times)):
                    if (i < window_len) or len(times) - i < window_len:
                        te_point = self.temp_df[chord].values[i]
                        ne_point = self.dens_df[chord].values[i]
                    else:
                        te_point = np.median(roll_avg_te[i-int(window_len/2):i+int(window_len/2)])
                        ne_point = np.median(roll_avg_ne[i-int(window_len/2):i+int(window_len/2)])

                    # Put into array.
                    roll_med_te[i] = te_point
                    roll_med_ne[i] = ne_point

                # Finally put the filtered data for this chord into the filtered array.
                te_filt[chord] = roll_med_te
                ne_filt[chord] = roll_med_ne

            # Don't forget to transpose the data again.
            te_filt = te_filt.T
            ne_filt = ne_filt.T

        elif method == 'median_average':

            # Transpose just so each row is the time series data. Will transpose
            # back at the end.
            te_filt = np.zeros(self.temp_df.values.shape).T
            ne_filt = np.zeros(self.temp_df.values.shape).T
            times = self.temp_df.index.values

            for chord in range(0, te_filt.shape[0]):

                # First calculate the rolling average.
                roll_avg_te = np.zeros(len(times))
                roll_avg_ne = np.zeros(len(times))
                for i in range(0, len(times)):
                    if (i < window_len) or (len(times) - i < window_len):
                        te_point = self.temp_df[chord].values[i]
                        ne_point = self.dens_df[chord].values[i]
                    else:
                        te_point = np.median(self.temp_df[chord].values[i-int(window_len/2):i+int(window_len/2)+1])
                        ne_point = np.median(self.dens_df[chord].values[i-int(window_len/2):i+int(window_len/2)+1])

                    # Put into array.
                    roll_avg_te[i] = te_point
                    roll_avg_ne[i] = ne_point

                # Then same thing, just do median of the rolling average array.
                roll_med_te = np.zeros(len(times))
                roll_med_ne = np.zeros(len(times))
                for i in range(0, len(times)):
                    if (i < window_len) or (len(times) - i < window_len):
                        te_point = self.temp_df[chord].values[i]
                        ne_point = self.dens_df[chord].values[i]
                    else:
                        te_point = np.average(roll_avg_te[i-int(window_len/2):i+int(window_len/2)+1])
                        ne_point = np.average(roll_avg_ne[i-int(window_len/2):i+int(window_len/2)+1])

                    # Put into array.
                    roll_med_te[i] = te_point
                    roll_med_ne[i] = ne_point

                # Finally put the filtered data for this chord into the filtered array.
                te_filt[chord] = roll_med_te
                ne_filt[chord] = roll_med_ne

            # Don't forget to transpose the data again.
            te_filt = te_filt.T
            ne_filt = ne_filt.T

        # Store the filtered data.
        self.temp_df_filt = pd.DataFrame(te_filt, columns=self.temp_df.columns, index=self.temp_df.index)
        self.dens_df_filt = pd.DataFrame(ne_filt, columns=self.temp_df.columns, index=self.temp_df.index)

        # Plot the data to see how it matches.
        if plot_it:
            if self.system == 'core':
                fig, axs = plt.subplots(4, 5, sharex=True, figsize=(15, 10))
                for i in range(0, 20):
                    ax = axs.flatten()[i]
                    ax.plot(self.temp_df_unfiltered.index, self.temp_df_unfiltered[i], 'k-')
                    ax.plot(self.temp_df_filt.index, self.temp_df_filt[i], 'r-')
                    ax.set_xlim([0, 6000])
                    ax.set_ylim([0, 120])
                    ax.annotate(str(i), (0.9, 0.9), xycoords='axes fraction')
                fig.tight_layout()
                fig.show()
            elif self.system == 'divertor':
                fig, axs = plt.subplots(2, 4, sharex=True, figsize=(12, 5))
                for i in range(0, 8):
                    ax = axs.flatten()[i]
                    ax.plot(self.temp_df_unfiltered.index, self.temp_df_unfiltered[i], 'k-')
                    ax.plot(self.temp_df_filt.index, self.temp_df_filt[i], 'r-')
                    ax.set_xlim([0, 6000])
                    ax.set_ylim([0, 120])
                    ax.annotate(str(i), (0.9, 0.9), xycoords='axes fraction')
                fig.tight_layout()
                fig.show()
