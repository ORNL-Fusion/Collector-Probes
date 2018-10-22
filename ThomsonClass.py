import numpy   as np
import pandas  as pd
import MDSplus as mds
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors
import matplotlib as mpl


print("Use ThomsonClass object as:")
print("  ts = ThomsonClass.ThomsonClass(176343, 'divertor')")
print("  ts.load_ts()")
print("  ts.map_to_efit(times=np.linspace(2000, 5000, 5), ref_time=2000)")
print("  ts.heatmap()\n")
print("Note: ref_time must be in times. It will be the drawn profile. Can also")
print("just use for load_ts to load the TS data for whatever purpose.")

class ThomsonClass:
    """
    This object will be written to hold a set of functions, as well as
    dictionaries (or DataFrames) of the Thomson Scattering (TS) data. There
    will also be functions that can make 1D and 2D plots of the data.
    """

    def __init__(self, shot, system):
        self.shot   = shot
        self.system = system
        self.conn = None
        self.ts_dict = {}

    def __repr__(self):
        return "ThomsonClass Object\n  Shot:   " + str(self.shot) + "\n  System: "+ str(self.system)

    def load_ts(self, tunnel=True):
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

        # Get the data for each node, and put it in a dictionary of X and Y values.
        # For 1D data, like "Z", the X will be the channel number, and the Y will
        # be the Z coordinate.
        for node in nodes:
            try:
                # Make the path to the node.
                path = base + "." + node
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
                        print("Getting data from node: " + path)
                        try:
                            data = conn.get(path).data()
                            data_dim = conn.get("dim_of(" + path + ")").data()
                            data_dict = {"Time":data_dim, subnode.lower():data}
                            self.ts_dict[node.lower() + "." + subnode.lower()] = data_dict
                        except (mds.MdsIpException, mds.TreeNODATA):
                            print("  Node has no data.")

            # This error is returned if the node is empty. Catch it.
            except (mds.MdsIpException, mds.TreeNODATA):
                print("  Node has no data.")

            # For compatablity with some of the older shots.
            except mds.TreeNNF:
                print("  Node not found.")

    def load_gfile_mds(self, shot, time, tree="EFIT01", exact=False,
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
                    print('Warning: ' + tree + ' does not exactly contain time \
                          %.2f' %time + ' the closest time is ' + str(time0))
                    print('Fetching time slice ' + str(time0))
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

    def map_to_efit(self, times=None, ref_time=None, average_ts=5):
        """
        This function uses the gfiles of each time in times to find the location
        of each TS chord relative to the X-point in polar coordinates, (d, theta).
        By doing this for a swept strike point and mapping each (d, theta) back
        to a reference frame, 2D profiles of Te and ne can be obtained.

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

        # Store times for later use in plotting function.
        self.times = np.array(times)
        self.ref_time = ref_time

        if ref_time in times:

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

            # Will create a DataFrame where the index is the psin location and the data
            # is either the temp/dens at the requested time or averaged over the
            # times, if multiple times passed.

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
            for time in times:
                print('Loading gfile (' + str(count) + '/' + str(len(times)) + ')...')
                gfile = self.load_gfile_mds(shot=self.shot, time=time, connection=self.conn, verbal=True)

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
                self.avg_temp_df = self.temp_df.iloc[idx - average_ts : idx + average_ts].mean()
                self.avg_dens_df = self.dens_df.iloc[idx - average_ts : idx + average_ts].mean()

                # Append this to our ref_df.
                self.ref_df['Te at ' + str(time)] = self.avg_temp_df
                self.ref_df['Ne at ' + str(time)] = self.avg_dens_df

                count += 1

        else:
            print("Error: ref_time not in times.")

    def heatmap(self, te_clip=1e10, ne_clip=2e25, rlim_min=1.2, rlim_max=1.7,
                zlim_min=-1.4, zlim_max=-0.9, te_lims=None, ne_lims=None):
        """
        Function to produce the 2D maps of Te and ne on the pooloidal cross
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

            # For reference in messing with the plotting limits.
            print('Te (min/max): ({:.2f}/{:.2f})'.format(tes.min(), tes.max()))
            print('Ne (min/max): ({:.2e}/{:.2e})'.format(nes.min(), nes.max()))

            # Perform clipping so it isn't skewed toward higher values.
            tes = np.clip(tes, 0, te_clip)
            nes = np.clip(nes, 0, ne_clip)

            fig = plt.figure()

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
                sp_rs = np.linspace(sp_r0, 1.5, 10)
                b = sp_z1 - m * sp_r1
                sp_zs = m * sp_rs + b
                ax.plot(sp_rs, sp_zs, 'k')


            # Te plot.
            ax1 = fig.add_subplot(121)
            plot_shot(fig, ax1)
            #cont1 = ax1.contourf(rs, zs, tes, 6, cmap='Reds', locator=ticker.LogLocator(), vmin=1, vmax=100)
            if te_lims is None:
                cont1 = ax1.contourf(rs, zs, tes, levels=10, cmap='inferno')
            else:
                cont1 = ax1.contourf(rs, zs, tes, levels=np.linspace(te_lims[0], te_lims[1], 10), cmap='inferno')
            ax1.set_title('Te (eV)')
            cbar1 = fig.colorbar(cont1)
            cbar1.ax.set_ylabel('Te (eV)', size=24)

            # Ne plot.
            ax2 = fig.add_subplot(122)
            plot_shot(fig, ax2)
            #cont2 = ax2.contourf(rs, zs, nes, 6, cmap='Blues', locator=ticker.LogLocator(), vmin=10e19, vmax=10e21)
            if ne_lims is None:
                cont2 = ax2.contourf(rs, zs, nes, levels=10, cmap='viridis')
            else:
                cont2 = ax2.contourf(rs, zs, nes, levels=np.linspace(ne_lims[0], ne_lims[1], 10), cmap='viridis')
                                 #levels=[1e19, 2.5e19, 5e19, 7.5e19, 1e20, 2.5e20, 5e20, 7.5e20, 1e21],
                                 #cmap='Blues', norm=colors.LogNorm())
            ax2.set_title('ne (m-3)')
            cbar2 = fig.colorbar(cont2)
            cbar2.ax.set_ylabel('ne (m-3)', size=24)

            fig.tight_layout()
            fig.show()
        else:
            print("Error: ref_time not one of the times in map_to_efit.")
