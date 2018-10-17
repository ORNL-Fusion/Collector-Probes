import numpy   as np
import pandas  as pd
import MDSplus as mds
import matplotlib.pyplot as plt
from scipy import interpolate
from matplotlib import ticker, cm, colors
import matplotlib as mpl


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

        # Until someone desires otherwise, default to the BLESSED revision.
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
        nodes = ["CALCMASK", "CDPOLYBOX", "CDPULSE", "CHANNEL", "CHI_MAX", "CHI_MIN", "DCDATA", "DENSITY", "DENSITY_E",
                 "DETOPT", "DCPEDESTAL", "DETOPT", "FITDATA", "FITTHRESHOLD", "FRACCHI",
                 "INIT_NE", "INIT_TE", "ITMICRO", "LFORDER", "LPROF", "MAXFITS", "PHI",
                 "PLDATA", "PLERROR", "PLPEDESTAL", "PLPEDVAR", "R", "REDCHISQ", "SUPOPT", "TEMP",
                 "TEMP_E", "THETA", "TIME", "Z"]

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

    def load_gfile_mds(self, shot, time, tree="EFIT01", exact=False, connection=None, tunnel=True, verbal=True):
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
                raise RuntimeError(tree + ' does not exactly contain time %.2f' %time + '  ->  Abort')
            else:
                if verbal:
                    print('Warning: ' + tree + ' does not exactly contain time %.2f' %time + ' the closest time is ' + str(time0))
                    print('Fetching time slice ' + str(time0))
                time = time0

        # store data in dictionary
        g = {'shot': shot, 'time': time}

        # get header line
        header = connection.get(base + 'ECASE').data()[k]

        # get all signals, use same names as in read_g_file
        translate = {'MW': 'NR', 'MH': 'NZ', 'XDIM': 'Xdim', 'ZDIM': 'Zdim', 'RZERO': 'R0',
                     'RMAXIS': 'RmAxis', 'ZMAXIS': 'ZmAxis', 'SSIMAG': 'psiAxis', 'SSIBRY': 'psiSep',
                     'BCENTR': 'Bt0', 'CPASMA': 'Ip', 'FPOL': 'Fpol', 'PRES': 'Pres',
                     'FFPRIM': 'FFprime', 'PPRIME': 'Pprime', 'PSIRZ': 'psiRZ', 'QPSI': 'qpsi',
                     'NBBBS': 'Nlcfs', 'LIMITR': 'Nwall'}
        for signal in translate:
            g[translate[signal]] = connection.get(base + signal).data()[k]

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
        for item in ['Xdim', 'Zdim', 'R0', 'R1', 'RmAxis', 'ZmAxis', 'psiAxis', 'psiSep', 'Bt0', 'Ip']:
            g[item] = np.round(np.float64(g[item]), 7)

        # convert single arrays (float32) to double arrays (float64)
        for item in ['Fpol', 'Pres', 'FFprime', 'Pprime', 'psiRZ', 'qpsi', 'lcfs', 'wall']:
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

    def map_to_efit(self, times=None, average_ts=5):
        """
        This function will use EFIT to map the machine coordinates to psin and
        rho coordinates. Really easy to do with pandas.

        times: A list of times to average over and map to. Can be a single float
               or list.
        """

        # Store times for later use in plotting function.
        self.times = np.array(times)

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

        # Current implementation makes the first time the reference frame.
        ref_flag = True
        for time in times:
            print('Loading gfile (' + str(count) + '/' + str(len(times)) + ')...')
            gfile = self.load_gfile_mds(shot=self.shot, time=time, connection=self.conn, verbal=True)

            # Create grid of R's and Z's.
            #Rs, Zs = np.meshgrid(gfile['R'], gfile['Z'])

            # Z and R of magnetic axis (where omp is), in m.
            #Z_axis = gfile['ZmAxis']
            #R_axis = gfile['RmAxis']

            # Z's and R's of the separatrix, in m.
            Zes = np.copy(gfile['lcfs'][:, 1])
            Res = np.copy(gfile['lcfs'][:, 0])

            # Only want right half of everything, otherwise we're asking too much
            # from the interpolation coming up. It should be accurate to assume
            # the TS system is to the right of the magnetic axis, simplifying
            # the interpolation. IGNORE FOR NOW. MIGHT BE UNNECESSARY.
            #Rs_trunc = Rs > R_axis

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

    def heatmap(self, ref_time, te_clip=80, ne_clip=2e20):

        # First make sure ref_time is in self.times.
        if ref_time in self.times:
            # Get the gfile for plotting the walls and such.
            gfile = self.load_gfile_mds(shot=self.shot, time=ref_time, connection=self.conn, verbal=True)

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
            #tes = np.clip(tes, 0, te_clip)
            #nes = np.clip(nes, 0, ne_clip)

            fig = plt.figure()

            # Te plot.
            ax1 = fig.add_subplot(121)
            ax1.plot(gfile['wall'][:, 0], gfile['wall'][:, 1], 'k')
            ax1.plot(gfile['lcfs'][:, 0], gfile['lcfs'][:, 1], 'k-')
            #cont1 = ax1.contourf(rs, zs, tes, 6, cmap='Reds', locator=ticker.LogLocator(), vmin=1, vmax=100)
            cont1 = ax1.contourf(rs, zs, tes, levels=[0,10,20,30,40,50,60,70,80,90,100], cmap='Reds')
            ax1.set_xlim([1.2, 1.7])
            ax1.set_ylim([-1.4, -0.9])
            ax1.set_xlabel('R (m)')
            ax1.set_ylabel('Z (m)')
            ax1.set_title('Te (eV)')
            cbar1 = fig.colorbar(cont1)
            cbar1.ax.set_ylabel('Te (eV)', size=24)

            # Ne plot.
            ax2 = fig.add_subplot(122)
            ax2.plot(gfile['wall'][:, 0], gfile['wall'][:, 1], 'k')
            ax2.plot(gfile['lcfs'][:, 0], gfile['lcfs'][:, 1], 'k-')
            #cont2 = ax2.contourf(rs, zs, nes, 6, cmap='Blues', locator=ticker.LogLocator(), vmin=10e19, vmax=10e21)
            cont2 = ax2.contourf(rs, zs, nes,
                                 levels=[1e19, 2.5e19, 5e19, 7.5e19, 1e20, 2.5e20, 5e20, 7.5e20, 1e21],
                                 cmap='Blues', norm=colors.LogNorm())
            ax2.set_xlim([1.2, 1.7])
            ax2.set_ylim([-1.4, -0.9])
            ax2.set_xlabel('R (m)')
            ax2.set_ylabel('Z (m)')
            ax2.set_title('ne (m-3)')
            cbar2 = fig.colorbar(cont2)
            cbar2.ax.set_ylabel('ne (m-3)', size=24)

            fig.tight_layout()
            fig.show()
        else:
            print("Error: ref_time not one of the times in map_to_efit.")
