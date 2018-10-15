import numpy   as np
import pandas  as pd
import MDSplus as mds
import matplotlib.pyplot as plt
from scipy.interpolate import Rbf


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

    def load_gfile_mds(shot, time, tree="EFIT01", exact=False, connection=None, tunnel=True):
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
        print("\nLoading gfile:")
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

    def map_to_efit(self, times=None):
        """
        This function will use EFIT to map the machine coordinates to psin and
        rho coordinates. Really easy to do with pandas.

        times: A list of times to average over and map to. Can be a single float
               or list.
        """

        # Pull these into DataFrames, a logical and easy way to represent the
        # data. Initially the rows are each a TS chord, and the columns are at
        # each time.

        temp_df = pd.DataFrame(columns=ts.ts_dict['temp']['X'],    data=ts.ts_dict['temp']['Y'])
        dens_df = pd.DataFrame(columns=ts.ts_dict['density']['X'], data=ts.ts_dict['density']['Y'])
        temp_df.index.name   = 'Chord'
        temp_df.columns.name = 'Time (ms)'
        dens_df.index.name   = 'Chord'
        dens_df.columns.name = 'Time (ms)'

        # Transpose the data so each row is at a specific time, and the columns
        # are the chords.
        temp_df = temp_df.tranpose()
        dens_df = dens_df.tranpose()

        # Will create a DataFrame where the index is the psin location and the data
        # is either the temp/dens at the requested time or averaged over the
        # times, if multiple times passed.

        # Load gfile(s).
        psin_df = pd.DataFrame()
        for time in times:
            gfile = load_gfile_mds(self.shot, time, connection=self.conn)

            # Create grid of R's and Z's.
            Rs, Zs = np.meshgrid(gfile['R'], gfile['Z'])

            # Z and R of magnetic axis (where omp is), in m.
            Z_axis = gfile['ZmAxis']
            R_axis = gfile['RmAxis']

            # Z's and R's of the separatrix, in m.
            Zes = np.copy(gfile['lcfs'][:, 1][13:-17])
            Res = np.copy(gfile['lcfs'][:, 0][13:-17])

            # Only want right half of everything, otherwise we're asking too much
            # from the interpolation coming up. It should be accurate to assume
            # the TS system is to the right of the magnetic axis, simplifying
            # the interpolation.
            Rs_trunc = Rs > R_axis

            # Interpolation functions of psin(R, Z) and R(psin, Z).
            f_psin = interpolate.Rbf(Rs[Rs_trunc], Zs[Rs_trunc], gfile['psiRZn'][Rs_trunc])
            f_Romp = interpolate.Rbf(gfile['psiRZn'][Rs_trunc], Zs[Rs_trunc], Rs[Rs_trunc], epsilon=0.00001)
            #f_Rs   = interpolate.interp1d(Zes, Res, assume_sorted=False)

            # We know the (R, Z) of each chord, so use f_psin to get the psin of
            # it, then store in DataFrame. Will do the averaging and such after
            # the loop for the sake of being modular. For now, just get the data.
            rs = self.ts_dict['r']['Y']; zs = self.ts_dict['z']['Y']
            psins = f_psin(r, z)
            psin_df[time] = psins
