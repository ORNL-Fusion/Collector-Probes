# Python3 script to pull Thomson Scattering data from MDSplus. See function
# at bottom called run_script for general use.
#
# Author: Shawn Zamperini

import MDSplus as mds
import numpy as np
import scipy.interpolate as scinter
import math


def load_gfile_mds(shot, time, tree="EFIT01", exact=False, connection=None, tunnel=True):
    """
    This is scavenged from the load_gfile_d3d script on the EFIT repository,
    except updated it to run on python3. Also took out the write2file part since
    I never use it.
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


def get_ts(shot, system="core", tunnel=True):
    """
    Function to get all the data from the Thomson Scattering "BLESSED" tree
    on atlas.

    shot:   The shot number
    system: Either core, divertor  or tangential.
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

    # Open the tree of the shot we want.
    tree = conn.openTree("d3d", shot)

    # Until someone desires otherwise, default to the BLESSED revision.
    base = "\\D3D::TOP.ELECTRONS.TS.BLESSED"

    # Specify which system we want.
    print("Thomson system: " + system)
    if system is "core":
        base = base + ".CORE"
    elif system is "divertor":
        base = base + ".DIVERTOR"
    elif system is "tangential":
        base = base + ".TANGENTIAL"

    # Dictionary to hold all the TS data. To be returned at the end.
    ts_dict = {}

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
            ts_dict[node.lower()] = data_dict

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
                        ts_dict[node.lower() + "." + subnode.lower()] = data_dict
                    except (mds.MdsIpException, mds.TreeNODATA):
                        print("  Node has no data.")

        # This error is returned if the node is empty. Catch it.
        except (mds.MdsIpException, mds.TreeNODATA):
            print("  Node has no data.")

    ts_dict["system"] = system
    ts_dict["shot"]   = shot
    return ts_dict


def assign_psins(ts_dict, tmin=2500, tmax=5000, tstep=100, tree="EFIT01"):
    """
    This function is here to calculate psins of each chord. It will put them in
    the dictionary passed to it (ts_dict). There will be psin values for each
    chord, and for each chord at each time in the range (tmin, tmax, tstep).

    ts_dict: The TS dictionary returned from the function get_ts above.
    tmin: The minimum time range to start calculating what the psins were.
    tmax: The maximum time range.
    tstep: Time step to take between tmin and tmax.
    """

    # Get the times at which the Te and ne data was taken for each chord.
    times = ts_dict["temp"]["X"]
    range_wanted = np.arange(tmin, tmax, tstep)

    # Get index where the times start and end for the requested range.
    min_index = np.argmax(times>tmin)
    max_index = np.argmax(times>tmax)

    print("Time range used: [" + str(times[min_index]) + ", " + str(times[max_index]) + "]")

    # Get the actual min and max times. Fill array of times for gfiles to request.
    ts_min_time = times[min_index]
    ts_max_time = times[max_index]
    times_for_gfile = np.arange(ts_min_time, ts_max_time, tstep)
    num_of_times = times_for_gfile.size

    # Parameters needed from ts_dict.
    shot = ts_dict["shot"]
    ts_Rs = ts_dict["r"]["Y"]
    ts_Zs = ts_dict["z"]["Y"]
    channels = ts_dict["channel"]["Y"]
    temp = ts_dict["temp"]["Y"]
    dens = ts_dict["density"]["Y"]

    # Create connection here then pass it to gfile function so you don't have to
    # create a new connection every time. Saves a lot of time.
    conn = mds.Connection("localhost")

    # 2D array to hold all the psin values. The index of each row corresponds to
    # a certain channel/chord.
    all_psins = np.zeros((num_of_times, len(channels)))
    all_tes   = np.zeros((num_of_times, len(channels)))
    all_nes   = np.zeros((num_of_times, len(channels)))
    all_omps  = np.zeros((num_of_times, len(channels)))
    row_index = 0

    # Get the gfiles one time slice at a time.
    for time in times_for_gfile:

        # Load gfile. Cast to int since gfiles don't have decimal times.
        gfile = load_gfile_mds(shot, int(time), tree=tree, connection=conn)

        # Create grid of R's and Z's.
        Rs, Zs = np.meshgrid(gfile['R'], gfile['Z'])
        Z_axis = gfile['ZmAxis']
        R_axis = gfile['RmAxis']
        Zes = np.copy(gfile['lcfs'][:, 1][13:-17])
        Res = np.copy(gfile['lcfs'][:, 0][13:-17])
        Rs_trunc = Rs > R_axis

        # Interpolation functions of psin(R, Z) and R(psin, Z).
        f_psiN = scinter.Rbf(Rs[Rs_trunc], Zs[Rs_trunc], gfile['psiRZn'][Rs_trunc])
        f_Romp = scinter.Rbf(gfile['psiRZn'][Rs_trunc], Zs[Rs_trunc], Rs[Rs_trunc], epsilon=0.00001)
        f_Rs = scinter.interp1d(Zes, Res, assume_sorted=False)

        # R value of separatrix at omp.
        rSep_omp = f_Rs(Z_axis)

        # Temporary array to hold psins.
        psins = np.array([])
        tes   = np.array([])
        nes   = np.array([])
        rminrsep_omp = np.array([])
        for index in range(0, len(channels)):

            # The R and Z of each channel.
            tmp_R = ts_Rs[index]
            tmp_Z = ts_Zs[index]

            # This is the psin at this specific time.
            tmp_psin = f_psiN(tmp_R, tmp_Z)

            # R at the omp of this psin.
            tmp_Romp = f_Romp(tmp_psin, Z_axis)

            # This is the corresponding Te at this time.
            te_index = np.argmax(times>time)
            tmp_te   = temp[index][te_index]
            tmp_ne   = dens[index][te_index]

            # Add it to psins/tes, where in order it is the channels.
            psins = np.append(psins, tmp_psin)
            tes = np.append(tes, tmp_te)
            nes = np.append(nes, tmp_ne)

            # The R-Rsep_omp position.
            tmp_romp = tmp_Romp - rSep_omp
            rminrsep_omp = np.append(rminrsep_omp, tmp_romp)

        # Put row into all the psins. So a 2D array, where each row is the psins/Tes
        # of every channel at a time.
        all_psins[row_index] = psins
        all_tes[row_index]   = tes
        all_nes[row_index]   = nes
        all_omps[row_index]  = rminrsep_omp

        row_index = row_index + 1

    def organize_2d(all_vals):
        """
        Function to reorganize the 2D array daya of Te, ne and psin into an
        array where each row corresponds to the data from a specific chord
        (normally each row is the data from all the chords at a specific time).
        """
        all_vals_org     = np.zeros((len(channels), num_of_times))
        avg_vals_org     = np.zeros(len(channels))
        avg_vals_org_err = np.zeros(len(channels))
        row_index   = 0
        for chord in range(0, len(channels)):
            chord_vals = np.array([])
            for vals_at_time in all_vals:
                chord_vals_at_time = vals_at_time[chord]
                chord_vals = np.append(chord_vals, chord_vals_at_time)
            avg_chord_vals     = np.mean(chord_vals)
            avg_chord_vals_err = np.std(chord_vals)
            avg_vals_org[row_index]     = avg_chord_vals
            avg_vals_org_err[row_index] = avg_chord_vals_err
            all_vals_org[row_index]     = chord_vals
            row_index = row_index + 1
        return all_vals_org, avg_vals_org, avg_vals_org_err

    # Get the organized 2D arrays of the psin, Te and ne data.
    all_psins_org, avg_psins_org, avg_psins_org_err = organize_2d(all_psins)
    all_tes_org,   avg_tes_org,   avg_tes_org_err   = organize_2d(all_tes)
    all_nes_org,   avg_nes_org,   avg_nes_org_err   = organize_2d(all_nes)
    all_omps_org,  avg_omps_org,  avg_omps_org_err  = organize_2d(all_omps)

    # Put these bad boys into the dictionary under psin to keep it organized.
    psins_dict = {"channels":channels,           "all_psins":all_psins_org,
                  "avg_psins":avg_psins_org,     "avg_Tes":avg_tes_org,
                  "avg_nes":avg_nes_org,         "avg_Tes_err":avg_tes_org_err,
                  "avg_nes_err":avg_nes_org_err, "avg_psins_err":avg_psins_org_err,
                  "all_Tes":all_tes_org,         "all_nes":all_nes_org,
                  "all_omps":all_omps_org,       "avg_omps":avg_omps_org,
                  "avg_omps_err":avg_omps_org_err,
                  "times":times_for_gfile}

    ts_dict["psins"] = psins_dict
    return ts_dict


def run_script(shot, system, tunnel=True, psins=True, tmin=2500, tmax=5000, tstep=250, tree="EFIT01"):
    """
    This is the master function and likely the one most user will want. It will
    return a dictionary of all the TS info. The data is organized as dictionaries
    inside this return dictionary. In each data dict, there is an X and a Y. X
    is whatever the dimension data is (most the time just the channel/chord number),
    and Y is the actual data.
    The temp and density data is a bit more complex: the X is the time, and the
    Y is a 2D array where each row is all the data for a particular chord. For
    example, to get time vs. Te for chord #6, you would plot the X array
    against the 6 row of the Y array.

    shot:   The shot TS data is requested for.
    system: Either "core", "tangential" or "divertor".
    tunnel: Use True if not running inside the DIII-D network (i.e. locally). This
            requires an SSH connection linking atlas.gat.com to your localhost.
            An example command doing this is:
            ssh -Y -p 2039 -L 8000:atlas.gat.com:8000 username@cybele.gat.com
            where this is run in a new terminal.
    psins:  Set to True if you want to know the psin values of the Te and ne
            measurements. This can take time, so if you don't need these values
            just set it to False.
    tmin:   The minimum time range to start calculating what the psins were.
    tmax:   The maximum time range to start calculating what the psins were.
    tstep:  Time step to take between tmin and tmax.
    tree:   Whichever EFIT tree you want to pull from.
    """

    # Get the TS data from MDSplus.
    ts_dict = get_ts(shot, system, tunnel)

    # Calculate and include the psin data if desired.
    if psins:
        ts_dict = assign_psins(ts_dict, tmin, tmax, tstep, tree)

    return ts_dict
