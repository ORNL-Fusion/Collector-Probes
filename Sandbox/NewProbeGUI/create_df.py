import pull_mdsplus   as pull
import pandas         as pd
import numpy          as np
import meas_locations as geo
import MDSplus        as mds
import itertools
from scipy import interpolate


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

    returns:    The requested gfile as a distionary.
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

def rbs_into_df(number, probe, start=2500, end=5000, step=500, remote=True, verbal=False):
    """
    Pulls RBS data from the MDSplus tree 'dp_probes' and puts it into a
    DataFrame ready for analysis. Require ssh to R2D2 if remote.

    number: Probe number.
    probe:  One of AD, AU, BD, BU, CD, CU.
    start:  Start of time that will be analyzed (i.e. the first gfile loaded).
    end:    End of time for analysis (i.e. the last gfile loaded).
    step:   Time step for the above.
    remote: Set to True if using outside the DIII-D network.

    returns: A DataFrame formatted and ready to be filled with data (R-Rsep,
             R-Rsep_omp, etc.)
    """

    # Create array of times to be sampled.
    times = np.arange(start, end, step)

    # Create MDSplus connection to R2D2.
    if verbal:
        print("Connecting to r2d2...")
    if remote:
        server = 'localhost'
    else:
        server = 'r2d2.gat.com'
    conn = pull.thin_connect(number, server=server)
    if verbal:
        print("Connection made to r2d2.")

    # Get shots probe was in for and Rprobe.
    shots  = pull.pull_shots(conn, probe, verbal=verbal)
    rprobe = pull.pull_rprobe(conn, probe, probe_corr=True, verbal=verbal)

    # Then pull the RBS data.
    rbs_dict = pull.pull_all_rbs(conn, number, probe, server=server, verbal=verbal)

    # Now prepare the DataFrame. Will have set of data at each time, at each
    # shot. So essentially len(times)*len(shots) DataFrames stacked together.
    rbs_df = pd.DataFrame(rbs_dict)

    # Want 'locs' as an index.
    rbs_df.set_index('locs', inplace=True)

    # Create set of DataFrames, len(times) of them, to be 'stacked' on top of each other.
    rbs_df = pd.concat(list(itertools.repeat(rbs_df, len(times))), keys=times, names=['times'])

    # Now do it again, except with shots.
    rbs_df = pd.concat(list(itertools.repeat(rbs_df, len(shots))), keys=shots, names=['shots'])

    if verbal:
        print("Data from r2d2 pulled.")

    return rbs_df, rprobe

def fill_in_rbs_df(rbs_df, probe, rprobe, remote=True, verbal=False):
    """
    Takes the rbs_df from above and fill it in with R-Rsep, R-Rsep_omp, etc. It
    returns all if it, so that it may then be averaged and get the std. dev. of
    after all the data colloction has taken place. Requires ssh to atlas if remote.

    rbs_df: The DataFrame returned from rbs_into_df.
    probe:  One of AD, AU, BD, BU, CD, CU.
    rprobe: Radial position of probe tip returned from rbs_into_df.
    remote: Set to to True if using outside DIII-D network.

    returns: Filled in rbs_df.
    """

    if verbal:
        print("Analyzing atlas relevant data...")

    # Get the shots, times and locs from the rbs_df index.
    shots = np.unique(rbs_df.index.get_level_values('shots').values)
    times = np.unique(rbs_df.index.get_level_values('times').values)
    locs  = np.unique(rbs_df.index.get_level_values('locs').values)

    # Extra columns to be filled out.
    rbs_df['R-Rsep (cm)']     = pd.Series()
    rbs_df['R-Rsep omp (cm)'] = pd.Series()
    rbs_df['Psin']            = pd.Series()
    rbs_df['R (cm)']          = pd.Series()

    # Establish the Z to be used depending on the probe.
    if   probe[0] == 'A': Z_probe = -0.188
    elif probe[0] == 'B': Z_probe = -0.1546
    elif probe[0] == 'C': Z_probe = -0.2054
    else: print("Error in probe entry.")

    # Create MDSplus connection to atlas.
    if verbal:
        print("Connecting to atlas...")
    if remote:
        server = 'localhost'
    else:
        server='atlas.gat.com'
    conn = mds.Connection(server)
    if verbal:
        print("Connection made to atlas.")

    for shot in shots:
        for time in times:
            # Load gfile.
            gfile = load_gfile_mds(shot, time, connection=conn, tunnel=True)

            # Create grid of R's and Z's.
            Rs, Zs = np.meshgrid(gfile['R'], gfile['Z'])

            # Z and R of magnetic axis (where omp is), in m.
            Z_axis = gfile['ZmAxis']
            R_axis = gfile['RmAxis']

            # Z's and R's of the separatrix, in m.
            Zes = np.copy(gfile['lcfs'][:, 1][13:-17])
            Res = np.copy(gfile['lcfs'][:, 0][13:-17])

            # Only want right half of everything.
            Rs_trunc = Rs > R_axis

            # Interpolation functions of psin(R, Z) and R(psin, Z).
            f_psin = interpolate.Rbf(Rs[Rs_trunc], Zs[Rs_trunc], gfile['psiRZn'][Rs_trunc])
            f_Romp = interpolate.Rbf(gfile['psiRZn'][Rs_trunc], Zs[Rs_trunc], Rs[Rs_trunc], epsilon=0.00001)
            f_Rs   = interpolate.interp1d(Zes, Res, assume_sorted=False)

            # R of the separatrix at each probe Z in cm.
            Rsep     = f_Rs(Z_probe) * 100.0
            Rsep_omp = f_Rs(Z_axis)  * 100.0

            # Get R of each location along the probe in cm, then R-Rsep.
            R_locs   = geo.calc_R_meas(rprobe, locs, probe)
            RminRsep = R_locs - Rsep

            # Get the corresponding psins of each location along the probe.
            psin_locs = f_psin(R_locs / 100.0, np.full((len(R_locs),), Z_probe))

            # Calculate R_loc at the omp, then R-Rsep omp.
            R_locs_omp   = f_Romp(psin_locs, np.full((len(psin_locs),), Z_axis)) * 100.0
            RminRsep_omp = R_locs_omp - Rsep_omp

            # Finally store all these in the corresponding part of the DataFrame.
            rbs_df.loc[shot].loc[time]['R-Rsep (cm)']     = pd.Series(RminRsep,     index=rbs_df.loc[shot].loc[time].index)
            rbs_df.loc[shot].loc[time]['R-Rsep omp (cm)'] = pd.Series(RminRsep_omp, index=rbs_df.loc[shot].loc[time].index)
            rbs_df.loc[shot].loc[time]['Psin']            = pd.Series(psin_locs,    index=rbs_df.loc[shot].loc[time].index)
            rbs_df.loc[shot].loc[time]['R (cm)']          = pd.Series(R_locs,       index=rbs_df.loc[shot].loc[time].index)

    return rbs_df

def rbs_df_stats(rbs_df, verbal=False):
    """
    Computes the average of each data point at each location along the probe.

    rbs_df:  DataFrame returned from the above 'fill_in_rbs_df'.

    returns: DataFrame of averages at each location for each time during each
             shot.
    """

    if verbal:
        print("Aggregating statistics over all shots and times...")

    # First get how many locations there are.
    locs  = np.unique(rbs_df.index.get_level_values('locs').values)
    nlocs = locs.size

    # The DataFrames that will hold our results.
    rbs_stat_df = pd.DataFrame()
    err_df      = pd.DataFrame()

    # To understand the indexing here, it'd be best to get this DataFrame into
    # a terminal and see how it works in there. It shows the beauty of pandas.
    for idx in range(0, nlocs):
        # Get the mean values at each location.
        rbs_stat_df = rbs_stat_df.append(rbs_df[idx::nlocs].mean(axis=0), ignore_index=True)
        # Get the standard deviations at each location.
        err_df      = err_df.append(rbs_df[idx::nlocs].std(axis=0), ignore_index=True)

    # Rename columns to appropriate names. The last two are already errors in rbs_stat_df,
    # so std. dev. of them isn't really a thing. Trash them.
    err_df.columns = ['Psin Error', 'trash1', 'R-Rsep Error (cm)',
                      'R-Rsep omp Error (cm)', 'trash2', 'trash3']
    err_df.drop(['trash1', 'trash2', 'trash3'], axis=1, inplace=True)

    # Put into one DataFrame to return.
    rbs_stat_df = rbs_stat_df.join(err_df)

    # Add locs in just because.
    rbs_stat_df['Distance from Tip (cm)'] = locs

    # Fix a couple column names. Sort Columns.
    rbs_stat_df.rename(columns={'areal':'W Areal Density (1e15 W/cm2)',
                                'areal_err':'W Areal Density Error (1e15 W/cm2)'},
                                inplace=True)
    rbs_stat_df = rbs_stat_df.sort_index(axis=1)

    return rbs_stat_df

def get_rbs(number, probe, start=2500, end=5000, step=500, remote=True, verbal=False):
    """
    This function can be considered a wrapper for the above functions. Use this
    function to get the rbs_df that you most likely want.

    number: The probe number.
    probe:  One of AD, AU, BD, BU, CD, CU.
    start:  Start time for time range to average R-Rsep and such over (inclusive).
    end:    End time for time range to average R-Rsep and such over (inclusive).
    step:   Time step for the above.
    remote: Set to True if accessing from outside DIII-D network.
    verbal: Set to True if you want feedback as program runs.

    returns: DataFrame with plasma coordinates of each location along the probe.
    """

    if remote:
        input("SSH link r2d2 to localhost. Press enter to continue...")
    # Get r2d2 data df, and rprobe (with correction already applied).
    rbs_df, rprobe = rbs_into_df(number, probe, start, end, step, remote, verbal)

    if remote:
        input("SSH link atlas to localhost. Press enter to continue...")
    # Get atlas related data.
    rbs_df = fill_in_rbs_df(rbs_df, probe, rprobe, remote, verbal)

    # Do statistics, creating much smaller df.
    stat_df = rbs_df_stats(rbs_df, verbal)

    return stat_df