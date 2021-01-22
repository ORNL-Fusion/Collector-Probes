# This script pulls the divertor langmuir probes and puts them into a dictionary,
# among other things. It is essentially a python translation of the Matlab script
# get_lp. Just import this script and run the function get_dict_of_lps(shot).
#
# Author: Shawn Zamperini
import MDSplus as mds
import numpy   as np
import pandas  as pd
import matplotlib.pyplot as plt


def get_mds_active_probes(shot, tunnel=True):
    """
    Get the probes that were active during the shot. Used in main function.

    shot: the shot you want
    """

    # MDSplus connection to atlas where the data is store on the "LANGMUIR" tree.
    if tunnel:
        conn = mds.Connection("localhost")
    else:
        conn = mds.Connection('atlas.gat.com')
    conn.openTree("LANGMUIR", shot)

    tmin = conn.get("\LANGMUIR::TOP.TMIN").data()
    tmax = conn.get("\LANGMUIR::TOP.TMAX").data()
    runid = conn.get("\LANGMUIR::TOP.RUNID").data()

    mds_index = []
    found_probes = []
    for mds_pnum in range(1,85):

        # Make sure probe name is in correct formart: 001, 002, ... , 084, 085.
        if mds_pnum < 10:
            probe = "00" + str(mds_pnum)
        else:
            probe = "0" + str(mds_pnum)

        # The complete path name to the lp. PNUM is the probe number, which does
        # not match its number in mdsplus (001 - 085).
        pname = "\LANGMUIR::TOP.PROBE_" + probe + ".PNUM"

        # Get the actual probe number if it is there. Not all MDS probes are used.
        try:
            check_pnum = conn.get(pname).data()
        except:
            pass
            #print "No data in probe " + str(probe) + "."

        # It will be '0' or blank if the MDS entry isn;t used. Otherwise it will
        # have the actual probe number in it.
        if check_pnum > 0:
            #print("Probe " + str(check_pnum) + " is MDS probe " + str(mds_pnum))
            mds_index.append(mds_pnum)
            found_probes.append(check_pnum)

    number_of_probes = len(found_probes)
    print("Found data for " + str(number_of_probes) + " probes.")

    # Store in dictionary and return it.
    active = {}
    active["tmin"] = tmin
    active["tmax"] = tmax
    active["runid"] = runid
    active["probes"] = found_probes
    active["mds_index"] = mds_index

    return active

def get_mds_lp_data(shot, mds_index, tunnel=True):
    """
    Get LP data for a single probe. Used in main function.

    shot: the shot you want
    mds_index: a number 1-85 that corresponds to the mds node. These do not
      match the probe number (which is PNUM).
    """

    # MDS connection required through atlas tunnel.
    if tunnel:
        conn = mds.Connection("localhost")
    else:
        conn = mds.Connection("atlas.gat.com")
    conn.openTree("LANGMUIR", shot)

    # Use correct form of probe name.
    if mds_index < 10:
        probe = "00" + str(mds_index)
    else:
        probe = "0" + str(mds_index)

    pname = "\LANGMUIR::TOP.PROBE_" + probe

    # All the data stored in a dictionary. All the data is in the subtree
    # indicated in pname. Just specify the node and grab the data.
    lp_data = {}
    lp_data["time"]       = conn.get(pname + ":TIME").data()
    lp_data["rprobe"]     = conn.get(pname + ":R").data()
    lp_data["zprobe"]     = conn.get(pname + ":Z").data()
    lp_data["label"]      = conn.get(pname + ":LABEL").data()
    lp_data["ntimes"]     = conn.get(pname + ":NTIMES").data()
    lp_data["pnum"]       = conn.get(pname + ":PNUM").data()
    lp_data["isat"]       = conn.get(pname + ":ISAT").data()
    lp_data["jsat"]       = conn.get(pname + ":JSAT").data()
    lp_data["temp"]       = conn.get(pname + ":TEMP").data()
    lp_data["dens"]       = conn.get(pname + ":DENS").data()
    lp_data["pot"]        = conn.get(pname + ":POT").data()
    lp_data["psin"]       = conn.get(pname + ":PSIN").data()
    lp_data["angle"]      = conn.get(pname + ":ANGLE").data()
    lp_data["area"]       = conn.get(pname + ":AREA").data()
    lp_data["delrsepout"] = conn.get(pname + ":DELRSEPOUT").data()
    lp_data["delrsepin"]  = conn.get(pname + ":DELRSEPIN").data()
    lp_data["delzsepout"] = conn.get(pname + ":DELZSEPOUT").data()
    lp_data["delzsepin"]  = conn.get(pname + ":DELZSEPIN").data()
    lp_data["csq"]        = conn.get(pname + ":CSQ").data()
    lp_data["res_err"]    = conn.get(pname + ":RES_ERR").data()
    lp_data["heatflux"]   = conn.get(pname + ":HEATFLUX").data()
    lp_data["pnum"]       = conn.get(pname + ":PNUM").data()

    # Include an estimate of the ground or SOL current.
    with np.errstate(all='ignore'):
        lp_data["ground_j"]   = lp_data["jsat"] * (1 - np.exp(-lp_data["pot"]/lp_data["temp"]))

    #print "Data stored for probe " + str(lp_data["pnum"]) + " (MDS index " + str(mds_index) + ")."

    return lp_data


def get_dict_of_lps(shot, tunnel=True):
    """
    Run this function to get the Langmuir probe data in a dictionary
    of dictionaries. Each entry will be all the probe data in the form
    of a dictionary (it just sound confusing in word it isn't really
    that weird).

    shot: the shot you want the data for.
    """

    # Get a dictionary with the probe active during this shot.
    active = get_mds_active_probes(shot, tunnel=tunnel)
    print("")

    # Get a dictionary of each probe data, then store it all in one big dictionary.
    lps = {}
    for mds_index in active["mds_index"]:
        lp_data = get_mds_lp_data(shot, mds_index, tunnel=tunnel)
        probe_name = "probe " + str(lp_data["pnum"])
        lps[probe_name] = lp_data
        print("Data stored for " + str(probe_name) + " (MDS index " + str(mds_index) + ").")

    return lps

def plot_lps(shot, tmin, tmax, xtype='rminrsep', xlim=None, filter='median',
             bins=5, tunnel=True, csv_path=None):
    """
    Plot LP data, with optional filtering applied.

    shot (int): Shot you want data for.
    tmin (float): Start time for data.
    tmax (float): End time for data.
    xtype (str): X-axis for plot. One of "rminrsep", "psin" or "time".
    xlim ([list, float]): X limits for the plot. Entered as list or tuple,
      e.g. (-0.99, 1.4).
    filter (str): One of "median" or "average". How to treat the binned LP data.
    bins (int): Number of bins to divide the LP data up into. Divides an LP
      signal in time up into number of bins and ten performs the filtering on
      each one.
    tunnel (bool): Whether to tunnel through atlas or not. If True, require ssh
      linking atlas ot localhost.
    csv_path (str): Optional path to save data to as a csv file.
    """

    # Load lp data.
    lps = get_dict_of_lps(shot, tunnel)

    # Arrays to hold final values for plotting.
    all_te       = np.array([])
    all_ne       = np.array([])
    all_jsat     = np.array([])
    all_x        = np.array([])
    all_ground   = np.array([])

    # Output arrays.
    pnames_filt   = np.array([])
    x_filt        = np.array([])
    te_filt       = np.array([])
    ne_filt       = np.array([])
    jsat_filt     = np.array([])
    heatflux_filt = np.array([])
    r_filt        = np.array([])
    z_filt        = np.array([])
    ground_filt   = np.array([])

    # Go through one probe at a time to get data for plotting.
    for key in lps.keys():

        # Get times and restrict to given time range.
        times = lps[key]['time']
        idx = np.logical_and(times > tmin, times < tmax)
        times = times[idx]

        # Get Te, ne and jsat. Also get R-Rsep out.
        te       = lps[key]['temp'][idx]
        ne       = lps[key]['dens'][idx]
        jsat     = lps[key]['jsat'][idx]
        heatflux = lps[key]['heatflux'][idx]
        ground   = lps[key]['ground_j'][idx]

        if xtype == 'rminrsep':
            x = lps[key]['delrsepout'][idx]
        elif xtype == 'psin':
            x = lps[key]['psin'][idx]
        elif xtype == 'time':
            x = lps[key]['time'][idx]
        else:
            print("Error in xtype entry. Must be either rminrsep or psin.")

        # Put into one array so we can sort them all by rminrsep, low -> high.
        probe_data = np.array((x, te, ne, jsat, heatflux, ground))
        probe_data = probe_data[:, np.argsort(probe_data[0])]

        # Divide the data up into the number of 'bins', then take the filter of each bin.
        bin_size = len(probe_data[0]) / bins
        #print("Bin size = {:.2f} --> Using integer = {}".format(bin_size, int(bin_size)))
        bin_size = int(bin_size)

        for bin in range(0, bins):
            tmp_x        = probe_data[0][bin_size*bin:bin_size*(bin+1)]
            tmp_te       = probe_data[1][bin_size*bin:bin_size*(bin+1)]
            tmp_ne       = probe_data[2][bin_size*bin:bin_size*(bin+1)]
            tmp_jsat     = probe_data[3][bin_size*bin:bin_size*(bin+1)]
            tmp_heatflux = probe_data[4][bin_size*bin:bin_size*(bin+1)]
            tmp_ground   = probe_data[5][bin_size*bin:bin_size*(bin+1)]

            # Apply the preferred filter.
            if filter == 'median':
                #filter = np.median
                filter = np.nanmedian
            elif filter == 'average':
                filter = np.mean
                #filter = np.nanmean

            # Filter and add to the output arrays to be plotted.
            x_filt        = np.append(x_filt,        filter(tmp_x))
            te_filt       = np.append(te_filt,       filter(tmp_te))
            ne_filt       = np.append(ne_filt,       filter(tmp_ne))
            jsat_filt     = np.append(jsat_filt,     filter(tmp_jsat))
            heatflux_filt = np.append(heatflux_filt, filter(tmp_heatflux))
            ground_filt   = np.append(ground_filt,   filter(tmp_ground))

            # Assign probe names so we can identify these data points later.
            pnames_filt   = np.append(pnames_filt, key)
            r_filt        = np.append(r_filt, lps[key]['rprobe'])
            z_filt        = np.append(z_filt, lps[key]['zprobe'])


    # Enumerate the pnames so they can be used for color selection.
    pnames_enum = np.array(list(enumerate(np.unique(pnames_filt))))

    # General plotting commands. These are the "Tableau 20" colors as RGB.
    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
                 (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
                 (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
                 (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
                 (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

    # A nice looking font.
    plt.rcParams['font.family'] = 'serif'

    # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
    for i in range(len(tableau20)):
        r, g, b = tableau20[i]
        tableau20[i] = (r / 255., g / 255., b / 255.)

    # Create figure.
    fig = plt.figure(figsize=(15,7))

    # Function for plotting.
    def plot_ax(fig, x, y, ylabel, ax_num, high_y, xlabel, xlim, legend=False):
        ax = fig.add_subplot(ax_num)

        # For each data point assign correct color.
        for i in range(0, len(x)):
            for pnames_pair in pnames_enum:
                if pnames_pair[1] == pnames_filt[i]:
                    color = int(pnames_pair[0])
                    label = pnames_pair[1]

            ax.plot(x[i], y[i], '^', ms=10, color=tableau20[color], label=label.title())
            ax.set_xlabel(xlabel, fontsize=18)
            ax.set_ylabel(ylabel, fontsize=18)
            ax.axvline(0.0, linestyle='--', color='k')
            ax.set_xlim(xlim)
            ax.set_ylim([0, high_y])

            # Process to remove duplicate legend entries.
            if legend:
                handles, labels = plt.gca().get_legend_handles_labels()
                newLabels, newHandles = [], []
                for handle, label in zip(handles, labels):
                  if label not in newLabels:
                    newLabels.append(label)
                    newHandles.append(handle)
                ax.legend(newHandles, newLabels, framealpha=0.5)

    if xtype == 'rminrsep':
        xlabel = 'R-Rsep (m)'
        if xlim is None:
            xlim = [-0.05, 0.1]
    elif xtype == 'psin':
        xlabel = 'Psin'
        if xlim is None:
            xlim = [0.98, 1.1]
    elif xtype == 'time':
        xlabel = 'Time (ms)'


    plot_ax(fig, x_filt, te_filt, 'Te (eV)', 131, 50, xlabel, xlim, legend=True)
    plot_ax(fig, x_filt, ne_filt, 'ne (cm-3)', 132, 10e13, xlabel, xlim)
    plot_ax(fig, x_filt, jsat_filt, 'jsat (A/cm2)', 133, 100, xlabel, xlim)

    fig.tight_layout()
    fig.show()

    lp_dict = {xtype:x_filt, 'Te (eV)':te_filt, 'ne (cm-3)':ne_filt,
           'jsat (A/cm2)':jsat_filt, 'heatflux (W/cm2)':heatflux_filt,
           'ground_filt':ground_filt, 'pnames':pnames_filt, 'R':r_filt, 'Z':z_filt}

    if csv_path != None:
        df = pd.DataFrame(lp_dict)
        df.to_csv(csv_path)

    return lp_dict
