import MDSplus as mds
import numpy   as np


def thin_connect(shot, tree='dp_probes', server='r2d2.gat.com', verbal=False):
    """
    Return thin connection to DP_PROBES MDSplus tree on R2D2.
    shot: The probe number (not a real DIII-D shot!).
    """
    if verbal:
        print("Retrieving MDSplus connection...")
    conn = mds.Connection(server)
    conn.openTree(tree, shot)
    return conn

def pull_rprobe(conn, probe, probe_corr=True, verbal=False):
    """
    Pull radial position of probe in cm.
    probe_corr: The probes were actually inserted 1.6 cm further out, apply
                this correction here.
    """
    if verbal:
        print("Retrieving rprobe...")
    path = '\\DP_PROBES::TOP.' + probe[0] + ':RPROBE'
    r_probe = conn.get(path).data()
    if probe_corr:
        r_probe += 1.6
    return r_probe

def pull_shots(conn, probe, verbal=False):
    """ Pull shots the probe was in for. """
    if verbal:
        print("Retrieving shots...")
    if probe in ['AU', 'AD', 'BU', 'BD', 'CU', 'CD']:
        path = '\\DP_PROBES::TOP.' + probe[0] + ':shots'
        shots = np.array(conn.get(path).data())
        return shots

def pull_rbs_loc(conn, probe, verbal=False):
    """ Pull locations of RBS data. Locations in mm, convert to cm here. """
    if verbal:
        print("Retrieving RBS locations...")
    if probe in ['AU', 'AD', 'BU', 'BD', 'CU', 'CD']:
        rbs_locs = np.array([])
        for run in range(25):
            if run < 10:
                run_str = '0' + str(run)
            else:
                run_str = str(run)
            path = '\\DP_PROBES::TOP.' + probe[0] + '.' + probe + '.RBS.RUN' + run_str + ':LOC'
            try:
                loc = conn.get(path).data() / 10.0
                rbs_locs = np.append(rbs_locs, loc)
                if verbal:
                    #print("  Found data for run " + str(run) + ".")
                    pass
            except:
                pass
        return rbs_locs
    else:
        print("Incorrect probe entry.")
        return None

def pull_rbs_areal(conn, probe, verbal=False):
    """ Pull W areal density at each location. Areal density in 1e15 W/cm2. """
    if verbal:
        print("Retrieving RBS areal density...")
    if probe in ['AU', 'AD', 'BU', 'BD', 'CU', 'CD']:
        rbs_areal = np.array([])
        for run in range(25):
            if run < 10:
                run_str = '0' + str(run)
            else:
                run_str = str(run)
            path = '\\DP_PROBES::TOP.' + probe[0] + '.' + probe + '.RBS.RUN' + run_str + ':W_AREAL'
            try:
                areal = conn.get(path).data()
                rbs_areal = np.append(rbs_areal, areal)
                if verbal:
                    #print("  Found data for run " + str(run) + ".")
                    pass
            except:
                pass
        return rbs_areal
    else:
        print("Incorrect probe entry.")
        return None

def pull_rbs_areal_err(conn, probe, verbal=False):
    """ Pull W areal density error. Units in 1e15 W/cm2. """
    if verbal:
        print("Retrieving RBS areal density error...")
    if probe in ['AU', 'AD', 'BU', 'BD', 'CU', 'CD']:
        rbs_areal_err = np.array([])
        for run in range(25):
            if run < 10:
                run_str = '0' + str(run)
            else:
                run_str = str(run)
            path = '\\DP_PROBES::TOP.' + probe[0] + '.' + probe + '.RBS.RUN' + run_str + ':W_AREAL_ERR'
            try:
                areal_err = conn.get(path).data()
                rbs_areal_err = np.append(rbs_areal_err, areal_err)
                if verbal:
                    #print("  Found data for run " + str(run) + ".")
                    pass
            except:
                pass
        return rbs_areal_err
    else:
        print("Incorrect probe entry.")
        return None

def pull_all_rbs(conn, shot, probe, verbal=False):
    """
    Pull RBS data and return it all in a dictionary.
    conn:  An MDSplus connection. Get from the above 'thin_connect' function.
    probe: One of AD, AU, BD, BU, CD, CU.
    """

    rbs_dict = {}
    #rbs_dict['rprobe']    = pull_rprobe(conn, probe, probe_corr=True, verbal=verbal)
    #rbs_dict['shots']     = pull_shots(conn, probe, verbal)
    rbs_dict['locs']      = pull_rbs_loc(conn, probe, verbal)
    rbs_dict['areal']     = pull_rbs_areal(conn, probe, verbal)
    rbs_dict['areal_err'] = pull_rbs_areal_err(conn, probe, verbal)

    return rbs_dict

def pull_lams(conn, shot, probe, verbal=False):
    """
    Pull the LAMS data and return it as a dictionary.
    """

    if verbal:
        print("Retrieving " + probe + str(shot) + " LAMS data...")
    path = '\\DP_PROBES::TOP.' + probe[0] + '.' + probe + '.LAMS'

    lams_dict = {}
    lams_dict['scandist_' + probe[1] + ' (mm)'] = conn.get(path + ':SCANDIST').data()
    for iso in ['180', '182', '183', '184', '186']:
        lams_dict['EF' + iso + '_' + probe[1]] = conn.get(path + ':EF' + iso).data()
        lams_dict['EF' + iso + '_err_' + probe[1]] = conn.get(path + ':EF' + iso + '_ERR').data()

    return lams_dict
