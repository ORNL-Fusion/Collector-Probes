# This script pulls the divertor langmuir probes and puts them into a dictionary,
# among other things. It is essentially a python2 translation of the Matlab script
# get_lp. Just import this script and run the function get_dict_of_lps(shot).
#
# Author: Shawn Zamperini

import MDSplus as mds


def get_mds_active_probes(shot):
    """
    Get the probes that were active during the shot. Used in main function.
        shot: the shot you want
    """

    # MDSplus connection to atlas where the data is store on the "LANGMUIR" tree.
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
            print "Probe " + str(check_pnum) + " is MDS probe " + str(mds_pnum)
            mds_index.append(mds_pnum)
            found_probes.append(check_pnum)

    number_of_probes = len(found_probes)
    print "Found data for " + str(number_of_probes) + " probes."

    # Store in dictionary and return it.
    active = {}
    active["tmin"] = tmin
    active["tmax"] = tmax
    active["runid"] = runid
    active["probes"] = found_probes
    active["mds_index"] = mds_index

    return active

def get_mds_lp_data(shot, mds_index):
    """
    Get LP data for a single probe. Used in main function.
        shot: the shot you want
        mds_index: a number 1-85 that corresponds to the mds node. These do not
            match the probe number (which is PNUM).
    """

    # MDS connection required through atlas tunnel.
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

    #print "Data stored for probe " + str(lp_data["pnum"]) + " (MDS index " + str(mds_index) + ")."

    return lp_data


def get_dict_of_lps(shot):
    """
    Run this function to get the Langmuir probe data in a dictionary
    of dictionaries. Each entry will be all the probe data in the form
    of a dictionary (it just sound confusing in word it isn't really
    that weird).

    shot: the shot you want the data for.
    """

    # Get a dictionary with the probe active during this shot.
    active = get_mds_active_probes(shot)
    print ""

    # Get a dictionary of each probe data, then store it all in one big dictionary.
    lps = {}
    for mds_index in active["mds_index"]:
        lp_data = get_mds_lp_data(shot, mds_index)
        probe_name = "probe " + str(lp_data["pnum"])
        lps[probe_name] = lp_data
        print "Data stored for " + str(probe_name) + " (MDS index " + str(mds_index) + ")."

    return lps
