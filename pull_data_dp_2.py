# Functions to pull data from MDSplus collector probe tree. The
# user functions are at the bottom. These return all the above 
# info as a dictionary.



import MDSplus as mds

def thin_connect(shot, tree='dp_probes'):
	conn = mds.Connection('r2d2.gat.com')
	conn.openTree(tree, shot)
	tree = mds.Tree(tree, shot)
	return conn

def get_tree(shot, tree='dp_probes'):
	tree = mds.Tree(tree, shot)
	return tree


# RBS relevant data.

# Note this requires the tree, not the connection.
def pull_rbs_raw(tree, probe, run):
	if probe in ['AU', 'AD', 'BU', 'BD', 'CU', 'CD']:
		if run < 10:
			run_str = '0' + str(run)
		else:
			run_str = str(run)
		path = '\\DP_PROBES::TOP.' + probe[0] + '.' + probe + '.RBS.RUN' + run_str + ':SIGNAL'
		sig_node = tree.getNode(path)
		raw_data = sig_node.getData().raw_of()
		return raw_data

def pull_rbs_wCounts(conn, probe, run):
	if probe in ['AU', 'AD', 'BU', 'BD', 'CU', 'CD']:
		if run < 10:
			run_str = '0' + str(run)
		else:
			run_str = str(run)
		path = '\\DP_PROBES::TOP.' + probe[0] + '.' + probe + '.RBS.RUN' + run_str + ':W_COUNTS'
		wcounts = conn.get(path).data()
		return wcounts
	else:
		print "Incorrect probe entry."

def pull_shots(conn, probe):
	if probe in ['AU', 'AD', 'BU', 'BD', 'CU', 'CD']:
		path = '\\DP_PROBES::TOP.' + probe[0] + ':shots'
		shots = conn.get(path).data()
		return shots

def pull_rbs_loc(conn, probe, run):
	if probe in ['AU', 'AD', 'BU', 'BD', 'CU', 'CD']:
		if run < 10:
			run_str = '0' + str(run)
		else:
			run_str = str(run)
		path = '\\DP_PROBES::TOP.' + probe[0] + '.' + probe + '.RBS.RUN' + run_str + ':LOC'
		try:
			loc = conn.get(path).data()
			return loc
		except:
			print "Data for run " + str(run) + " not available."
	else:
		print "Incorrect probe entry."

def pull_rbs_microcol(conn, probe, run):
	if probe in ['AU', 'AD', 'BU', 'BD', 'CU', 'CD']:
		if run < 10:
			run_str = '0' + str(run)
		else:
			run_str = str(run)
		path = '\\DP_PROBES::TOP.' + probe[0] + '.' + probe + '.RBS.RUN' + run_str + ':MICROCOL'
		micro = conn.get(path).data()
		return micro
	else:
		print "Incorrect probe entry."

def pull_rbs_areal(conn, probe, run):
	if probe in ['AU', 'AD', 'BU', 'BD', 'CU', 'CD']:
		if run < 10:
			run_str = '0' + str(run)
		else:
			run_str = str(run)
		path = '\\DP_PROBES::TOP.' + probe[0] + '.' + probe + '.RBS.RUN' + run_str + ':W_AREAL'
		areal = conn.get(path).data()
		return areal
	else:
		print "Incorrect probe entry."

def pull_rbs_ref_au_counts(conn):
	path = '\\DP_PROBES::TOP.RBS_REFS:AU_COUNTS'
	au_counts = conn.get(path).data()
	return au_counts

def pull_rbs_ref_au_microcol(conn):
	path = '\\DP_PROBES::TOP.RBS_REFS:AU_MICROCOL'
	au_microcol = conn.get(path).data()
	return au_microcol

def pull_rbs_ref_au_areal(conn):
	path = '\\DP_PROBES::TOP.RBS_REFS:AU_AREAL'
	au_areal = conn.get(path).data()
	return au_areal

def pull_rbs_rcx_au(conn):
	path = '\\DP_PROBES::TOP.RBS_RCX:AU_RCX'
	au_rcx = conn.get(path).data()
	return au_rcx


def pull_rbs_rcx_w(conn):
	path = '\\DP_PROBES::TOP.RBS_RCX:W_RCX'
	w_rcx = conn.get(path).data()
	return w_rcx




# ICPMS relvant data.

# This one doesn't work (yet).
def pull_standard_conc(conn, probe, loc, stan):
	if probe in ['AU', 'AD', 'BU', 'BD', 'CU', 'CD']:
		path = '\\DP_PROBES::TOP.' + probe[0] + '.' + probe + '.ICPMS.LOC' + str(loc) + '.STANDARDS.STANDARD' + str(stan) + ':CONC'
		stan_conc = conn.get(path).data()
		return stan_conc
	else:
		print "Incorrect probe entry."

def pull_conc(conn, probe, loc, spectrum):
	if probe in ['AU', 'AD', 'BU', 'BD', 'CU', 'CD']:
		path = '\\DP_PROBES::TOP.' + probe[0] + '.' + probe + '.ICPMS.LOC' + str(loc) + '.SPECTRUM' + str(spectrum) + ':CONC'
		conc = conn.get(path).data()
		return conc
	else:
		print "Incorrect probe entry."

def pull_icpms_pos(conn, probe, loc):
	if probe in ['AU', 'AD', 'BU', 'BD', 'CU', 'CD']:
		path = '\\DP_PROBES::TOP.' + probe[0] + '.' + probe + '.ICPMS.LOC' + str(loc) + ':POSITION'
		position = conn.get(path).data()
		return position
	else:
		print "Incorrect probe entry"

# Note these last two require the tree, not the connection.
def pull_icpms_signal(tree, probe, spectrum, loc):
	if probe in ['AU', 'AD', 'BU', 'BD', 'CU', 'CD']:
		path = '\\DP_PROBES::TOP.' + probe[0] + '.' + probe + '.ICPMS.LOC' + str(loc) + '.SPECTRUM' + str(spectrum) + ':DATA'
		sig_node = tree.getNode(path)
		raw_sig = sig_node.getData().raw_of()
		return raw_sig
	else:
		print "Incorrect probe entry"

def pull_standard_data(tree, probe, loc, stan):
	if probe in ['AU', 'AD', 'BU', 'BD', 'CU', 'CD']:
		path = '\\DP_PROBES::TOP.' + probe[0] + '.' + probe + '.ICPMS.LOC' + str(loc) + '.STANDARD:STANDARD' + str(stan)	
		sig_node = tree.getNode(path)
		raw_sig = sig_node.getData().raw_of()
		return raw_sig
	else:
		print "Incorrect probe entry."


# The functions you care about.

# Get all relevant rbs data in a dictionary. Not all probes have rbs data.
# 	probe        - one of AU, AD, BU, BD, CU, CD
#	probe_number - 1-37 (not all probes have data though)
# 	run          - 1-21 (not all runs have data)

def rbs_dict(probe, probe_number, run):
	if probe in ['AU', 'AD', 'BU', 'BD', 'CU', 'CD']:
		tree = get_tree(probe_number)
		conn = thin_connect(probe_number)
		rbs_dict = {'probe': probe, 'probe number': probe_number}

		# Which shots was the probe inserted for.
		rbs_dict['shots in for'] = pull_shots(conn, probe)
	
		# What is the location along the probe, in mm from the tip.
		rbs_dict['location'] = pull_rbs_loc(conn, probe, run)

		# W areal density in W/cm^2.
		rbs_dict['w_areal'] = pull_rbs_areal(conn, probe, run)

		# The raw rbs spectrum of how many counts per channel.
		rbs_dict['rbs raw spectrum'] = pull_rbs_raw(tree, probe, run)

		# The sum of channels 420-440, assumed to correspond to W.
		rbs_dict['rbs w counts'] = pull_rbs_wCounts(conn, probe, run)

		return rbs_dict
	
	else:
		print "Incorrect probe entry."

# Get all relevant icpms data.
# 	probe - one of AU, AD, BU, BD, CU, CD
#	probe_number - 1-37 (not all probes have data though)
# 	loc_number - 1-6, correspoinding to a section of the probe that was sampled.
# 	spectrum_number - 1-3, corresponding to the three concentrations used (i.e. 1 for 25 ppb, 2 for 50 ppb, 3 for 75 ppb)
# 	stan_number - 1-5, corresponding to which standard you want data for. Each standard has its own concentration.

def icpms_dict(probe, probe_number, loc_number, spectrum_number, stan_number):
	if probe in ['AU', 'AD', 'BU', 'BD', 'CU', 'CD']:
		tree = get_tree(probe_number)
		conn = thin_connect(probe_number)
		icpms_dict = {'probe': probe, 'probe number': probe_number}

		# Which shots the probe was inserted for.
		icpms_dict['shots in for'] = pull_shots(conn, probe)

		# The concentration of the particular spectrum.
		icpms_dict['conc'] = pull_conc(conn, probe, loc_number, spectrum_number)

		# The average position along the probe, from the tip, in mm.
		icpms_dict['position'] = pull_icpms_pos(conn, probe, loc_number)

		# Raw ICPMS spectrum used in another script to convert to counts per amu.
		icpms_dict['raw_spectrum'] = pull_icpms_signal(tree, probe, spectrum_number, loc_number)

		# Raw standard spectrum.
		icpms_dict['raw_standard'] = pull_standard_data(tree, probe, loc_number, stan_number)

		return icpms_dict

	else:
		print "Incorrect probe entry."


