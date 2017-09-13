from MDSplus import *
import openpyxl as xl

tree = Tree('dp_probes', -1,'NEW')

# Create A,B,C childs.
print "Creating A,B,C's..."
letters = ['A', 'B', 'C']
for letter in letters:
	path = '\\DP_PROBES::TOP' + '.' + letter
	tree.addNode(path, 'structure')

# Add shots in for for each probe.
print "Creating shots..."
for letter in letters:
	path = '\\DP_PROBES::TOP' + '.' + letter + ':' + 'shots'
	tree.addNode(path, 'numeric')

# Create child for each AU,AD,BU,... probe.
print "Creating probes..."
probes = ['A.AU', 'A.AD','B.BU', 'B.BD', 'C.CU', 'C.CD']
#probes = ['A', 'B', 'C']
for probe in probes:
	path = '\\DP_PROBES::TOP' + '.' + probe
	tree.addNode(path, 'structure')

# Add RBS child for each probe type.
print "Creating RBS..."
for probe in probes:
	path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'RBS'
	tree.addNode(path, 'structure')

# Create structure to hold "raw" RBS data, leter to be used in building the signals of RUN01-RUN21.
#print "Creating RBS raw..."
#for probe in probes:
#	path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'RBS' + '.' + 'RAW'
#	tree.addNode(path, 'structure')

# Create 21 "raw" numerics to hold arrays of RBS data, 512 entries each.
#print "Creating RBS raw entries..."
#for probe in probes:
#	for entry in range(1,22):
#		if entry < 10:
#			path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'RBS' + '.' + 'RAW' + ':' + 'RAW' + '0' + str(entry)
#		else:
#			path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'RBS' + '.' + 'RAW' + ':' + 'RAW' + str(entry)
#
#		tree.addNode(path, 'numeric')

# Add 21 runs to RBS (each will have 512 channels).
print "Creating RBS run structures..."
for probe in probes:
	for run in range(1,22):
		if run < 10:
			path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'RBS' + '.' + 'run' + '0' + str(run)
		else:
			path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'RBS' + '.' + 'run' + str(run)
		tree.addNode(path,'structure')

print "Creating RBS run information..."
for probe in probes:
	for run in range(1,22):
		if run < 10:
			run_str = '0' + str(run)
		else:
			run_str = str(run)

		path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'RBS' + '.' + 'run' + run_str + ':' + 'signal'
		tree.addNode(path, 'signal')

		path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'RBS' + '.' + 'run' + run_str + ':' + 'loc'
		tree.addNode(path, 'numeric')
		
		# W Counts is the sum of channels 420 - 440. 
		path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'RBS' + '.' + 'run' + run_str + ':' + 'w_counts'
		tree.addNode(path, 'numeric')
		#wc_node = tree.getNode(path)
		#expr = Data.compile("SUM(" + probe + ".RBS.RUN" + run_str + ":SIGNAL[420..440])")
		#wc_node.putData(expr)

		path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'RBS' + '.' + 'run' + run_str + ':' + 'microcol'
		tree.addNode(path, 'numeric')

		path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'RBS' + '.' + 'run' + run_str + ':' + 'w_areal'
		tree.addNode(path, 'numeric')

# Add the extra refs and RCX data for RBS.
print "Creating RBS refs and RCX nodes..."
print "Loading excel file... This could take a minute."
wb = xl.load_workbook("RBS_excel_file.xlsx", data_only=True)

path = '\\DP_PROBES::TOP.RBS_REFS'
tree.addNode(path, 'structure')
path = '\\DP_PROBES::TOP.RBS_REFS:AU_COUNTS'
tree.addNode(path, 'numeric')
sheet = wb['refs']
au_count = sheet['B516'].value
my_node = tree.getNode(path)
my_node.putData(au_count)

path = '\\DP_PROBES::TOP.RBS_REFS:AU_MICROCOL'
tree.addNode(path, 'numeric')
au_micro = sheet['B517'].value
my_node = tree.getNode(path)
my_node.putData(au_micro)

path = '\\DP_PROBES::TOP.RBS_REFS:AU_AREAL'
tree.addNode(path, 'numeric')
au_areal = sheet['B518'].value
my_node = tree.getNode(path)
my_node.putData(au_areal)


path = '\\DP_PROBES::TOP.RBS_RCX'
tree.addNode(path, 'structure')
path = '\\DP_PROBES::TOP.RBS_RCX:W_RCX'
tree.addNode(path, 'numeric')
sheet = wb['RCX']
w_rcx = sheet['I5'].value
my_node = tree.getNode(path)
my_node.putData(w_rcx)

path = '\\DP_PROBES::TOP.RBS_RCX:AU_RCX'
tree.addNode(path, 'numeric')
au_rcx = sheet['I6'].value
my_node = tree.getNode(path)
my_node.putData(au_rcx)



# Add ICPMS and LAMS for each probe type.
print "Creating ICPMS and LAMS..."
for probe in probes:
	path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'ICPMS'
	tree.addNode(path,'structure')
	path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'LAMS'
	tree.addNode(path,'structure')

# For each probe type (AU, AD, ...) add six locations in ICPMS.
print "Creating ICPMS locations..."
locations = ['loc1','loc2','loc3','loc4','loc5', 'loc6']
for probe in probes:
	for location in locations:
		path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'ICPMS' + '.' + location
		tree.addNode(path,'structure')
		path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'ICPMS' + '.' + location + ':' + 'position'
		tree.addNode(path,'numeric')

# Add mass spectrums and standards to each location in ICPMS.
print "Creating ICPMS spectrums..."
spectrums = ['spectrum1', 'spectrum2', 'spectrum3']
for probe in probes:
	for location in locations:
		path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'ICPMS' + '.' + location + '.' + 'standards'
		tree.addNode(path,'structure')
		for number in range(1,6):
			path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'ICPMS' + '.' + location + '.' + 'standards.standard' + str(number)
			tree.addNode(path, 'structure')
			path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'ICPMS' + '.' + location + '.' + 'standards.standard' + str(number) + ':data'
			tree.addNode(path, 'signal')
			path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'ICPMS' + '.' + location + '.' + 'standards.standard' + str(number) + ':conc'
			tree.addNode(path, 'numeric')
			path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'ICPMS' + '.' + location + '.' + 'standards.standard' + str(number) + ':coeff'
			tree.addNode(path, 'numeric')
			

		for spectrum in spectrums:
			path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'ICPMS' + '.' + location + '.' + spectrum  
			tree.addNode(path,'structure')

# Add concentration of each spectrum and the signal itself.
print "Creating ICPMS signals..."
for probe in probes:
	for location in locations:
		for spectrum in spectrums:
			path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'ICPMS' + '.' + location + '.' + spectrum + ':' + 'conc'
			tree.addNode(path,'numeric')
			path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'ICPMS' + '.' + location + '.' + spectrum + ':' + 'coeff'
			tree.addNode(path,'numeric')
			path = '\\DP_PROBES::TOP' + '.' + probe + '.' + 'ICPMS' + '.' + location + '.' + spectrum + ':' + 'data'
			tree.addNode(path,'signal')

# Save the tree
print "Saving tree..."
tree.write()
print "Done."
