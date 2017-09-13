import openpyxl as xl
from MDSplus import *
import sys

# Used to get location of .scn files. 
from Tkinter import Tk
from tkFileDialog import askopenfilename




def get_RBS(tree, letter_probes, shot):
	# Grab the RBS data. 
	print "\nLoading RBS Excel file... This may take a minute."
	rbs_file = xl.load_workbook("RBS_excel_file.xlsx", data_only=True)
	print "RBS Excel file loaded."
	rbs_probe_list = rbs_file.get_sheet_names()

	# Remove unecessary sheets.
	rbs_probe_list.remove('refs')
	rbs_probe_list.remove('RCX')
	rbs_probe_list.remove('SUMMARY')
	rbs_probe_list.remove('Sheet6')

	# Check if RBS data available for the selected probes.
	for letter_probe in letter_probes:
		tmp_name = letter_probe + 'U' + str(shot)
		if (tmp_name not in rbs_probe_list):
			print 'RBS data not available for ' + tmp_name + '.' 
		tmp_name = letter_probe + 'D' + str(shot)
		if (tmp_name not in rbs_probe_list):
			print 'RBS data not available for ' + tmp_name + '.' 

	# Collect data from Excel sheet and put them into a signal.
	for letter_probe in letter_probes:
		for u_or_d in ['U', 'D']:
			name = letter_probe + u_or_d + str(shot)

			# Pass through if there isn't RBS data.
			if (name not in rbs_probe_list): continue

			print "Assigning RBS data to " + name + " probe..."

			# Grab the corresponding RBS sheet.
			sheet = rbs_file.get_sheet_by_name(name)

			# Fill in run data, microcoul, w_counts and w_areal density.
			count = 0
			for row in 'BCDEFGHIJKLMNOPQRSTUV':
				count = count + 1
				if count < 10:
					count_str = '0' + str(count)
				else:
					count_str = str(count)
			
				# Run data.
				rbs_cells = sheet[row + '2': row + '513']
				rbs_vals = []
				for index in range(0,512):
					rbs_vals.append(rbs_cells[index][0].value)
				# If "NoneType" (i.e. blank cell), skip over.
				if (rbs_vals[0] is None): 
					print "Column " + row + " blank."			
					continue
			
				path = '\\DP_PROBES::TOP.' + letter_probe + '.' + letter_probe + u_or_d + '.RBS.RUN' + count_str + ':' + 'SIGNAL'	
				my_node = tree.getNode(path)
				#sig_expr = Data.compile("BUILD_SIGNAL($VALUE, BUILD_WITH_UNITS($1,'COUNTS'), \
				#			 BUILD_WITH_UNITS($2,'CHANNEL'))", rbs_vals, range(1,513))
				#my_node.putData(sig_expr)
				raw = Int32Array(rbs_vals)
				raw = raw.setUnits('Counts')
				dim = Int32Array(range(1,513))
				dim = dim.setUnits('Channel')
				sig = Signal('$VALUE', raw, dim)
				my_node.putData(sig)
				

				# W Counts data.
				wCount = sheet[row + '515'].value
				path = '\\DP_PROBES::TOP.' + letter_probe + '.' + letter_probe + u_or_d + '.RBS.RUN' + count_str + ':' + 'w_counts'
				my_node = tree.getNode(path)
				wCount = Int32(wCount)
				wCount = wCount.setUnits('Counts')
				my_node.putData(wCount)	

				# Microcoulomb data.
				microcol = sheet[row + '516'].value
				path = '\\DP_PROBES::TOP.' + letter_probe + '.' + letter_probe + u_or_d + '.RBS.RUN' + count_str + ':' + 'microcol'
				my_node = tree.getNode(path)
				my_node.putData(microcol)
	
				# W Areal Density
				w_areal = sheet[row + '517'].value
				w_areal_error = sheet[row + '518'].value
				path = '\\DP_PROBES::TOP.' + letter_probe + '.' + letter_probe + u_or_d + '.RBS.RUN' + count_str + ':' + 'w_areal'
				my_node = tree.getNode(path)
				w_areal = Float64(w_areal)
				w_areal = w_areal.setUnits('W/cm^2')
				w_areal_error = Float64(w_areal_error)
				w_areal - w_areal.setError(w_areal_error)
				#expr = Data.compile("BUILD_WITH_UNITS(BUILD_WITH_ERROR($1, $2), 'W/cm^2')", w_areal, w_areal_error)
				my_node.putData(w_areal)

				# Location
				loc = sheet[row + '525'].value
				path = '\\DP_PROBES::TOP.' + letter_probe + '.' + letter_probe + u_or_d + '.RBS.RUN' + count_str + ':' + 'loc'	
				my_node = tree.getNode(path)
				loc = Int32(loc)
				loc = loc.setUnits('mm')
				my_node.putData(loc)



def get_ICPMS(tree, letter_probes, shot):
		
	# Ask user which probe data is being inserted for.
	another = True
	while (another == True):
		while (True):
			print "Which probe is ICP-MS data being added for? Please select from the following: \nAD, AU, BD, BU, CD, CU"
			print "Enter 'q' to quit."
			probe = raw_input("--> ")
			if (probe == 'q'): break
			elif probe not in ['AD', 'AU', 'BD', 'BU', 'CD', 'CU']:
				print "Error: Incorrect entry. Please try again."
			else: break

		# Get the location of the ICPMS measurements for the samples.
		if (probe == 'q'): break

		locations = input("Enter in measured locations, separated by commas: ")
		concentrations = input("Enter in concentrations used for this probe, separated by commas: ")

		# Get the .scn files for each ppb at each location.
		conc_files_all = []
		for location in locations:
			conc_files = []
			for conc in concentrations:
				print "Select .scn file for " + str(conc) + " ppb at " + str(location) + " mm..."
				Tk().withdraw()
				filename = askopenfilename()
				conc_files.append(filename)
			conc_files_all.append(conc_files)		
		
		# Get the standard used for this probe.
		standards = []
		print "Select the five standard .scn files used."
		for value in range(1,6):
			print "Standard " + str(value) + "..."
			standards.append(askopenfilename())




		# Start filling in the tree. Starting with the locations.
		for number in range(1, len(locations)+1):
			print "Adding data for location " + str(location[number-1])
			path = '\\DP_PROBES::TOP.' + probe[0] + '.' + probe + '.ICPMS.LOC' + str(number) + ':POSITION'
			my_node = tree.getNode(path)
			my_node.putData(locations[number-1])

			# Then fill in concentration values.
			for sample in range(1, len(concentrations)+1):
				path = '\\DP_PROBES::TOP.' + probe[0] + '.' + probe + '.ICPMS.LOC' + str(number) + '.SPECTRUM' + str(sample) + ':CONC'
				my_node = tree.getNode(path)
				my_node.putData(concentrations[sample-1])

			# Then the .scn files.
			for m in conc_files_all:
				for n in m:
					print "Adding file: " + str(n)
					with open(n) as f:
						content = f.readlines()
					content = [x.strip() for x in content]
					counts = [float(x) for x in content[4:len(content)-2]]
					path = '\\DP_PROBES::TOP.' + probe[0] + '.' + probe + '.ICPMS.LOC' + str(number) + '.SPECTRUM' + str(sample) + ':DATA'
					my_node = tree.getNode(path)
					sig_expr = Data.compile("BUILD_SIGNAL($VALUE, BUILD_WITH_UNITS($1,'COUNTS'), \
								 BUILD_WITH_UNITS($2,'CHANNEL'))", counts, range(0,len(counts)))		
					my_node.putData(sig_expr)
			# Then the standard .scn files.
			count = 0
			for m in standards:
				count = count + 1
				print "Adding standard: " + str(m)
				with open(m) as f:
					content = f.readlines()
				content = [x.strip() for x in content]
				counts = [float(x) for x in content[4:len(content)-2]]
				path = '\\DP_PROBES::TOP.' + probe[0] + '.' + probe + '.ICPMS.LOC' + str(number) + '.STANDARDS.STANDARD' + str(count) + ':DATA'
				my_node = tree.getNode(path)
				sig_expr = Data.compile("BUILD_SIGNAL($VALUE, BUILD_WITH_UNITS($1,'COUNTS'), \
							 BUILD_WITH_UNITS($2,'CHANNEL'))", counts, range(0,len(counts)))		
				my_node.putData(sig_expr)


			print ""

		 	

		


		# Ask if user wants to select data for another probe.
		print "Insert data for another probe (y/n)?"
		answer = None
		while (answer not in ['y', 'n']):
			answer = raw_input("--> ")
			if (answer == 'y'):
				another = True
				break
			elif (answer == 'n'):
				another = False
				break
			else:
				print "Please answer (y/n)."


		
