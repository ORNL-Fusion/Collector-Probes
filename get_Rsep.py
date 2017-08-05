# In separate terminal use 
# ssh -Y -p 2039 -L 8000:atlas.gat.com:8000 username@cybele.gat.com
# where you substitute in your username. The MDSplus connection is
# then MDSplus.Connection('localhost').
import numpy as np
import scipy.interpolate as scinter

import EFIT.load_gfile_d3d_sz as loadg
import Geometry.meas_locations_v3_3 as geo

# Return dictionary of average R - Rsep value for each of the probes.
#
# shots     = list of shots in for. If only one shot still enter as list.
# r_probe   = radial position of probe holder tip (cm).
# location  = location along probe where you want the R - Rsep value (cm)
# startTime = start of time range to be averaged over (ms).
# endTime   = end of time range (ms).
#
def avg_Rsep(shots, r_probe, location, writeToFile=False, filename="IDidNotEnterAFilename.txt", startTime=2500, endTime=5000):
	time = startTime
	count = 0
	
	sums = {'ad':0, 'au':0, 'bd':0, 'bu':0, 'cd':0, 'cu':0}
	
	# Radial position of each probe.
	rad_pos = {
		'ad': geo.calc_R_measAD(r_probe, location),
		'au': geo.calc_R_measAU(r_probe, location),
		'bd': geo.calc_R_measBD(r_probe, location),
		'bu': geo.calc_R_measBU(r_probe, location),
		'cd': geo.calc_R_measCD(r_probe, location),
		'cu': geo.calc_R_measCU(r_probe, location)}
	
	for shot in shots:	
		time = startTime
		while time <= endTime:
			
			print("Shot:     " + str(shot))
			print("Location: " + str(location))
			print("Time:     " + str(time))
			
			# Lines from Zeke's code.
			parmDICT = loadg.read_g_file_mds(shot, time, Server='localhost', write2file=False)
			Rs, Zs = np.meshgrid(parmDICT['R'], parmDICT['Z'])
			Zes = np.copy(parmDICT['lcfs'][:, 1][13:-12])
			Res = np.copy(parmDICT['lcfs'][:, 0][13:-12])
			f_Rs = scinter.interp1d(Zes, Res, assume_sorted=False)
			
			# R_Sep for each z location of the three probes.
			rSep = {}
			rSep['a'] = f_Rs(-0.18)
			rSep['b'] = f_Rs(-0.1546)
			rSep['c'] = f_Rs(-0.2054)
			
			# rSep is in meters. Convert to cm.
			for value in rSep:
				rSep[value] = rSep[value] * 100.0
			
			# Sums for calculating avg R - Rsep.
			sums['ad'] += rad_pos['ad'] - rSep['a']
			sums['au'] += rad_pos['au'] - rSep['a']
			sums['bd'] += rad_pos['bd'] - rSep['b']
			sums['bu'] += rad_pos['bu'] - rSep['b']
			sums['cd'] += rad_pos['cd'] - rSep['c']
			sums['cu'] += rad_pos['cu'] - rSep['c']
			
			# Next time step.
			time += 500
			count += 1
			print "\n"
			
	# Calculate avg R - Rsep
	avg_RminRseps = {
		'location':location,
		'ad':float(sums['ad']) / count,
		'au':float(sums['au']) / count,
		'bd':float(sums['bd']) / count,
		'bu':float(sums['bu']) / count,
		'cd':float(sums['cd']) / count,
		'cu':float(sums['cu']) / count}
		
	#print avg_RminRseps['ad']
		
	# Save to a txt file.
	if writeToFile:
		#filename = raw_input('Enter filename: ')
		f = open(filename, 'a')
		f.write("Shot: " + str(shot) + "\n")
		for key, value in avg_RminRseps.iteritems():
			f.write(key + " ")
			f.write(str(value))
			f.write("\n")
		f.write("\n")
		f.close()
		
		
	return avg_RminRseps
		
		
		
		
			
			
