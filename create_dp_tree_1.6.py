from MDSplus import *
import sys
import openpyxl as xl
import get_data_modules_3 as mod


# Open up the model tree.
model = Tree('dp_probes', -1)

# Create new shot corresponding to the probe number, then assign it.
shot = int(sys.argv[1])

try:
	tree = Tree('dp_probes', shot)
except:
	model.createPulse(shot)
	tree = Tree('dp_probes', shot)

print "Assigning data to \'" + str(shot) + "\' probes."

# Fill in what shots each probe was in for.
letter_probes = raw_input("Enter which probes (separated by commas): ")
print letter_probes
for letter_probe in letter_probes:
	start_shot = input("\nStart shot for " + letter_probe + " probe: ")
	end_shot = input("End shot for " + letter_probe + " probe: ")
	shots_in_for = []
	for shot_tmp in range(start_shot, end_shot+1):
		shots_in_for.append(int(shot_tmp))

	path = '\\DP_PROBES::TOP.' + letter_probe + ':SHOTS'
	node = tree.getNode(path)
	node.putData(Uint32Array(shots_in_for))

mod.get_RBS(tree, letter_probes, shot)

#mod.get_ICPMS(tree, letter_probes, shot)
