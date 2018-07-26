# Collector-Probes

Scripts related to accessing and analyzing collector probe data from the 2016 Metal Rings Campaign. What follows is an example session  to pull the data into an Excel file, written in python3. Note that this reqiures the following packages: numpy, pandas, openpyxl, scipy and MDSplus.

```
In [1]: import get_cp_data as cp
```
A message will appear. Consider looking at help(cp) as well.
```
In [2]: probes = ['A35']
In [3]: probe_list = cp.get_cp_data(probes, time_start=2500, time_end=5000, time_step=500, remote=True)
SSH link r2d2 to localhost. Press enter to continue...
```
"remote" is a boolean asking if you are running this on your own machine or on the d3d network (i.e. iris). If it is set to True, you must ssh link your localhost to atlas or r2d2 when the script prompts you to (contact Andreas Wingen for r2d2 access). To link them to localhost, open a new terminal and type:
```
$ ssh -Y -p 2039 -L 8000:r2d2.gat.com:8000 username@cybele.gat.com
```
Where you can swap out r2d2 with atlas when needed, and username is your Cybele username. You will prompted to login to Cybele. Afterwards return to the python session in your initial terminal and press enter. The data stored on r2d2 will be loaded:
```
Connecting to r2d2... Connected.
Retrieving shots...
Retrieving rprobe...
Shots to be analyzed: [167536 167537]

Loading AU35 data...
Retrieving RBS locations...
...
Data from r2d2 pulled.

SSH link atlas to localhost. Press enter to continue...
```
Now go back to the ssh link terminal and logout, then ssh link again except this time swap r2d2 out with atlas. Return to the python terminal and press enter:
```
Connecting to atlas... Connected.
Analyzing atlas relevant data...

Loading gfile:
  Shot: 167536
  Tree: EFIT01
  Time: 3000
...
Save to Excel files (y/n)?
```
Enter 'y' to save to an Excel file. 

Listing of the relevant files here:

- **get_cp_data.py** - The master script for accessing the collector probe data.

- **create_df.py** - Used in get_cp_data.py. The meat of the analysis involving mapping the locations to fields lines and such.

- **meas_locations.py** - Used in get_cp_data.py. Accounts for the fact that the arm the probes are inserted in on is at a 13 degree angle, so this is just a bunch of trigonometry to put a location on a probe face into machine coordinates.

- **pull_mdsplus.py** - Used in get_cp_data.py. Functions related to pulling the RBS and LAMS data from the MDSplus server on r2d2.

- **get_ts.py** - Script to pull data from Thomson Scattering (core, divertor or tangential) into a python dictionary. In addition to pulling the data, it maps it to flux surfaces to get R-Rsep omp values. Check the function in it called "run_script" for instructions. Note this also requires an ssh link to atlas (see above). This script also includes a general purpose function for loading a gfile, used in mapping to flux surfaces and such.

- **get_lp.py** - Script to get the Langmuir probe data into a python dictionary. Essentially just a python translation of a script written by Jon Watkins with minor changes.

- **get_conn_length.py** - Script to get the poloidal connection length of a set of R values for a probe and plot it. Get the R values along a probe face from get_cp_data.py. Loads the gfiles using the function in get_ts.py.

- **get_avgZmag_and_rsep.py** - An old python2 script to get the average location of the magnetic axis Z coordinate and the average R value of the separatrix at the three Z locations of the probes. 

------------------------------------------------------------------------------------------------------------------
# Cloning the Repositories on Iris

Cloning the repositories to Iris is not as straightforward as cloning to your local machine. It requires adding a new
ssh key to your github account. The instructions located here worked for me:

[Adding new SSH Key Instructions](https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/)

where Step 1 (the one with xclip) is being run on Iris. You then copy that key to your account. Then in order to clone
the Collector-Probe repository, you run `git clone [Clone With SSH Link]` on Iris. 

------------------------------------------------------------------------------------------------------------------
  
Email any questions to Shawn Zamperini at zamp@utk.edu.
