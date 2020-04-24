"""
This script is meant to be copy/pasted into an OMFIT session, so that it can
pull out the "sliced" Thomson data generated from OMFITprofiles. Sliced means
data during an ELM has been excluded. Steps to use are as such:

1. Open up an OMFIT session on iris. This normally involves at the command line
    typing "module load omfit" and then "omfit".
2. Import OMFITprofiles. This is done via File -> Import module... -> OMFITprofiles.
    You should now see it on the home screen.
3. Double-click OMFITprofiles to bring up the GUI. Insert the shot and times you
    want. For times, you can use, for example, range(2000, 4510, 10) to get
    all the times between 2000-4500 in 10 ms increments.
4. I could try and explain the purpose of each tab, but you will probably be fine
    just clicking "Fetch + Slice + Fit" under Main. This uses the default values
    which probably are good enough. This could take a long time if you have alot
    of times.
5. Minimize the GUI and go back to the OMFIT home screen. In the bottom right
    you should see a blank white command box. Copy and paste this script into
    it, and press Execute.
6. Go back to your normal iris terminal. In your home directory should now be
    a file called "omfit_excel_file.xlsx". This is the file used in create_ts_omfit
    in oedge_plots. Copy this file over to your computer and use it as the input
    "omfit_path" in create_ts_omfit. WinSCP is nice for copying over files.
7. Pat yourself on the back for being a good person (good people only).
"""

import pandas as pd
import numpy  as np
import matplotlib.pyplot as plt
from os.path import expanduser


def uncert_var(arr):
    val = np.array([a.n for a in arr])
    err = np.array([a.s for a in arr])
    return val, err

# Times entered in. Also have times tiled 54 times, one set for each of the 54 channels.
time = OMFIT['OMFITprofiles']['OUTPUTS']['SLICE']['TS']['time'].data
time_out = np.tile(time, 54)

# Data that needs to be repeated.
subsystem = np.repeat(OMFIT['OMFITprofiles']['OUTPUTS']['SLICE']['TS']['subsystem'].data, len(time))
channel   = np.repeat(OMFIT['OMFITprofiles']['OUTPUTS']['SLICE']['TS']['channel'].data, len(time))

# Data that just needs to be flattened.
psin = OMFIT['OMFITprofiles']['OUTPUTS']['SLICE']['TS']['psi_n'].data.flatten()
r    = OMFIT['OMFITprofiles']['OUTPUTS']['SLICE']['TS']['R'].data.flatten()
z    = OMFIT['OMFITprofiles']['OUTPUTS']['SLICE']['TS']['Z'].data.flatten()

# Data that needs to be flattened but has error data in it
te, te_err = uncert_var(OMFIT['OMFITprofiles']['OUTPUTS']['SLICE']['TS']['T_e'].data.flatten())
ne, ne_err = uncert_var(OMFIT['OMFITprofiles']['OUTPUTS']['SLICE']['TS']['n_e'].data.flatten())

# Bundle everything up into a dataframe, save to Excel file.
all_data = np.vstack((subsystem, channel, time_out, psin, r, z, te, te_err, ne, ne_err))
df = pd.DataFrame(all_data.T, columns=('subsystem', 'channel', 'time', 'psin', 'r', 'z', 'te', 'te_err', 'ne', 'ne_err'))
home = expanduser("~")
df.to_excel(home + '/omfit_excel_file.xlsx')

# Plot to make sure it looks reasonable.
x = df['time'].unique()
y = df[df['channel']=='TS_divertor_r-1_4']['te'].values
fig, ax = plt.subplots()
ax.plot(x, y)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Te (eV)')
fig.tight_layout()
fig.show()
