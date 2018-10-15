import matplotlib.pyplot  as plt
import pandas             as pd
import numpy              as np
from mpl_toolkits.mplot3d import Axes3D


# Location of Excel file with "TOTAL DEPOSITION" table copy/pasted into it.
filename = '/mnt/c/Users/Shawn/Documents/d3d_work/test3_totdep.xlsx'

# Load the data into a Dataframe, drop the last row that has junk in it, and
# then rename the columns that is a string '0.1' to a float 0.0 (not sure why
# this happens).
df = pd.read_excel(filename, index_col=0)
df.drop(df.columns[-1], axis=1, inplace=True)
df.rename({'0.1':0.0}, axis=1, inplace=True)

# Create the meshgrid for use in the 3D plots (X = poloidal, Y = radial?
# Z = counts). Z is multiplied by -1 because they're counts of erosion, so
# deposition shows up as a negative number.
X, Y = np.meshgrid(df.columns, df.index)
Z = df.values * -1

# Plotting commands.
def plot_3d(angle1, angle2):
    font = {'fontsize':18}
    fig = plt.figure(figsize=((15,8)))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, cmap='Reds')
    ax.set_xlabel('\nPoloidal (m)', font)
    ax.set_ylabel('\n\nRadial? Distance \nalong probe? (m)', font)
    ax.set_zlabel('\nDeposition counts', font)
    ax.view_init(angle1, angle2)
    fig.tight_layout()
    fig.show()

plot_3d(35, -85)
plot_3d(35, -32)
