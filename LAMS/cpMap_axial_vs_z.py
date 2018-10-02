# -*- coding: utf-8 -*-
"""
Created on Wed May  9 12:15:54 2018

@author: jduran2
"""
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import openpyxl as xl
import scipy
from scipy.interpolate import interpn
import os
from matplotlib.ticker import FormatStrFormatter

    
def putIntoArray(sheet, minRange, maxRange):
    cells = sheet[minRange:maxRange]
    cells = np.array(cells)
    cells = np.reshape(cells, cells.size)
    values = [cell.value for cell in cells]
    values = np.transpose(values)
    return values    

os.chdir(r'C:\Users\jduran2\Dropbox (ORNL)\UTK\SCGSR\LAMS Data\Probe Maps') # pick the working directory


file_name = "CD06_Map.xlsx"
data_min = 2
data_max = 2223

wb = xl.load_workbook(file_name, data_only=True) #change this to a user input string for which file to analyze
worksheets = wb.get_sheet_names()
sheet = wb.get_sheet_by_name(worksheets[1])


#for U
data = putIntoArray(sheet, "C"+ str(data_min), "C" + str(data_max)) # Typically D for *dict and D for raw files
rmm = putIntoArray(sheet, "A"+ str(data_min), "A" + str(data_max))
Pmm = putIntoArray(sheet, "B"+ str(data_min), "B" + str(data_max))

#for D
#data = putIntoArray(sheet, "C"+ str(data_min), "C" + str(data_max))
#rmm = putIntoArray(sheet, "A"+ str(data_min), "A" + str(data_max))
#Pmm = putIntoArray(sheet, "B"+ str(data_min), "B" + str(data_max))

rmm=np.around(rmm,decimals=1)

z=data
x=rmm
y=Pmm
z_array=[]

#z2 = scipy.interpolate.interp2d(x, y, z, kind='linear')

#Use this for 2D contour plots
"""
from scipy.interpolate import griddata

X,Y = np.meshgrid(x,y)
Z = griddata((x,y),z,(x[None,:],y[:,None]), method='nearest', rescale='false')
#Z = interpn((x,y),z,(x[None,:],y[:,None]))
CP = plt.contourf(X,Y,Z,100,cmap=cm.jet)
plt.scatter(x, y, marker='o', s=5, zorder=10)

fsize=19
plt.title('CP Contour Plot', size=fsize)
plt.xlabel('r [mm]', size=fsize)
plt.ylabel('z [mm]', size=fsize)
cbar = plt.colorbar(CP)
cbar.ax.set_ylabel('Total W counts via LAMS', size=fsize)
plt.rcParams.update({'font.size': fsize})
"""
#Use this for 2Dcolormaps of Raw Data

import pandas as pd
import seaborn as sns

df = pd.DataFrame.from_dict(np.array([x,y,z]).T)
df.columns = ['Axial Location [mm]','z Location [mm]','Z_value']
df['Z_value'] = pd.to_numeric(df['Z_value'])
pivotted= df.pivot('z Location [mm]','Axial Location [mm]','Z_value')
ax=sns.heatmap(pivotted,cmap=cm.jet,cbar_kws={'label': 'Total W Intensity'}, vmax=30000) #vmax=###### for cbar limit 'Total W Intensity'
ax.invert_yaxis()
#ax.set(xlim=(0, 201)) use if data set is for 100mm and want 50mm

for ind, label in enumerate(ax.get_xticklabels()):
    if ind % 20 == 0:  # every 10th label is kept
        label.set_visible(True)
    else:
        label.set_visible(False)
        
for ind, label in enumerate(ax.get_yticklabels()):
    if ind % 2 == 0:  # every 10th label is kept
        label.set_visible(True)
    else:
        label.set_visible(False)

font_size = 20

for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
              ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(font_size)

ax.figure.axes[-1].yaxis.label.set_size(font_size)
ax.tick_params(labelsize=font_size)
cax = plt.gcf().axes[-1]
cax.tick_params(labelsize=font_size)



#ax.figure.savefig("output.png")

#Use this for 3Dcolormaps
"""
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.2)
#ax.zaxis.set_scale('log')

for angle in range(0, 360):
    ax.view_init(90, angle)
    plt.draw()
    plt.pause(.001)

plt.show()
"""

