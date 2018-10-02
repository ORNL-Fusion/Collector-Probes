# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 16:59:31 2018

@author: jduran2
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

# Choose directory and file name
os.chdir(r'C:\Users\jduran2\Dropbox (ORNL)\UTK\SCGSR\LAMS Data\Probe Maps')

file = 'AU09_Map2_dict.xlsx'


# Read data from excel file using pandas
data = pd.read_excel(file,
sheetname=1,
header=0,
index_col=False,
keep_default_na=True
)


# Specify value to locate in Radial data column and pull associated columns in indexed rows
#findvalue = 7.1119
findvalue = data['Radial [mm]'][36]

index_vals = data['Radial [mm]'][data['Radial [mm]']==findvalue].index.values
                 
df = data.iloc[index_vals]
df['NormTotal'] = df.loc[:,'total']/df['total'].max()

# Plot desired data and control font sizes
fig, ax = plt.subplots()

#ax.errorbar(df['Poloidal [mm]'], df['NormTotal'], yerr=None,fmt='.', c='b', label = 'LAMS') #Use for D side of Probe
ax.errorbar(df['z [mm]'], df['NormTotal'], yerr=None,fmt='.', c='b', label = 'LAMS')#Use for U side of Probe

ax.set_ylabel('Normalized W Intensity')
ax.set_xlabel('Z Location [mm]')

font_size = 15

for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
              ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(font_size)