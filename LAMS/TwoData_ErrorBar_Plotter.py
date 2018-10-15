# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 16:26:37 2018

@author: jduran2
"""

import numpy as np
import openpyxl as xl
import os
import matplotlib.pyplot as plt

def putIntoArray(sheet, minRange, maxRange):
    cells = sheet[minRange:maxRange]
    cells = np.array(cells)
    cells = np.reshape(cells, cells.size)
    values = [cell.value for cell in cells]
    values = np.transpose(values)
    return values

os.chdir(r'C:\Users\jduran2\Dropbox (ORNL)\UTK\SCGSR\LAMS Data\CP_LAMS_Dict') # pick the working directory
    
file_name = 'AU32_dict.xlsx'
#file_name = input('What is the name of the file to analyze i.e. AU##.xlsx : ')
        
wb = xl.load_workbook(file_name, data_only=True) #change this to a user input string for which file to analyze
worksheets = wb.get_sheet_names()
sheet = wb.get_sheet_by_name(worksheets[1])

#RBS_Dist = putIntoArray(sheet, "A2", "A11")
#RBS_Data = putIntoArray(sheet, "B2", "B11")
#RBS_Err = putIntoArray(sheet, "C2", "C11")
#LAMS_Dist = putIntoArray(sheet, "D2", "D202")
#LAMS_Data = putIntoArray(sheet, "E2", "E202")
#LAMS_Err = putIntoArray(sheet, "F2", "F202")

#model_x = putIntoArray(sheet, "C2", "C1001")
#model_shelf = putIntoArray(sheet, "D2", "D1001")
#model_floor = putIntoArray(sheet, "E2", "E1001")

#exp_x = putIntoArray(sheet, "K2", "K17")
#exp_x_err = putIntoArray(sheet, "L2", "L17")
#exp_floor = putIntoArray(sheet, "M2", "M17")
#exp_floor_err = putIntoArray(sheet, "O2", "O17")
#exp_shelf = putIntoArray(sheet, "N2", "N17")
#exp_shelf_err = putIntoArray(sheet, "P2", "P17")

scan_dist = putIntoArray(sheet, "A2", "A401")
LAMS = putIntoArray(sheet, "D2", "D401")
LAMS_err = putIntoArray(sheet, "E2", "E401")
RBS_x = putIntoArray(sheet, "G2", "G21")
RBS = putIntoArray(sheet, "J2", "J21")
RBS_err = putIntoArray(sheet, "K2", "K21")


fig, ax = plt.subplots() 

#ax.errorbar(model_x,model_shelf, label='93% W-182')
#ax.errorbar(model_x,model_floor, label="Natural W")
#ax.errorbar(exp_x, exp_floor, xerr=exp_x_err, yerr=exp_floor_err,fmt='.', c='k', label = 'Cocktail Tests')
#ax.errorbar(exp_x, exp_shelf, xerr=exp_x_err, yerr=exp_shelf_err,fmt='.', c='k')
#ax1.fill_between(W_RBS_x, RBS_merr, RBS_perr, facecolor='r', interpolate=True)

ax.errorbar(scan_dist, LAMS, yerr=LAMS_err,fmt='.', c='b', label = 'LAMS')
ax.errorbar(RBS_x, RBS, yerr=RBS_err,fmt='x', c='r', label = 'RBS')
ax.plot(RBS_x,RBS,'xr-')

#ax.set_xlim(xmin=-1,xmax=0)
#ax.set_ylim(ymin=0,ymax=1)
ax.set_ylabel('Normalized Intensity')
ax.set_xlabel('Axial Location [mm]')


#ax1.tick_params('y', colors='r', size=14)

font_size = 15

for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
              ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(font_size)
    
ax.legend(loc='center right', prop={'size': font_size})

RminRsep_min = 90.41303
RminRsep_max = 181.8124
#Add space to the bottom
fig.subplots_adjust(bottom=0.2)
#Define axis
ax2 = ax.twiny()
ax2.yaxis.set_visible(False) # hide the yaxis
# Move twinned axis ticks and label from top to bottom
ax2.xaxis.set_ticks_position("bottom")
ax2.xaxis.set_label_position("bottom")
# Offset the twin axis below the host
ax2.spines["bottom"].set_position(("axes", -0.22))
ax2.set_xlim(RminRsep_min,RminRsep_max)
ax2.set_xlabel('R-Rsep [mm]', size=15)

for item in ([ax2.title, ax2.xaxis.label, ax2.yaxis.label] +
              ax2.get_xticklabels() + ax2.get_yticklabels()):
    item.set_fontsize(font_size)