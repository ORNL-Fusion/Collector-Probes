# -*- coding: utf-8 -*-

"""
Created on Fri Jan 12 13:23:32 2018

@author: jduran2

Workflow should include the following
    - Find excel sheet
    - Open First sheet in excel document
    - Go to line 5 and column 1 to acquire time data
    - repeat for isotope data 
    - select a starting point for calculations
    - truncate data past 100mm
    *- calculate scan time
    - scan distance
    - select a background region and subtract counts from data
    - use an if statement to force negative counts to zero
    - use this final data version to propagate error
    - calculate enrichment fractions
    - calculate data totals to compare with RBS

"""

import numpy as np
import openpyxl as xl
import csv
import os

#Find excel sheet of interest
#Use the full directory with \\ in the file name
#write a function or script to 
#ExcelPath = 'C:\\Users\\jduran2\\Dropbox\\SCGSR\\LAMS Data\\A32\\AU32.xlsx'
#xl = pd.ExcelFile(ExcelPath)

# Take an Excel sheet and a range of cells (i.e. "H3":"H64") and return an array of the values.
def putIntoArray(sheet, minRange, maxRange):
    cells = sheet[minRange:maxRange]
    cells = np.array(cells)
    cells = np.reshape(cells, cells.size)
    values = [cell.value for cell in cells]
    values = np.transpose(values)
    return values

# Open the plunging LP workbook and grab the sheet with all the data.

os.chdir(r'C:\Users\jduran2\Dropbox (ORNL)\UTK\SCGSR\LAMS Data\Probe Maps') # pick the working directory


file_name = input('What is the name of the file to analyze i.e. AU##.xlsx : ')

bkg_min = input('Input starting row for background range: ')
bkg_max = input('Input ending row for background range: ')
data_min = input('Input beginning row of LAMS data: ') #some user input for min data range
data_max = input('Input ending row of LAMS data: ') #some user input for max data  range or calculate  based on min data  range assuming same scan rate

wb = xl.load_workbook(file_name, data_only=True) #change this to a user input string for which file to analyze
worksheets = wb.get_sheet_names()
sheet = wb.get_sheet_by_name(worksheets[0])
LAMS_Time = putIntoArray(sheet, "A"+ data_min, "A" + str(data_max))
LAMS_180 = putIntoArray(sheet, "B"+ data_min, "B" + str(data_max))
LAMS_182 = putIntoArray(sheet, "C"+ data_min, "C" + str(data_max))
LAMS_183 = putIntoArray(sheet, "D"+ data_min, "D" + str(data_max))
LAMS_184 = putIntoArray(sheet, "E"+ data_min, "E" + str(data_max))
LAMS_186 = putIntoArray(sheet, "F"+ data_min, "F" + str(data_max))
bkg_180 = putIntoArray(sheet, "B"+ bkg_min, "B" + bkg_max)
bkg_182 = putIntoArray(sheet, "C"+ bkg_min, "C" + bkg_max)
bkg_183 = putIntoArray(sheet, "D"+ bkg_min, "D" + bkg_max)
bkg_184 = putIntoArray(sheet, "E"+ bkg_min, "E" + bkg_max)
bkg_186 = putIntoArray(sheet, "F"+ bkg_min, "F" + bkg_max)

#create scan time array
ScanTime = []
ScanTime.append(0)

for row in range(1, len(LAMS_Time)-1): #start from 1, to leave out row 0
        ScanTime.append(ScanTime[row-1] + LAMS_Time[row] - LAMS_Time[row-1])
        
#create scan distance array
ScanDist = []
for row in range(0, len(ScanTime)): #start from 0
        ScanDist.append(ScanTime[row]*0.5)
        
#create average background value in gas blank

avg_bkg_sub_180 = np.mean(bkg_180)
avg_bkg_sub_182 = np.mean(bkg_182)
avg_bkg_sub_183 = np.mean(bkg_183)
avg_bkg_sub_184 = np.mean(bkg_184)
avg_bkg_sub_186 = np.mean(bkg_186)
        
#create background subtracted data
bkg_sub_180 = []
bkg_sub_182 = []
bkg_sub_183 = []
bkg_sub_184 = []
bkg_sub_186 = []
final_180 = []
final_182 = []
final_183 = []
final_184 = []
final_186 = []

for row in range(0, len(ScanTime)): #start from 0
        bkg_sub_180.append(LAMS_180[row]-avg_bkg_sub_180)
        bkg_sub_182.append(LAMS_182[row]-avg_bkg_sub_182)
        bkg_sub_183.append(LAMS_183[row]-avg_bkg_sub_183)
        bkg_sub_184.append(LAMS_184[row]-avg_bkg_sub_184)
        bkg_sub_186.append(LAMS_186[row]-avg_bkg_sub_186)

        
#create a function for generating final arrays/lists
def final_W_data(sizing_data, bkg_sub_data):
    final_data = []
    for row in range(0, len(sizing_data)): #start from 0
        if bkg_sub_data[row] <= 0:
            final_data.append(0)
        else:
            final_data.append(bkg_sub_data[row])   
    return final_data
        
final_180 = final_W_data(ScanTime, bkg_sub_180)
final_182 = final_W_data(ScanTime, bkg_sub_182)
final_183 = final_W_data(ScanTime, bkg_sub_183)        
final_184 = final_W_data(ScanTime, bkg_sub_184)
final_186 = final_W_data(ScanTime, bkg_sub_186)

#create total W signal data
total_W = []
for row in range(0, len(ScanTime)): #start from 0
        total_W.append(final_180[row] + final_182[row] + final_183[row] + final_184[row] + final_186[row])
        
#create signal error based on detector counting statistics or root(N) and error propagation
err180 = []
err182 = []
err183 = []
err184 = []
err186 = []
err_total = []

for row in range(0, len(ScanTime)): #start from 0
        err180.append(np.sqrt(final_180[row]))
        err182.append(np.sqrt(final_182[row]))
        err183.append(np.sqrt(final_183[row]))
        err184.append(np.sqrt(final_184[row]))
        err186.append(np.sqrt(final_186[row]))

for row in range(0, len(ScanTime)):
        err_total.append(np.sqrt(np.power(err180[row],2)+np.power(err182[row],2)
        +np.power(err183[row],2)+np.power(err184[row],2)+np.power(err186[row],2)))
        
#create enrichment fractions
EF180 = []
EF182 = []
EF183 = []
EF184 = []
EF186 = []

for row in range(0, len(ScanTime)): #start from 0
        EF180.append(np.divide(final_180[row],total_W[row]))
        EF182.append(np.divide(final_182[row],total_W[row]))
        EF183.append(np.divide(final_183[row],total_W[row]))
        EF184.append(np.divide(final_184[row],total_W[row]))
        EF186.append(np.divide(final_186[row],total_W[row]))
        
#create enrichment fraction error based on detector counting statistics and error propagation
EF180err = []
EF182err = []
EF183err = []
EF184err = []
EF186err = []

for row in range(0, len(ScanTime)): #start from 0
        EF180err.append(np.multiply(EF180[row],np.sqrt(np.power(np.divide(err180[row],final_180[row]),2) + np.power(np.divide(err_total[row],total_W[row]),2))))
        EF182err.append(np.multiply(EF182[row],np.sqrt(np.power(np.divide(err182[row],final_182[row]),2) + np.power(np.divide(err_total[row],total_W[row]),2))))
        EF183err.append(np.multiply(EF183[row],np.sqrt(np.power(np.divide(err183[row],final_183[row]),2) + np.power(np.divide(err_total[row],total_W[row]),2))))
        EF184err.append(np.multiply(EF184[row],np.sqrt(np.power(np.divide(err184[row],final_184[row]),2) + np.power(np.divide(err_total[row],total_W[row]),2))))
        EF186err.append(np.multiply(EF186[row],np.sqrt(np.power(np.divide(err186[row],final_186[row]),2) + np.power(np.divide(err_total[row],total_W[row]),2))))

#create a function for removing nan values
def final_W_err_convert(sizing_data, EF_err):
    final_EF_err = []
    for row in range(0, len(sizing_data)): #start from 0
        if EF_err[row] == 'nan':
            final_EF_err.append(0)
        else:
            final_EF_err.append(EF_err[row])   
    return final_EF_err
        
EF180err = final_W_err_convert(ScanTime, EF180err)
EF182err = final_W_err_convert(ScanTime, EF182err)
EF183err = final_W_err_convert(ScanTime, EF183err)        
EF184err = final_W_err_convert(ScanTime, EF184err)
EF186err = final_W_err_convert(ScanTime, EF186err)
       
W_Dict = {'ScanTime': ScanTime, 'ScanDist': ScanDist, 'W180': final_180, '180err': err180, 
    'W182': final_182, '182err': err182, 'W183': final_183, '183err': err183, 'W184': final_184,
    '184err': err184, 'W186': final_186, '186err': err186, 'Total W': total_W, 'Total W err': err_total,
    'EF180':EF180, 'EF180err':EF180err,'EF182':EF182, 'EF182err':EF182err,'EF183':EF183, 'EF183err':EF183err,
    'EF184':EF184, 'EF184err':EF184err,'EF186':EF186, 'EF186err':EF186err}

newstr = file_name.replace(".xlsx", "_")

with open(newstr +'dict.csv', 'w', newline='') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(W_Dict.keys())
    writer.writerows(zip(*W_Dict.values()))
    

###############################################################################
#Plotting
###############################################################################

#import matplotlib.pyplot as plt
#
#
#W180merr = [i - j for i, j in zip(EF180, EF180err)]
#W180perr = [i + j for i, j in zip(EF180, EF180err)]
#W182merr = [i - j for i, j in zip(EF182, EF182err)]
#W182perr = [i + j for i, j in zip(EF182, EF182err)]
#W183merr = [i - j for i, j in zip(EF183, EF183err)]
#W183perr = [i + j for i, j in zip(EF183, EF183err)]
#W184merr = [i - j for i, j in zip(EF184, EF184err)]
#W184perr = [i + j for i, j in zip(EF184, EF184err)]
#W186merr = [i - j for i, j in zip(EF186, EF186err)]
#W186perr = [i + j for i, j in zip(EF186, EF186err)]
#
###############################################
#
##Create the plot
#fig, ax1 = plt.subplots() 
#
##Plot axis 1 data
##ax1.scatter(W_RBS_x, W_RBS, s=10, c='r', marker="^", label='RBS')
##ax1.errorbar(W_RBS_x,W_RBS, yerr=W_RBS_err, c='r', zorder=60)
##ax1.fill_between(W_RBS_x, RBS_merr, RBS_perr, facecolor='r', interpolate=True)
#ax1.set_xlim(xmin=0,xmax=100)
#ax1.set_ylabel('W Concentration [1e15 W/cm^2]', size=14, color='r')
#ax1.set_xlabel('Distance Along Collector Probe [mm]', size=14)
#ax1.tick_params('y', colors='r', size=14)
#
#
##Plot axis 2 data
#ax2 = ax1.twinx()
##ax2.scatter(scanDist, W180R, s=10, c='b', marker="o", label='W180')
##ax2.scatter(scanDist, W182R, s=10, c='g', marker="o", label='W182')
##ax2.scatter(scanDist, W183R, s=10, c='r', marker="o", label='W183')
##ax2.scatter(scanDist, W184R, s=10, c='c', marker="o", label='W184')
##ax2.scatter(scanDist, W186R, s=10, c='m', marker="o", label='W186')
#ax2.plot(ScanDist, EF180, c='b', label='W180')
#ax2.fill_between(ScanDist, W180merr, W180perr, facecolor='b', interpolate=True, zorder=10)
#ax2.plot(ScanDist, EF183, c='y', label='W183')
#ax2.fill_between(ScanDist, W183merr, W183perr, facecolor='y', interpolate=True, zorder=40)
#ax2.plot(ScanDist, EF184, c='c', label='W184')
#ax2.fill_between(ScanDist, W184merr, W184perr, facecolor='c', interpolate=True, zorder=30)
#ax2.plot(ScanDist, EF186, c='m', label='W186')
#ax2.fill_between(ScanDist, W186merr, W186perr, facecolor='m', interpolate=True, zorder=20)
#ax2.plot(ScanDist, EF182, c='g', label='W182')
#ax2.fill_between(ScanDist, W182merr, W182perr, facecolor='g', interpolate=True, zorder=50) #plot this last because zorder doesn't work
#ax2.set_xlim(xmin=0,xmax=100)
#ax2.set_ylim(ymin=.0001,ymax=1)
#ax2.set_xlim(xmin=0,xmax=100)
#ax2.set_ylabel('Enrichment Fraction',size=14)
#ax2.set_yscale('log') #set the scale to logarithmic
#ax2.tick_params('y', size=14)
#ax2.legend(bbox_to_anchor=(1.2, 1), loc=2, borderaxespad=0.)
#
##Add a second x-axis below the first
##Define your axis limits based on the trendline of RBSx vs. R-Rsepx
#
#RminRsep_min = float(input('Input the RminRsep axis minimum value: '))
#RminRsep_max = float(input('Input the RminRsep axis maximum value: '))
##Add space to the bottom
#fig.subplots_adjust(bottom=0.2)
##Define axis
#ax3 = ax1.twiny()
#ax3.yaxis.set_visible(False) # hide the yaxis
## Move twinned axis ticks and label from top to bottom
#ax3.xaxis.set_ticks_position("bottom")
#ax3.xaxis.set_label_position("bottom")
## Offset the twin axis below the host
#ax3.spines["bottom"].set_position(("axes", -0.2))
#ax3.set_xlim(RminRsep_min,RminRsep_max)
#ax3.set_xlabel('R-Rsep [mm]', size=14)
#
##Plot axis 1 data AGAIN for plotting order....
##ax1.scatter(W_RBS_x, W_RBS, s=10, c='r', marker="^", label='RBS')
#ax5 = ax1.twiny()
##ax5.errorbar(W_RBS_x,W_RBS, yerr=W_RBS_err, c='r', linewidth=3, zorder=60)
##ax1.fill_between(W_RBS_x, RBS_merr, RBS_perr, facecolor='r', interpolate=True)
#ax5.set_xlim(xmin=0,xmax=100)
#ax5.set_ylabel('W Concentration [1e15 W/cm^2]', size=14, color='r')
#ax5.set_xlabel('Distance Along Collector Probe [mm]', size=14)
#ax5.tick_params('y', colors='r', size=14)
#ax5.xaxis.set_visible(False) # hide the xaxis
#
#
##Add top ticks
#ax4 = ax1.twiny()
#ax4.set_xticklabels([])
#
#plt.show()
#
#        
