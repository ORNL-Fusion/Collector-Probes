# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 10:04:19 2018

@author: J. Duran and S. Zamp

--Collector Probe W Analysis and Mapping--

"""

# set imports
import os
import pandas as pd
import numpy as np

#import matplotlib.pyplot as plt

# define file directory and file name
os.chdir(r'C:\Users\jduran2\Documents\Python Scripts\Python Test Directory')
file_name = 'AU27_Map_100um.xlsx'

####################################################
# User Inputs
bkg_min = 0
bkg_max = 80
scan_length = 50 # mm
scan_speed_um = 500 # um/sec

####################################################

# import data from excel file into pandas dataframe
df = pd.read_excel(file_name,
sheetname = 1,
header = 0,
skiprows = 3,
index_col = False,
keep_default_na = True
)

# calculate total W signature and add to df
df['Total W'] = df['W180']+df['W182']+df['W183']+df['W184']+df['W186']


# set background range of rows, average background, and subtract from remaining data
#bkg_min = input('Input starting row for background range: ') # user input option
#bkg_max = input('Input ending row for background range: ') # user input option
bkg_min = bkg_min
bkg_max = 1 + bkg_max
bkg_df = df.iloc[bkg_min:bkg_max]
mean_bkg_180 = bkg_df["W180"].mean()
mean_bkg_182 = bkg_df["W182"].mean()
mean_bkg_183 = bkg_df["W183"].mean()
mean_bkg_184 = bkg_df["W184"].mean()
mean_bkg_186 = bkg_df["W186"].mean()
mean_bkg_TotalW = bkg_df["Total W"].mean()
df['W180'] = df['W180'] - mean_bkg_180
df['W182'] = df['W182'] - mean_bkg_182
df['W183'] = df['W183'] - mean_bkg_183
df['W184'] = df['W184'] - mean_bkg_184
df['W186'] = df['W186'] - mean_bkg_186

  
# re-calculate total W signature and add to df
df['Total W'] = df['W180']+df['W182']+df['W183']+df['W184']+df['W186']


# calculate signal error based on detector counting statistics or root(N) and error propagation
for row in range(0, len(df['Time [Sec]'])): #start from 0
        df.loc[row,'W180 error'] = np.sqrt(df.loc[row,'W180'])
        df.loc[row,'W182 error'] = np.sqrt(df.loc[row,'W182'])
        df.loc[row,'W183 error'] = np.sqrt(df.loc[row,'W183'])
        df.loc[row,'W184 error'] = np.sqrt(df.loc[row,'W184'])
        df.loc[row,'W186 error'] = np.sqrt(df.loc[row,'W186'])
        df.loc[row,'Total W error'] = np.sqrt(np.power(df.loc[row,'W180 error'],2)+np.power(df.loc[row,'W182 error'],2)
+np.power(df.loc[row,'W183 error'],2)+np.power(df.loc[row,'W184 error'],2)+np.power(df.loc[row,'W186 error'],2))
        df.loc[row,'W180 EF'] = np.divide(df.loc[row,'W180'],df.loc[row,'Total W'])
        df.loc[row,'W182 EF'] = np.divide(df.loc[row,'W182'],df.loc[row,'Total W'])
        df.loc[row,'W183 EF'] = np.divide(df.loc[row,'W183'],df.loc[row,'Total W'])
        df.loc[row,'W184 EF'] = np.divide(df.loc[row,'W184'],df.loc[row,'Total W'])
        df.loc[row,'W186 EF'] = np.divide(df.loc[row,'W186'],df.loc[row,'Total W'])
        df.loc[row,'W180 EF error'] = np.multiply(df.loc[row,'W180 EF'],np.sqrt(np.power(np.divide(df.loc[row,'W180 error'],df.loc[row,'W180']),2) + np.power(np.divide(df.loc[row,'Total W error'],df.loc[row,'Total W']),2)))
        df.loc[row,'W182 EF error'] = np.multiply(df.loc[row,'W182 EF'],np.sqrt(np.power(np.divide(df.loc[row,'W182 error'],df.loc[row,'W182']),2) + np.power(np.divide(df.loc[row,'Total W error'],df.loc[row,'Total W']),2)))
        df.loc[row,'W183 EF error'] = np.multiply(df.loc[row,'W183 EF'],np.sqrt(np.power(np.divide(df.loc[row,'W183 error'],df.loc[row,'W183']),2) + np.power(np.divide(df.loc[row,'Total W error'],df.loc[row,'Total W']),2)))
        df.loc[row,'W184 EF error'] = np.multiply(df.loc[row,'W184 EF'],np.sqrt(np.power(np.divide(df.loc[row,'W184 error'],df.loc[row,'W184']),2) + np.power(np.divide(df.loc[row,'Total W error'],df.loc[row,'Total W']),2)))
        df.loc[row,'W186 EF error'] = np.multiply(df.loc[row,'W186 EF'],np.sqrt(np.power(np.divide(df.loc[row,'W186 error'],df.loc[row,'W186']),2) + np.power(np.divide(df.loc[row,'Total W error'],df.loc[row,'Total W']),2)))


# create a function for removing nan values
#def final_W_err_convert(sizing_data, EF_err):
#    final_EF_err = []
#    for row in range(0, len(sizing_data)): #start from 0
#        if EF_err[row] == 'nan':
#            final_EF_err.append(0)
#        else:
#            final_EF_err.append(EF_err[row])   
#    return final_EF_err
#        
#EF180err = final_W_err_convert(ScanTime, EF180err)
#EF182err = final_W_err_convert(ScanTime, EF182err)
#EF183err = final_W_err_convert(ScanTime, EF183err)        
#EF184err = final_W_err_convert(ScanTime, EF184err)
#EF186err = final_W_err_convert(ScanTime, EF186err)


# set scan length and scan speed
#scan_length = 50 # mm
axial_scan_data_length = scan_length*4
#scan_speed_um = 500 # um/sec
scan_speed_mm = scan_speed_um/1000 #mm/sec


# create scan time np.array
ScanTime = np.array([])
ScanTime = np.append(ScanTime,0)

for row in range(1, len(df['Time [Sec]'])-1): #start from 1, to leave out row 0
        ScanTime = np.append(ScanTime,ScanTime[row-1] + df['Time [Sec]'][row] - df['Time [Sec]'][row-1])

# create scan distance np.array
ScanDist = np.array([])
for row in range(0, len(ScanTime)): #start from 0
        ScanDist = np.append(ScanDist,ScanTime[row]*scan_speed_mm)
        
Axial_Location = np.array(ScanDist[0:axial_scan_data_length])


# filter data based on 'Total W' column
map_df = pd.DataFrame()
Axial_Location_Final = np.array([])
bkg_mult = 10 # may modify to fit each dataset
i=0
j=0 
k=0
z_Location = np.array([])

while i < len(df) and i<len(df)-5:
    if all(val> bkg_mult*mean_bkg_TotalW for val in df.iloc[[i,i+1,i+2,i+3,i+4,i+5]]['Total W'].values):
        map_df = map_df.append(df.iloc[i+1:i+axial_scan_data_length+1])
        i = i + 201
        Axial_Location_Final = np.append(Axial_Location_Final, Axial_Location)
        z_Location = np.append(z_Location,j)
        j = j + 0.1 # z Location separation distance
        k = k + 1
       # print('YES')
    else:
        i = i + 1

# add z Location [mm] column to dataframe
z_Location_Final = np.array([])
for i in range(0,k):
    z_Location_Final = np.append(z_Location_Final,np.repeat(z_Location[i],200))


# add an if statement to flip the z data with arr[::-1] when U side of probe vs normal when z
if file_name[1] == 'D':
    z_Location_Final = z_Location_Final
elif file_name[1] =='U':
    z_Location_Final = z_Location_Final[::-1]    
else:
    print('Input file_name in non-standard format')
                                                 

map_df.insert(loc=0, column='z Location [mm]', value=z_Location_Final)

# add Axial Location [mm] column to front of dataframe 'map_df'
map_df.insert(loc=0, column='Axial Location [mm]', value=Axial_Location_Final)

newstr = file_name.replace(".xlsx", "_MapData.xlsx")
    
# Export to excel file....
writer = pd.ExcelWriter(newstr, engine='xlsxwriter')
map_df.to_excel(writer, sheet_name='MapData')
writer.save()
# Plotting.... separate scripts unique to each user