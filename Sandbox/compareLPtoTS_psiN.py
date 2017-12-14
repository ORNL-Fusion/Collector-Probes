# Will need to create tunnel through Cybele to access atlas. On-site the
# normal command works. Need to test off-site.

# OMFITprofiles data in netCDF file format.
from netCDF4 import Dataset
# Plunging LP data in Excel sheet.
import openpyxl as xl
# loadgfile needed to get the psiN values.
import EFIT.load_gfile_d3d as loadg
# MDSplus needed to connect to atlas for the gfile.
import MDSplus as mds

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import curve_fit


# Some constants.
zOfLP = -0.18

# Take an Excel sheet and a range of cells (i.e. "H3":"H64") and return an array of the values.
def putIntoArray(sheet, minRange, maxRange):
    cells = sheet[minRange:maxRange]
    cells = np.array(cells)
    cells = np.reshape(cells, cells.size)
    values = [cell.value for cell in cells]
    values = np.transpose(values)
    return values

def returnTeArrAtTime(TeVar, times, time):
    index, = np.where(times==time)[0]
    TeArray = Te[:][index]
    return TeArray

# Import the script then run this function.
def runScript():
    # Open up the netCDF file root group.
    rootgrp = Dataset("OMFITprofiles_167192_FIT.nc", "r", format="NETCDF4")
    # Grab the psi_n, Te and time Variable objects.
    psi_nVar = rootgrp["psi_n"]
    TeVar    = rootgrp["T_e"]
    timesVar  = rootgrp["time"]
    # Now get the actual data from the Variable objects.
    psi_n = psi_nVar[:]
    times = timesVar[:]
    Te = TeVar[:][0]

    # Open the plunging LP workbook and grab the sheet with all the data.
    wb = xl.load_workbook("recLPdata.xlsx", data_only=True)
    sheet = wb.get_sheet_by_name("Sheet1")
    LPradii = putIntoArray(sheet, "B3", "B64")
    LPTes = putIntoArray(sheet, "D3", "D64")
    LPtimes = putIntoArray(sheet, "A3", "A64")

    # Convert to radii from cm to m.
    LPradii = [radius/100.0 for radius in LPradii]

    # Convert the LPradii at the corresponding time into a psi_N value.
    MDSplusConn = mds.Connection("atlas.gat.com")
    # Array to hold PsiN values.
    PsiNarr = np.array([])
    TStempsArr = np.array([])
    for radius in LPradii:
        # Get the index of where the requested radius is at, then get the corresponding time.
        index, = np.where(LPradii==radius)[0]
        LPtimeRequested = LPtimes[index]

        # Need to FIX. TS netCDF file only goes up to 1999. Just need to do it again on OMFIT.
        if LPtimeRequested == 2000:
            continue

        # Load the g-file.
        print "Loading g-file: (" + str(index) + "/" + str(len(LPradii) - 1) + ")"
        parameterDict = loadg.read_g_file_mds(167192, LPtimeRequested,
                                              connection=MDSplusConn,
                                              write2file=False,
                                              tree="EFIT01")
        print ""

        # Create a grid out of the R,Z values.
        Rs, Zs = np.meshgrid(parameterDict['R'], parameterDict['Z'])
        # R,Z of the magnetic axis.
        RmAxis = parameterDict["RmAxis"]
        ZmAxis = parameterDict["ZmAxis"]
        # The R's and Z's of the LCFS.
        lcfsZs = np.copy(parameterDict['lcfs'][:, 1][13:-12])
        lcfsRs = np.copy(parameterDict['lcfs'][:, 0][13:-12])

        # Function to interpolate along the LCFS. Give it a Z, and it gives you Rsep.
        fLcfs = interpolate.interp1d(lcfsZs, lcfsRs, assume_sorted=False)
        # Only grab the right half of the R's (don't care about the inner wall side).
        RsTrunc = Rs > RmAxis
        # Function to interpolate the psiN values. Give it an R and Z and it gives you psiN at that location.
        fPsiN = interpolate.Rbf(Rs[RsTrunc], Zs[RsTrunc],
                                parameterDict['psiRZn'][RsTrunc], function='linear')

        # Now we give fPsiN the radius and Z of the LP probe to the psiN's along the probe.
        tmpPsiN = fPsiN(radius, zOfLP)
        PsiNarr = np.append(PsiNarr, tmpPsiN)

    # Put it all in a dictionary and return it.
    LPdict = {}
    LPdict["LPradii"] = LPradii
    LPdict["LPtemps"] = LPTes
    LPdict["LPpsiNs"] = PsiNarr
    LPdict["TStemps"] = Te
    LPdict["TSpsiNs"] = psi_n

    return LPdict

def plotLPpsiN(LPdict=None):
    if LPdict == None:
        LPdict = runScript()

    psiNs = LPdict["LPpsiNs"]
    temps = LPdict["LPtemps"][:61]
    TS_Te = LPdict["TStemps"][:210]
    # Only use TS up to 1.05 since it stops agreeing with the LP after that.
    TS_psiN = LPdict["TSpsiNs"][:210]

    # Put all data into one array, then exponential fit to it all to fill in the gap.
    # Need to flip LP data since it's backwards.
    allPsi = np.append(TS_psiN, np.flip(psiNs,0))
    allTemp = np.append(TS_Te, np.flip(temps,0))

    def exp_fit(x, a, b, c):
        return a * np.exp(-b * x) + c

    popt, pcov = curve_fit(exp_fit, allPsi[200:], allTemp[200:])
    plt.plot(allPsi[200:], exp_fit(allPsi[200:], *popt), 'r-', label='Exp. fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))

    plt.plot(psiNs, temps, label="LP Te")
    plt.plot(TS_psiN, TS_Te, label = "TS Te")
    plt.xlim([1,1.4])
    plt.ylim([0,50])
    plt.legend()
    plt.xlabel(r"$\phi_N$")
    plt.ylabel("Temperature (eV)")
    plt.title("SOL Edge Temperatures Shot #167192")
    plt.show()

    popt, pcov = curve_fit(exp_fit, allPsi[200:], allTemp[200:])
    plt.plot(allPsi[200:], exp_fit(allPsi[200:], *popt), 'r-', label='Exp. fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))

    plt.plot(allPsi, allTemp, label="All Data")
    plt.xlim([1,1.4])
    plt.ylim([0,50])
    plt.legend()
    plt.xlabel(r"$\phi_N$")
    plt.ylabel("Temperature (eV)")
    plt.title("SOL Edge Temperatures Shot #167192")
    plt.show()
