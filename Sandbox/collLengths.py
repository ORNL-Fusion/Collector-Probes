import openpyxl as xl
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import math
#import shawn_model as model

wb = xl.load_workbook("recLPdata.xlsx", data_only=True)
dataSheet = wb.get_sheet_by_name("Sheet1")

# Mass of Deuterium in eV s^2 m^-2
massD = 2.01 * 931.49 * 10**6 / ((3*10**8)**2.0)
massW = 183.84 * 931.49 * 10**6 / ((3*10**8)**2.0)
#massW = 15.999 * 931.49 * 10**6 / ((3*10**8)**2.0)
massElec = 0.511 * 10**6 / ((3*10**8)**2.0)
massElec_kg = 9.11*10**(-31)
massD_kg = 2.01 * 1.66*10**(-27)
massW_kg = 183.84 * 1.66*10**(-27)
# Z of each element.
zD = 1.0
zE = 1.0
zW = 5.5
# Perpendicular diffusion in m^2 s^-1
dPerp = 1.0
# Probe widths in m
aSize = 3.0 / 100.0
bSize = 1.0 / 100.0
cSize = 0.5 / 100.0
# Electric charge
elec = 1.609*10**(-19)
# Vacuum permittivity
epsilon = 8.85*10**(-12)
# Coulumb's constant[den * 0.03 for den in densitiesElec]
coulConst = 8.99 * 10**9


def returnArray(lowRange, highRange):
    cells = dataSheet[lowRange:highRange]
    cells = np.transpose(cells)
    cells = np.reshape(cells, cells.size)
    values = [cell.value for cell in cells]
    return values

def putN2dict(shotANDplunge):
    # Get the correct cells.
    if shotANDplunge == "167192.1":
        timeLow = "A3"
        timeHigh = "A64"
        densLow = "C3"
        densHigh = "C64"
        tempLow = "D3"
        tempHigh = "D64"
        rMinRsepLow = "G3"
        rMinRsepHigh = "G64"
    elif shotANDplunge == "167192.2":
        timeLow = "I3"
        timeHigh = "I55"
        densLow = "K3"
        densHigh = "K55"
        tempLow = "L3"
        tempHigh = "L55"
        rMinRsepLow = "O3"
        rMinRsepHigh = "O55"
    elif shotANDplunge == "167193.1":
        timeLow = "Q3"
        timeHigh = "Q61"
        densLow = "S3"
        densHigh = "S61"
        tempLow = "T3"
        tempHigh = "T61"
        rMinRsepLow = "W3"
        rMinRsepHigh = "W61"
    elif shotANDplunge == "167193.2":
        timeLow = "Y3"
        timeHigh = "Y48"
        densLow = "AA3"
        densHigh = "AA48"
        tempLow = "AB3"
        tempHigh = "AB48"
        rMinRsepLow = "AE3"
        rMinRsepHigh = "AE48"
    elif shotANDplunge == "167194.1":
        timeLow = "AG3"
        timeHigh = "AG71"
        densLow = "AI3"
        densHigh = "AI71"
        tempLow = "AJ3"
        tempHigh = "AJ71"
        rMinRsepLow = "AM3"
        rMinRsepHigh = "AM71"
    elif shotANDplunge == "167194.2":
        timeLow = "AO3"
        timeHigh = "AO67"
        densLow = "AQ3"
        densHigh = "AQ67"
        tempLow = "AR3"
        tempHigh = "AR67"
        rMinRsepLow = "AU3"
        rMinRsepHigh = "AU67"
    elif shotANDplunge == "167195.1":
        timeLow = "AW3"
        timeHigh = "AW60"
        densLow = "AY3"
        densHigh = "AY60"
        tempLow = "AZ3"
        tempHigh = "AZ60"
        rMinRsepLow = "BC3"
        rMinRsepHigh = "BC60"
    elif shotANDplunge == "167195.2":
        timeLow = "BE3"
        timeHigh = "BE59"
        densLow = "BG3"
        densHigh = "BG59"
        tempLow = "BH3"
        tempHigh = "BH59"
        rMinRsepLow = "BK3"
        rMinRsepHigh = "BK59"
    #else:
        #return print("Incorrect entry.")

    times = returnArray(timeLow, timeHigh)
    dens  = returnArray(densLow, densHigh)
    temps = returnArray(tempLow, tempHigh)
    rmins = returnArray(rMinRsepLow, rMinRsepHigh)

    # Go from 10^18 m^-3 to just m^-3.
    for index in range(0, len(dens)):
        if dens[index] is None:
            continue
        else:
            dens[index] = dens[index] * 10**18

    # Plasma sound speed assuming Te = Ti.
    soundSpeeds = [(temp*2 / massD)**0.5 for temp in temps]

    valuesDict = {"Shot.plunge" :shotANDplunge,
                  "Times"       :times,
                  "Densities"   :dens,
                  "Temperatures":temps,
                  "RminRSeps"   :rmins,
                  "Sound Speeds":soundSpeeds}

    return valuesDict


def avgAllPlunges():
    p167192n1 = putN2dict("167192.1")
    p167192n2 = putN2dict("167192.2")
    p167193n1 = putN2dict("167193.1")
    p167193n2 = putN2dict("167193.2")
    p167194n1 = putN2dict("167194.1")
    p167194n2 = putN2dict("167194.2")
    p167195n1 = putN2dict("167195.1")
    p167195n2 = putN2dict("167195.2")

    plungeList = [p167192n1, p167192n2, p167193n1, p167193n2,
                  p167194n1, p167194n2, p167195n1, p167195n2]

    # 167195.2 reaches the furthest in. Will use this out to it's outer bound (15.923),
    # and 167193.2 for the rest of the way (up until 16.124).
    # OR
    # R-Rsep range with each contributing is 15.775 - 4.953 cm. This encompasses all
    # the '2' probe data.

    # Make Te interpolations for each.
    te_fits = []
    for plunge in plungeList:
        fit = interpolate.interp1d(plunge["RminRSeps"], plunge["Temperatures"])
        te_fits.append(fit)

    # Make ne interpolations for each.
    ne_fits = []
    for plunge in plungeList:
        fit = interpolate.interp1d(plunge["RminRSeps"], plunge["Densities"])
        ne_fits.append(fit)

    # Find average of Te and ne between 4.953 - 15.775 cm. Put into list.
    avg_dict = {}
    avg_dict["RminRSeps"] = np.arange(5.077, 15.775, 0.1)

    tmp_list = []
    for r in np.arange(5.077, 15.775, 0.1):
        tmp_sum = 0
        #print(r)
        for fit in te_fits:
            tmp_sum = tmp_sum + fit(r)
        avg = tmp_sum / 8.0
        tmp_list.append(avg)
    avg_dict["Temperatures"] = tmp_list

    ne_list = []
    for r in np.arange(5.077, 15.775, 0.1):
        tmp_sum = 0
        #print(r)
        for fit in ne_fits:
            tmp_sum = tmp_sum + fit(r)
        avg = tmp_sum / 8.0
        ne_list.append(avg)
    avg_dict["Densities"] = ne_list

    # Do an average of the sound speeds as well.
    sound_speeds = []
    for temp in avg_dict["Temperatures"]:
        cs = (temp*2 / massD)**0.5
        sound_speeds.append(cs)
    avg_dict["Sound Speeds"] = sound_speeds

    return avg_dict


def collectionLengths(valuesDict):
    collectionLengthsA = [speed * aSize**2 / (8 * dPerp) for speed in valuesDict["Sound Speeds"]]
    collectionLengthsB = [speed * bSize**2 / (8 * dPerp) for speed in valuesDict["Sound Speeds"]]
    collectionLengthsC = [speed * cSize**2 / (8 * dPerp) for speed in valuesDict["Sound Speeds"]]

    valuesDict["A Collection Lengths"] = collectionLengthsA
    valuesDict["B Collection Lengths"] = collectionLengthsB
    valuesDict["C Collection Lengths"] = collectionLengthsC

    return valuesDict


def calcCollisionLengths(valuesDict):
    const = elec**4 / (4 * 3.1415 * epsilon)**2
    densities = valuesDict["Densities"]
    temps = valuesDict["Temperatures"]

    mfps = []
    for index in range(0, len(densities)):
        if temps[index] is None:
            mfps.append(None)
        elif densities[index] is None:
            mfps.append(None)
        else:
            tmp = 10**16 * temps[index]**2 / densities[index]
            mfps.append(tmp)

    valuesDict["Mean Free Paths"] = mfps

    return valuesDict

def calcDeflTime(valuesDict):

    # Electron density and temperatures.
    densitiesElec = valuesDict["Densities"]
    temps = valuesDict["Temperatures"]
    # Need an estimate of the carbon and deuterium density.
    densitiesCarbon = [den * 0.03 for den in densitiesElec]
    densitiesDeuterium = [den * 0.03 for den in densitiesElec]
    # Then use these to estimate the tugnsten density where using particle
    # balance we say: ne = nD + nC/Z_C + nW/Z_W.


    # Is speed sqrt(2kT/m) or the sound speed? <-- for tungsten
    # It is the first one.
    speedsW = []
    for index in range(0, len(temps)):
        if temps[index] is None:
            speedsW.append(None)
        else:
            tmp = (2 * temps[index] / massW)**(0.5)
            speedsW.append(tmp)

    # Same question for electrons
    speedsElec = []
    for index in range(0, len(temps)):
        if temps[index] is None:
            speedsElec.append(None)
        else:
            tmp = (2 * temps[index] / massElec)**(0.5)
            #print(tmp)
            speedsElec.append(tmp)

    # And same questions for ions.
    speedsD = []
    for index in range(0, len(temps)):
        if temps[index] is None:
            speedsD.append(None)
        else:
            tmp = (2 * temps[index] / massD)**(0.5)
            speedsD.append(tmp)

    # The lfs which are just me or mi / 2kT
    lfs_e = []
    lfs_i = []
    for index in range(0, len(temps)):
        if temps[index] is None:
            lfs_e.append(None)
            lfs_i.append(None)
        else:
            lfs_e.append((massElec / (2 * temps[index]))**(0.5))
            lfs_i.append((massD / (2 * temps[index]))**(0.5))

    #print(lfs_e)

    # Calculate alpha. Since zE = zD alpha will be the same for elec or D+.
    alphas = []
    for index in range(0, len(temps)):
        if temps[index] is None:
            alphas.append(None)
        elif densitiesElec[index] is None:
            alphas.append(None)
        else:
            # Need the mean Z of tungsten.
            tmp = temps[index]
            meanZofW = model.meanZofTungsten(tmp)
            alpha = 3 / (2 * meanZofW * zE * elec**3) * (temps[index]**3 / (3.1415 * densitiesElec[index]))**(1/2)
            alphas.append(alpha)
            #print(math.log(alpha))

    #print(alphas)

    # Calculate p0. This is under the assumption tungsten can be considered
    # stationary compared to ions and electrons.
    p0s_e = []
    p0s_i = []
    for index in range(0, len(speedsElec)):
        if speedsElec[index] is None:
            p0s_e.append(None)
            p0s_i.append(None)
        else:

            densElec = densitiesElec[index]
            meanZofTungsten = model.meanZofTungsten(evalTemp=temps[index], density=densElec)
            p0_e = coulConst * zE * meanZofTungsten * elec**2 / (massElec_kg * speedsElec[index]**2)
            p0_i = coulConst * zD * meanZofTungsten * elec**2 / (massD_kg * speedsD[index]**2)
            #print(p0_e)
            p0s_e.append(p0_e)
            p0s_i.append(p0_i)

    #print(p0s_e)
    # Need to interpolate for a function for psi - G. Table 5.2.
    xValues = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
               1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
               2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0,]
    psiMinGValues = [0.000, 0.750, 0.149, 0.221, 0.292, 0.358, 0.421, 0.480, 0.534, 0.584,
               0.629, 0.669, 0.706, 0.736, 0.766, 0.791, 0.813, 0.832, 0.849, 0.863,
               0.876, 0.920, 0.944, 0.959, 0.969, 0.980, 0.986, 0.990, 0.992, 0.995]
    f = interpolate.interp1d(xValues, psiMinGValues)
    psiMinGs_e = []
    psiMinGs_i = []
    for index in range(0, len(speedsW)):
        if speedsW[index] is None:
            psiMinGs_e.append(None)
            psiMinGs_i.append(None)
        else:
            lfTimesW_e = lfs_e[index] * speedsW[index]
            lfTimesW_i = lfs_i[index] * speedsW[index]
            psiMinGs_e.append(f(lfTimesW_e)[()])
            psiMinGs_i.append(f(lfTimesW_i)[()])

            #print(f(lfTimesW_e))

    #print(psiMinGs_e)

    # Then the deflection time is just (assuming ne = ni):
    deflecTimes_e = []
    deflecTimes_i = []
    for index in range(0, len(densitiesElec)):
        if densitiesElec[index] is None:
            deflecTimes_e.append(None)
            deflecTimes_i.append(None)
        elif temps[index] is None:
            deflecTimes_e.append(None)
            deflecTimes_i.append(None)
        else:
            #print(p0s_e[index])
            #time_e = (8 * 3.1415 * densitiesElec[index] * speedsW[index] * p0s_e[index]**2 * psiMinGs_e[index] * math.log(alphas[index]))**(-1)
            #time_i = (8 * 3.1415 * densitiesElec[index] * speedsW[index] * p0s_i[index]**2 * psiMinGs_i[index] * math.log(alphas[index]))**(-1)
            # Assume ln(alpha) = 15
            time_e = (8 * 3.1415 * densitiesElec[index] * speedsW[index] * p0s_e[index]**2 * psiMinGs_e[index] * 15)**(-1)
            time_i = (8 * 3.1415 * densitiesElec[index] * speedsW[index] * p0s_i[index]**2 * psiMinGs_i[index] * 15)**(-1)

            deflecTimes_e.append(time_e)
            deflecTimes_i.append(time_i)

    valuesDict["Elec Defl Times"] = deflecTimes_e
    valuesDict["Duet Defl Times"] = deflecTimes_i

    # Length to do a 90 degree deflection
    collisionLength_e = []
    collisionLength_i = []
    for index in range(0, len(deflecTimes_e)):
        if deflecTimes_e[index] is None:
            collisionLength_e.append(None)
            collisionLength_i.append(None)
        else:
            # Either thermal speed or sound speed
            tmp_e = deflecTimes_e[index] * speedsW[index]
            tmp_i = deflecTimes_i[index] * speedsW[index]
            #tmp_e = deflecTimes_e[index] * valuesDict["Sound Speeds"][index]
            #tmp_i = deflecTimes_i[index] * valuesDict["Sound Speeds"][index]
            collisionLength_e.append(tmp_e)
            collisionLength_i.append(tmp_i)

    valuesDict["Elec Collision Length"] = collisionLength_e
    valuesDict["Deut Collision Length"] = collisionLength_i
    #print(collisionLength_i)

    return valuesDict


def plotCollLengths(valuesDict):
    rminrseps = valuesDict["RminRSeps"]
    aCollLengths = valuesDict["A Collection Lengths"]
    bCollLengths = valuesDict["B Collection Lengths"]
    cCollLengths = valuesDict["C Collection Lengths"]
    deutCollisionLength = valuesDict["Stacey Lengths"]
    #elecCollisionLength = valuesDict["Elec Collision Length"]
    #deutCollisionLengthSpit = valuesDict["Deut Collision Length"]

    plt.semilogy(rminrseps, aCollLengths, label="A2")
    plt.semilogy(rminrseps, bCollLengths, label="B2")
    plt.semilogy(rminrseps, cCollLengths, label="C2")
    #plt.semilogy(rminrseps, deutCollisionLengthSpit, label=r'$L_{D^+-W}$ Collision Length Spit')
    plt.semilogy(rminrseps, deutCollisionLength, label=r'$L_{D^+-W}$ Collision Length')
    plt.xlabel("R-Rsep (cm)")
    plt.ylabel("Perturbation Length (m)")
    plt.title("'2' Probes Perturbation Lengths")
    plt.legend()
    plt.show()


def staceyDeflLengths(valuesDict):
    reduced_mass = massW_kg * massD_kg / (massW_kg + massD_kg)
    print reduced_mass
    ln_alpha = 15

    lp_temps = valuesDict["Temperatures"]
    lp_dens = valuesDict["Densities"]

    # Time between deflections in the CM.
    defl_times_cm = np.array([])
    for index in range(0, len(lp_temps)):
        tmp_temp = lp_temps[index]
        tmp_dens = lp_dens[index]
        tmp_time_cm = 2 * 3.1415 * math.sqrt(reduced_mass) * epsilon**2 * (3 * tmp_temp)**(1.5) / (tmp_dens * (zD * elec * zW * elec)**2 * ln_alpha)
        defl_times_cm = np.append(defl_times_cm, tmp_time_cm)

    # To make unit works need to multiply by (1.609x10^18)^3/2
    defl_times_cm = defl_times_cm * (elec**(1.5))

    valuesDict["Stacey Times CM"] = defl_times_cm

    # But we want the time between deflections in the LAB frame.
    defl_times_lab = defl_times_cm * (massW / massD)
    valuesDict["Stacey Times Lab"] = defl_times_lab

    # Electron temperatures.
    temps = valuesDict["Temperatures"]

    # Speed of tungsten assuming Tw = Ti = Te.
    fudge_factor = 1.0
    speedsW = []
    for index in range(0, len(temps)):
        if temps[index] is None:
            speedsW.append(None)
        else:
            tmp = (2 * temps[index] * fudge_factor / massW)**(0.5)
            speedsW.append(tmp)
    valuesDict["W Speeds"] = speedsW

    # Deflection length is just speed * time between deflections.
    defl_lengths = np.array([])
    for index in range(0, len(defl_times_lab)):
        tmp_time = defl_times_lab[index]
        tmp_speed = speedsW[index]
        tmp_length = tmp_time * tmp_speed
        defl_lengths = np.append(defl_lengths, tmp_length)

    valuesDict["Stacey Lengths"] = defl_lengths

    return valuesDict

def runScript():
    # Get the average plunging LP data.
    valuesDict = avgAllPlunges()

    # Get the collection lengths of the probes.
    valuesDict = collectionLengths(valuesDict)

    # Get the Spitzer deflection lengths
    #valuesDict = calcDeflTime(valuesDict)

    # Get the Stacey collision lengths.
    valuesDict = staceyDeflLengths(valuesDict)

    plotCollLengths(valuesDict)

    return valuesDict
