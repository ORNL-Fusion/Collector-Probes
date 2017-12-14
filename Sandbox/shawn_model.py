# This program determines the total flux to a collector probe (maxFlux). The
# basis is net_flux = maxFlux - loss_flux, where loss_flux is due to sputtering.


from __future__ import print_function
import numpy as np
import openpyxl as xl
import collLengths as coll
import scipy.integrate as integrate
import math
import atomic.atomic as atomic
from scipy import interpolate
import matplotlib.pyplot as plt


# Define constants.
# Probe widths in m
aSize = 3.0 / 100.0
bSize = 1.0 / 100.0
cSize = 0.5 / 100.0
# Mass of tungsten in eV s^2 m^-2
massW = 183.84 * 931.49 * 10**6.0 / ((3*10**8.0)**2.0)
# Mass of Deuterium
massD = 2.01 * 931.49 * 10**6 / ((3*10**8)**2.0)
# Charge state of tungsten and deuterium
chargeW = 74.0
chargeD = 1.0
# Perpendicular diffusion in m^2 s^-1
dPerp = 1.0
# Approx U_0 as the heat of sublimation of tungsten in eV
u0 = 8.68
# Alpha used in yield calculation.
alpha = 3.44
# Sound speed for Ti=Te=25 eV for reference
ref_cs = ((25 + 25) / massD)**(0.5)
# Time of shot in s
timeOfShot = 5.0 * 25

valuesDict = coll.avgAllPlunges()
rminrseps  = valuesDict["RminRSeps"]
temps = valuesDict["Temperatures"]
densities = valuesDict["Densities"]

# Open up the Excel file with all the 2 probe data.
wb = xl.load_workbook("A2.xlsx")
A2_sheet = wb.get_sheet_by_name("A2")
B2_sheet = wb.get_sheet_by_name("B2")
C2_sheet = wb.get_sheet_by_name("C2")


def meanZofTungsten(evalTemp, LPtemps=temps, density=1e19):
    """ Uses the atomic library to calculate the mean charge state of tungsten
            using the temperature and density data from the LP's."""

    ad = atomic.element("tungsten")
    eq = atomic.CollRadEquilibrium(ad)
    y = eq.ionisation_stage_distribution(LPtemps, density)

    f = interpolate.interp1d(LPtemps, y.mean_charge())
    meanChargeAtT = f(evalTemp)

    return meanChargeAtT

def meanZofCarbon(evalTemp, LPtemps=temps, density=1e19):
    """ Uses the atomic library to calculate the mean charge state of carbon
            using the temperature and density data from the LP's."""

    ad = atomic.element("carbon")
    eq = atomic.CollRadEquilibrium(ad)
    y = eq.ionisation_stage_distribution(LPtemps, density)

    f = interpolate.interp1d(LPtemps, y.mean_charge())
    meanChargeAtT = f(evalTemp)

    return meanChargeAtT


# Function for the yield of sputtering for D+ on W from Was textbook.
def Y(energy):
    yield_func = 0.528 * alpha * chargeD * (massD / (u0*(massD + massW))) * 0.059 * energy ** (1/3)
    return yield_func


# The flux at the probe as a function of energy assuming Maxwellian distribution
# and a constant Deuterium density in the flux tube (so n(E) -> n).
def fluxD(energy, Ti, ne):
    flux_func = ref_cs * ne * 2 * (energy/3.1415)**0.5 * (1/float(Ti))**(1.5) * math.exp(-energy/Ti)
    return flux_func


# The flux of W off the probe due to sputtering. sputt_flux = yield * flux of dueterium.
def sputt_flux(ne=10**18, Ti=25.0, Te=25.0):
    # Sputtering energy threshold of tungsten oxide in eV. Note pure W is 160 eV.
    eThresh = 65
    soundSpeed = ((float(Te) + float(Ti)) / massD)**0.5

    # Use lambda function for use in integrate,
    func = lambda E: 0.528 * alpha * chargeD * (massD / (u0*(massD + massW))) * 0.059 * (E+3*Ti) ** (1.0/3.0) * soundSpeed * ne * 2 * (E/3.1415)**0.5 * (1/float(Ti))**(1.5) * math.exp(-E/Ti)
    ans, err = integrate.quad(func, eThresh, np.inf)

    print("Sputtered Flux:      " + str(ans))
    #print("Sputtered Flux Error: " + str(err/ans * 100) + "%")

    return ans


# The loss_flux is that which is sputtered and does NOT return to the probe. It
# is assumed if the sputtered W ionizes in the flux tube it will return to the
# probe.
def loss_flux(ne=10**18, Ti=25.0, Te=25.0, probe="A"):
    # Use corresponding size for desired probe.
    if probe=="A":
        size = aSize
    elif probe=="B":
        size = bSize
    elif probe=="C":
        size = cSize
    else:
        print("Incorrect probe entry. Should be either A, B, or C.")

    # Get the ionization rate coefficient for a specific temperature.
    ad = atomic.element('tungsten')
    temperatureRange = np.logspace(0,4,100)
    S = ad.coeffs['ionisation']
    f = interpolate.interp1d(temperatureRange, S(0, temperatureRange, ne))
    coeff = f(Te)

    # Initial speed entering the flux tube approx. v0.
    soundSpeed = ((float(Te) + float(Ti)) / massD)**0.5
    v0 = 0.5 * soundSpeed

    # Calculate lamda_ionization.
    lambda_iz = v0 * (ne * coeff)**(-1)

    # Fraction ionized in the flux tube (i.e. it will return to the probe)
    frac = 1 - math.exp(-size / lambda_iz)
    print("Fraction Ionized:    " + str(frac))

    # Thus the fraction lost is 1-frac of the sputtered flux.
    fracFluxLost = (1 - frac) * sputt_flux(ne=ne, Ti=Ti, Te=Te)
    print("Flux Lost:           " + str(fracFluxLost))

    return fracFluxLost


# net_flux is defined as maxFlux - loss_flux. It can also be approximated
# as net_flux = areal density of W / timeOfShot. This may be too rough an
# approximation though.
def net_flux(probe="AD"):
    # Choose correct Excel file sheet
    if probe[0]=="A":
        sheet = A2_sheet
    elif probe[0]=="B":
        sheet = B2_sheet
    elif probe[0]=="C":
        sheet = C2_sheet

    # Extract the cells then extract the values from them.
    if probe=="AD":
        rminrsep_cells = sheet["A2":"A20"]
        areal_cells    = sheet["C2":"C20"]
    elif probe=="AU":
        rminrsep_cells = sheet["H2":"H20"]
        areal_cells    = sheet["I2":"I20"]
    elif probe=="BD":
        rminrsep_cells = sheet["F2":"F22"]
        areal_cells    = sheet["C2":"C22"]
    elif probe=="BU":
        rminrsep_cells = sheet["H2":"H22"]
        areal_cells    = sheet["I2":"I22"]
    elif probe=="CD":
        rminrsep_cells = sheet["F2":"F22"]
        areal_cells    = sheet["C2":"C22"]
    elif probe=="CU":
        rminrsep_cells = sheet["H2":"H22"]
        areal_cells    = sheet["I2":"I22"]
    else:
        print("Incorrect probe entry. Must be AD, AU, BD, BU, CD or CU.")

    rminrsep_cells = np.transpose(rminrsep_cells)
    areal_cells    = np.transpose(areal_cells)
    rminrsep_cells = np.reshape(rminrsep_cells, rminrsep_cells.size)
    areal_cells    = np.reshape(areal_cells, areal_cells.size)
    rminrsep = [cell.value for cell in rminrsep_cells]
    areal    = [cell.value for cell in areal_cells]

    # Convert arealD units from 10^15 cm^-2 to just m^-2.
    areal = [value*10**19 for value in areal]

    # A first order approximation assumes short shots such that net flux = areal density / time of shot
    tmp_net = []
    for index in range(0, len(areal)):
        tmp = areal[index] / timeOfShot
        #print ("Net Flux: " + str(tmp) + "\n")
        tmp_net.append(tmp)

    net_dict = {"rminrsep":rminrsep, "net_flux":tmp_net, "areal":areal}

    return net_dict


# max_flux is the total flux of W to the probe. It can be determined from
# max_flux = loss_flux + net_flux.
def max_flux(probe="AD"):
    # Get the net flux dictionary, give the same rminrsep to the max flux dict.
    net_dict = net_flux(probe=probe)
    max_dict = {}
    max_dict["rminrsep"] = net_dict["rminrsep"]

    # Interpolations from the LP's for rminrsep vs. temp and density
    f_LP_temps = interpolate.interp1d(rminrseps, temps)
    f_LP_dens = interpolate.interp1d(rminrseps, densities)

    # Fill in the max_dict with max_flux = net flux + loss flux
    max_dict["max_flux"] = []
    max_dict["densityW"] = []
    max_dict["loss_flux"] = []
    max_dict["net_flux"] = []
    for index in range(0, len(max_dict["rminrsep"])):
        tmp_rmrs = max_dict["rminrsep"][index]

        # If the rminrsep is out of range of the interpolation, just put a zero and continue.
        if tmp_rmrs > max(rminrseps):
            max_dict["max_flux"].append(0)
            max_dict["densityW"].append(0)
            max_dict["net_flux"].append(0)
            max_dict["loss_flux"].append(0)

            continue
        if tmp_rmrs  < min(rminrseps):
            max_dict["max_flux"].append(0)
            max_dict["densityW"].append(0)
            max_dict["net_flux"].append(0)
            max_dict["loss_flux"].append(0)
            continue

        # Temporary values of the temp/density at the specified rminrsep.
        tmp_temp = f_LP_temps(tmp_rmrs)
        tmp_dens = f_LP_dens(tmp_rmrs)
        tmp_net  = net_dict["net_flux"][index]

        # Put the net flux in the dict as well.
        max_dict["net_flux"].append(tmp_net)

        # Get the loss flux at the specified parameters.
        tmp_loss = loss_flux(ne=tmp_dens, Ti=tmp_temp, Te=tmp_temp, probe=probe[0])

        # Put the loss_flux in the dict as well.
        max_dict["loss_flux"].append(tmp_loss)

        # Add net flux and loss flux and put into max flux dict.
        max_dict["max_flux"].append(tmp_net + tmp_loss)

        # We can estimate the tungsten density in the corresponding flux tube under
        # our assumptions as densityW = max_flux / soundSpeed.
        soundSpeed = ((float(tmp_temp) + float(tmp_temp)) / massD)**0.5
        tmp_wdens = (tmp_net + tmp_loss) / soundSpeed
        max_dict["densityW"].append(tmp_wdens)
        print("Percent W in Flux Tube: " + str(tmp_wdens / tmp_dens * 100) + "%")

        if tmp_net == 0:
            continue
        else:
            print("Percent Flux Lost:   " + str(tmp_loss / tmp_net * 100.0) + "%")

        print("Net Flux:            " + str(tmp_net))


        print("\n")

    # Put the W areal density in max_dict just because.
    max_dict["areal"] = net_dict["areal"]

    return max_dict


def plotFluxes(probe):
    myDict = max_flux(probe)

    x = myDict["rminrsep"]
    max = myDict["max_flux"]
    net = myDict["net_flux"]
    loss = myDict["loss_flux"]

    plt.plot(x, max, label="Real Flux Onto Probe")
    plt.plot(x, net, label="Net Flux (Real - Loss)")
    plt.plot(x, loss, label="Loss Flux Off Probe Due to Sputtering")
    plt.legend()
    plt.xlabel(r"${\rm R - R_{sep}\ (cm)}$")
    plt.ylabel(r"${\rm Flux\ (cm^{-2} s^{-1})}$")
    plt.title("Comparision of Fluxes for " + probe + "2")
    plt.show()
