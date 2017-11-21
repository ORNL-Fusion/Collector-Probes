import numpy as np
import openpyxl as xl
import collLengths as coll
import scipy.integrate as integrate
import math

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

valuesDict = coll.putN2dict("167192.1")
rminrseps  = valuesDict["RminRSeps"]
temps = valuesDict["Temperatures"]
densities = valuesDict["Densities"]

# Function for the yield of sputtering for D+ on W from Was textbook.
def Y(energy):
    yield_func = 0.528 * alpha * chargeD * (massD / (u0*(massD + massW))) * 0.059 * energy ** (1/3)
    return yield_func


# The flux at the probe as a function of energy assuming Maxwellian distribution
# and a constant Deuterium density in the flux tube (so n(E) -> n).
def fluxD(energy, Ti, ne):
    flux_func = ref_cs * ne * 2 * (energy/3.1415)**0.5 * (1/Ti)**(3/2) * math.exp(-energy/Ti)
    return flux_func


def sputt_flux(ne=10**18, Ti=25, soundSpeed=ref_cs):
    # Sputtering energy threshold of tungsten oxide in eV. Note pure W is 160 eV.
    eThresh = 65

    func = lambda E: 0.528 * alpha * chargeD * (massD / (u0*(massD + massW))) * 0.059 * (E+3*Ti) ** (1/3) * soundSpeed * ne * 2 * (E/3.1415)**0.5 * (1/Ti)**(3/2) * math.exp(-E/Ti)

    ans, err = integrate.quad(func, eThresh, 1000000)

    print(ans)
    print (err)

    return ans

def loss_flux(sputtFlux):
    frac = 
