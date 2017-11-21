import pull_data_dp as pull
import numpy as np
from scipy import interpolate

# Requires the atomic library for rate coefficients. Can be found at:
# https://github.com/cfe316/atomic
# Place a soft link in the Collector-Probe folder to the repository:
# ln -s /path/to/atomic/
import atomic.atomic as atomic

# Probe widths in m
aSize = 3.0 / 100.0
bSize = 1.0 / 100.0
cSize = 0.5 / 100.0
# Mass of tungsten in eV s^2 m^-2
massW = 183.84 * 931.49 * 10**6.0 / ((3*10**8.0)**2.0)
# Mass of Deuterium
massD = 2.01 * 931.49 * 10**6 / ((3*10**8)**2.0)
# Z of tungsten
chargeW = 74.0
# Perpendicular diffusion in m^2 s^-1
dPerp = 1.0

# ne in 10^18 m^-2
def stange_impur_model(Te=25, Ti=25, ne = 10):
    # Time of shot in s
    timeOfShot = 5.0

    # Plasma sound speed in m/s
    soundSpeed = ((Te + Ti)/massD)**(0.5)
    # Initial speed entering the flux tube approx.
    v0 = 0.5 * soundSpeed


    # Approx collection lengths of each probe in m
    aCollLength = aSize**2.0 * soundSpeed / (8.0 * dPerp)
    bCollLength = bSize**2.0 * soundSpeed / (8.0 * dPerp)
    cCollLength = cSize**2.0 * soundSpeed / (8.0 * dPerp)
    print "A collection length: " + str(aCollLength)
    print "B collection length: " + str(bCollLength)
    print "C collection length: " + str(cCollLength)

    # Impurity ion momentum collision frequency
    nu = 1.4*10**6 * (ne * chargeW**2)/(Te**(1.5)) * (massW / massD)

    # Get locations and areal density (W/cm^2)
    conn = pull.thin_connect(2, server='localhost')
    locations = []
    wAreal    = []
    for run in range(1, 21):
        # D-face data
        loc = pull.pull_rbs_loc(conn, 'AD', run)
        areal = pull.pull_rbs_areal(conn, 'AD', run)
        locations.append(loc)
        # Store as W/m^2
        wAreal.append(areal * 10**(4))
        #print type(loc)
        #print type(wAreal[run-1])
        #print len(wAreal)
        print "Obtained data for run " + str(run)

    # Value used in speed calculation
    alpha = chargeW * Te / (massW * aCollLength**2) + nu * soundSpeed / aCollLength
    # Speed calculation in m/s at the end of the collection length.
    wSpeed = nu * aCollLength / 4 * (-1+(1+8*alpha/(nu**2.0))**(0.5))

    print "Nu: " + str(nu)
    print "wSpeed: " + str(wSpeed)
    print "alpha: " + str(alpha)
    print "a/nu : " + str(alpha/nu**2)
    print "massW: " + str(massW)
    print "massD: " + str(massD)
    print "soundSpeed: " + str(soundSpeed)

    # First order flux approximation
    fluxes = []
    for spot in range(0, len(locations)):
        fluxes.append(wAreal[spot] / timeOfShot)
    # First order upstream density approximation (assumes const density along
    # flux tube).
    upDens = []
    for spot in range(0, len(fluxes)):
        upDens.append(fluxes[spot] / wSpeed)

    # Get the ionization rate coefficient for a specific temperature.
    ad = atomic.element('tungsten')
    temperatureRange = np.logspace(0,4,100)
    S = ad.coeffs['ionisation']
    f = interpolate.interp1d(temperatureRange, S(0, temperatureRange, ne))
    coeff = f(Te)
    print "Ionization Coeff: " + str(coeff)

    # Calculate lamda_ionization.
    lambda_iz = v0 * (ne*10**18 * coeff)**(-1)
    print "Lambda_iz : " + str(lambda_iz)

    # Estimation of tungsten density at end of flux tube.
    nWinf = 0.01 * ne * 10**18
    # Flux of W at end of flux tube.
    fluxWinf = nWinf * wSpeed
    # Estimation for short exposures.
    imp_dep = fluxWinf * timeOfShot

    # Need Te from LPs

    return imp_dep
