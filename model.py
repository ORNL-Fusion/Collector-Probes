import pull_data_dp as pull
import numpy as np

# Define constants
#
# Probe widths in m
aSize = 3.0 / 10.0
bSize = 1.0 / 10.0
cSize = 0.5 / 10.0

# ne in 10^18 m^-2
def impur_model(Te=25, Ti=25, ne = 10):
    # Mass of tungsten in eV s^2 m^-2
    massW = 183.84 * 931.49 * 10**6.0 / ((3*10**8.0)**2.0)
    # Mass of Deuterium
    massD = 2.01 * 931.49 * 10**6 / ((3*10**8)**2.0)
    # Z of tungsten
    chargeW = 74.0
    # Perpendicular diffusion in m^2 s^-1
    dPerp = 1.0
    # Time of shot in s
    timeOfShot = 5.0

    # Plasma sound speed in m/s
    soundSpeed = ((Te + Ti)/massD)**(0.5)

    # Approx collection lengths of each probe in m
    aCollLength = aSize**2.0 * soundSpeed / (8.0 * dPerp)
    bCollLength = bSize**2.0 * soundSpeed / (8.0 * dPerp)
    cCollLength = cSize**2.0 * soundSpeed / (8.0 * dPerp)
    print "A collection length: " + str(aCollLength)
    print "B collection length: " + str(bCollLength)
    print "C collection length: " + str(cCollLength)

    # Impurity ion momentum collision frequency
    nu = 1.4*10**6 * (ne * chargeW**2)/(Te**(1.5)) * (massW / massD)

    # Get locations and areal density
    conn = pull.thin_connect(2, server='localhost')
    locations = []
    wAreal    = []
    for run in range(1, 21):
        # D-face data
        loc = pull.pull_rbs_loc(conn, 'AD', run)
        areal = pull.pull_rbs_areal(conn, 'AD', run)
        locations.append(loc)
        wAreal.append(areal)
        #print type(loc)
        #print type(wAreal[run-1])
        #print len(wAreal)
        print "Obtained data for run " + str(run)

    # First order flux approximation
    fluxes = []
    #print wAreal
    for spot in range(0, len(locations)):
        fluxes.append(wAreal[spot] / timeOfShot)

    # Value used in speed calculation
    alpha = chargeW * Te / (massW * aCollLength**2) + nu * soundSpeed / aCollLength
    # Speed calculation
    wSpeed = nu * aCollLength / 4 * (-1+(1+8*alpha/(nu**2.0))**(0.5))
    print "Nu: " + str(nu)
    print "wSpeed: " + str(wSpeed)
    print "alpha: " + str(alpha)
    print "a/nu : " + str(alpha/nu**2)
    print "massW: " + str(massW)
    print "massD: " + str(massD)
    print "soundSpeed: " + str(soundSpeed)

    # First order upstream density approximation (assumes const density along
    # flux tube).
    upDens = []
    for spot in range(0, len(fluxes)):
        upDens.append(fluxes[spot] / wSpeed)

    return upDens
