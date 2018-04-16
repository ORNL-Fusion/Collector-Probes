# Written in python 2.7 since the atomic library is written in it :(.

from __future__ import print_function
import openpyxl as xl
import numpy as np
import atomic.atomic as atomic
import scipy.integrate as integrate
from scipy import interpolate
import math


# Probe widths in m
aSize = 3.0 / 100.0
bSize = 1.0 / 100.0
cSize = 0.5 / 100.0
# Mass of tungsten in eV s^2 m^-2
massW = 183.84 * 931.49 * 10**6.0 / ((3*10**8.0)**2.0)
# Mass of Deuterium in eV s^2 m^-2
massD = 2.01 * 931.49 * 10**6 / ((3*10**8)**2.0)
# The time the probe was in for shots.
timeOfShots = 5 * 25
# Z of of tungsten and deuterium
Z_W = 74.0
Z_D = 1.0
# Perpendicular diffusion in m^2 s^-1
dPerp = 1.0
# Approx U_0 as the heat of sublimation of tungsten in eV
u0 = 8.68
# Alpha used in yield calculation.
alpha = 3.44


class SputteringModel:
    def __init__(self, shot=167192):
        self.shot = shot
        self.net_dict   = {}
        self.loss_dict  = {}
        self.total_dict = {}

    def returnArray(self, dataSheet, lowRange, highRange):
        """
        Accepts range of cells from a sheet to convert into a normal
        1D numpy array.
        """
        cells = dataSheet[lowRange:highRange]
        cells = np.transpose(cells)
        cells = np.reshape(cells, cells.size)
        values = [cell.value for cell in cells]
        return values

    def calc_net_flux(self):
        # Open Excel workbook.
        try:
            wb = xl.load_workbook("A2.xlsx")
        except:
            print("Error loading workbook.")

        # Get the sheets for each probe.
        A2_sheet = wb.get_sheet_by_name("A2")
        B2_sheet = wb.get_sheet_by_name("B2")
        C2_sheet = wb.get_sheet_by_name("C2")

        # Create list of probe and sheet names to loop through.
        probe_list = ["AD", "AU", "BD", "BU", "CD", "CU"]

        # Extract the cells then extract the values from them.
        for probe in probe_list:
            if probe=="AD":
                rminrsep = self.returnArray(A2_sheet, "A2", "A20")
                areal    = self.returnArray(A2_sheet, "C2", "C20")
            elif probe=="AU":
                rminrsep = self.returnArray(A2_sheet, "H2", "H20")
                areal    = self.returnArray(A2_sheet, "I2", "I20")
            elif probe=="BD":
                rminrsep = self.returnArray(B2_sheet, "F2", "F22")
                areal    = self.returnArray(B2_sheet, "C2", "C22")
            elif probe=="BU":
                rminrsep = self.returnArray(B2_sheet, "H2", "H22")
                areal    = self.returnArray(B2_sheet, "I2", "I22")
            elif probe=="CD":
                rminrsep = self.returnArray(C2_sheet, "F2", "F22")
                areal    = self.returnArray(C2_sheet, "C2", "C22")
            elif probe=="CU":
                rminrsep = self.returnArray(C2_sheet, "H2", "H22")
                areal    = self.returnArray(C2_sheet, "I2", "I22")
            else:
                print("Incorrect probe entry. Must be AD, AU, BD, BU, CD or CU.")

            # Convert arealD units from 10^15 cm^-2 to just m^-2.
            areal = [value*10**19 for value in areal]

            # A first order approximation assumes short shots such that net flux = areal density / time of shot
            tmp_net = []
            for index in range(0, len(areal)):
                tmp = areal[index] / timeOfShots
                tmp_net.append(tmp)

            # Put the net flux for this probe into the dict.
            self.net_dict[probe] = {"rminrsep":rminrsep, "flux":tmp_net}

    def calc_loss_flux(self, shotANDplunge="167192.1"):
        """
        Calculates the flux sputtered from the probe that does not immediately
        return. It is assumed that if the sputtered tungsten does not ionize
        immediately, it will be lost. If it does reionize, it will return and thus
        not contribute to the loss flux.

        shotANDplunge: In the form of shot.plunge (as a string).
        """

        wb = xl.load_workbook("recLPdata.xlsx", data_only=True)
        dataSheet = wb.get_sheet_by_name("Sheet1")

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
        else:
            return print("Incorrect shot/plunge.")

        times = self.returnArray(dataSheet, timeLow, timeHigh)
        dens  = self.returnArray(dataSheet, densLow, densHigh)
        temps = self.returnArray(dataSheet, tempLow, tempHigh)
        rmins = self.returnArray(dataSheet, rMinRsepLow, rMinRsepHigh)

        # Go from 10^18 m^-3 to just m^-3.
        for index in range(0, len(dens)):
            if dens[index] is None:
                continue
            else:
                dens[index] = dens[index] * 10**18

        # Plasma sound speed assuming Te = Ti.
        sound_speeds = [(temp*2 / massD)**0.5 for temp in temps]

        self.shot_and_plunge = shotANDplunge
        self.times = times
        self.dens = dens
        self.temps = temps
        self.rmins = rmins
        self.sound_speeds = sound_speeds

        # The flux of W off the probe due to sputtering. sputt_flux = yield * flux of dueterium.
        def sputt_flux(ne, Ti, Te):
            # Sputtering energy threshold of tungsten oxide in eV. Note pure W is 160 eV.
            eThresh = 65
            soundSpeed = ((float(Te) + float(Ti)) / massD)**0.5

            # Use lambda function for use in integrate,
            func = lambda E: 0.528 * alpha * Z_D * (massD / (u0*(massD + massW))) * 0.059 * (E+3*Ti) ** (1.0/3.0) * soundSpeed * ne * 2 * (E/3.1415)**0.5 * (1/float(Ti))**(1.5) * math.exp(-E/Ti)
            ans, err = integrate.quad(func, eThresh, np.inf)

            #print("Sputtered Flux:      " + str(ans))
            #print("Sputtered Flux Error: " + str(err/ans * 100) + "%")

            return ans


        for probe in ["A", "B", "C"]:
            # Use corresponding size for desired probe.
            if probe=="A":
                size = aSize
            elif probe=="B":
                size = bSize
            elif probe=="C":
                size = cSize
            else:
                print("Incorrect probe entry. Should be either A, B, or C.")

            print("Calculating loss flux for " + probe + " probes...")

            flux_loss = []
            for index in range(0, len(self.temps)):
                Te = self.temps[index]
                ne = self.dens[index]
                cs = self.sound_speeds[index]

                # Approx. speed of W entering flux tube.
                v0 = 0.5 * cs

                # Get the ionization rate coefficient for a specific temperature.
                ad = atomic.element('tungsten')
                temperatureRange = np.logspace(0,4,100)
                S = ad.coeffs['ionisation']
                f = interpolate.interp1d(temperatureRange, S(0, temperatureRange, ne))
                coeff = f(Te)

                # Calculate lamda_ionization.
                lambda_iz = v0 * (ne * coeff)**(-1)

                # Fraction ionized in the flux tube (i.e. it will return to the probe)
                frac = 1 - math.exp(-size / lambda_iz)
                #print("Fraction Ionized:    " + str(frac))

                # Thus the fraction lost is 1-frac of the sputtered flux.
                Ti = Te
                fracFluxLost = (1 - frac) * sputt_flux(ne=ne, Ti=Ti, Te=Te)
                #print("Flux Lost:           " + str(fracFluxLost))

                flux_loss.append(fracFluxLost)

            self.loss_dict[probe] = {"rminrsep":self.rmins, "flux":flux_loss}

    def calc_total_flux(self):

        probe_list = ["AD", "AU", "BD", "BU", "CD", "CU"]

        for probe in probe_list:
            tmp_net_rmin = self.net_dict[probe]["rminrsep"]
            tmp_net_flux = self.net_dict[probe]["flux"]
            tmp_loss_rmin = self.loss_dict[probe[0]]["rminrsep"]
            tmp_loss_flux = self.loss_dict[probe[0]]["flux"]

            tmp_f_net = interpolate.interp1d(tmp_net_rmin, tmp_net_flux)
            tmp_f_loss = interpolate.interp1d(tmp_loss_rmin, tmp_loss_flux)

            rminrsep_max = min(max(tmp_net_rmin), max(tmp_loss_rmin))
            rminrsep_min = max(min(tmp_net_rmin), min(tmp_loss_rmin))

            total_rminrsep = []
            total_flux = []
            for index in range(0, len(tmp_net_rmin)):
                #print(probe + ": " + str(tmp_net_rmin[index]) + " > " + str(rminrsep_min) + " and < " + str(rminrsep_max))
                if tmp_net_rmin[index] > rminrsep_min and tmp_net_rmin[index] < rminrsep_max:
                    print("At R-Rsep: " + str(tmp_net_rmin[index]))
                    print("Net Flux:   " + str(tmp_net_flux[index]))
                    print("Loss Flux:  " + str(tmp_f_loss(tmp_net_rmin[index])))
                    print("Total Flux: " + str(tmp_net_flux[index] + tmp_f_loss(tmp_net_rmin[index])))
                    print()
                    tmp_total = tmp_net_flux[index] + tmp_f_loss(tmp_net_rmin[index])
                    total_flux.append(tmp_total)
                    total_rminrsep.append(tmp_net_rmin[index])

            self.total_dict[probe] = {"rminrsep":total_rminrsep, "flux":total_flux}
