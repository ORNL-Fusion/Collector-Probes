"""
Author: Shawn Zamperini
Email:  zamp@utk.edu

This is a collection of plotting functions for plotting the input/output data
from a 3DLIM run.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import netCDF4
from scipy.optimize import curve_fit


# Some plot properties to make them a bit nicer.
plt.ion()
plt.rcParams['font.family'] = 'serif'
fontsize = 12
ms = 2
lw = 5
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the tableau20 RGBs to numbers between (0,1) since this is how mpl accepts them.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)


class LimPlots:
    """
    Class object to store data and plotting routines for a 3DLIM run.
    """

    def __init__(self, ncpath):
        """
        Just intialize with the netCDF file. Will assume the lim and dat file
        are in the same folder.
        """

        # Just a default file for testing.
        if ncpath == 'test':
            ncpath = 'colprobe-z1-001e.nc'

        # Load in netcdf file.
        self.nc = netCDF4.Dataset(ncpath)

        # Load in the lim file.
        limpath = ncpath[:-2] + 'lim'
        with open(limpath) as f:
            self.lim = f.read()

        # load in the dat file.
        datpath = ncpath[:-2] + 'dat'
        with open(datpath) as f:
            self.dat = f.read()

        # Save file name.
        self.file = ncpath.split('/')[-1][:-3]

        # Save case name.
        self.case = ''
        for b in self.nc.variables['TITLE'][:].data:
            self.case = self.case + b.decode('UTF-8')

    def __repr__(self):

        message = 'LimPlots Object\n' + \
                  '  Case: ' + self.case + '\n' + \
                  '  File: ' + self.file + '\n'

        return message

    def summary(self):

        # Output dictionary to print out results easier.
        output = dict()

        # Case info.
        output['Case'] = self.case
        output['File'] = self.file

        # Time for run in hours.
        time = int(self.dat.split('TOTAL CPU TIME USED     (S)')[1].split('\n')[0])
        output['Time'] = str(time) + 's (' + format(time/3600, '.2f') + ' hours)'

        # Number of impurities followed.
        num = int(self.dat.split('NO OF IMPURITY IONS TO FOLLOW')[1].split('\n')[0])
        output['Ions Followed'] = "{:,}".format(num)

        # Find longest output for formatting.
        pad = 0
        for val in output.values():
            if len(str(val)) > pad:
                pad = len(str(val))

        # Printing commands.
        num_stars = 2 + 15 + 2 + pad
        print("\n" + "*" * num_stars)
        for key, val in output.items():
            print("* {:15}{:<{pad}} *".format(key, val, pad=pad))
        print("*" * num_stars)

    def get_dep_array(self):
        """
        Load the deposition arrays for the collector probes.

        To-Do
        - Add option to combine the arrays of multiple repeat runs.
        """

        # Only load it once. Keep track if it's already been loaded by trying
        # to see if it's been defined yet.
        try:
            self.dep_arr

        # Not defined, so load it.
        except AttributeError:

            # Create the deposition array for the initial file.
            dep_arr = np.array(self.nc.variables['NERODS3'][0] * -1)

            # Define dep_arr so next time you won't have to choose all the file
            # locations.
            self.dep_arr = dep_arr

        return self.dep_arr

    def centerline(self, log=False, fit_exp=False, plotnum=0):
        """
        Plot the ITF and OTF deposition along the centerlines on the same plot.

        log:       Option to make y axis a log scale.
        fit_exp:   Do an exponential fit onto the data and get the lambdas.

        To-Do
        - Add option so ITF/OTF is only over, say, first 5 cm.
        """

        #The deposition array.
        dep_arr = self.get_dep_array()

        # Location of each P bin, and its width.
        ps     = np.array(self.nc.variables['PS'][:].data)
        pwids  = np.array(self.nc.variables['PWIDS'][:].data)

        # Array of poloidal locations (i.e. the center of each P bin).
        pol_locs = ps - pwids/2.0

        # Distance cell centers along surface (i.e. the radial locations).
        rad_locs = np.array(self.nc.variables['ODOUTS'][:].data)

        # Get the centerline index (or closest to it).
        cline = np.abs(pol_locs).min()

        # Index the deposition array at the centerline for plotting.
        itf_x = rad_locs[np.where(rad_locs > 0.0)[0]]
        itf_y = dep_arr[np.where(pol_locs == cline)[0], np.where(rad_locs > 0.0)[0]]
        otf_x = rad_locs[np.where(rad_locs < 0.0)[0]] * -1
        otf_y = dep_arr[np.where(pol_locs == cline)[0], np.where(rad_locs < 0.0)[0]]

        # Plotting commands.
        if plotnum == 0:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            ax = self.master_fig.axes[plotnum-1]

        # Option for a log axis.
        if log:
            ax.semilogy(itf_x*100, itf_y, '-', label='ITF', ms=ms, color=tableau20[6])
            ax.semilogy(otf_x*100, otf_y, '-', label='OTF', ms=ms, color=tableau20[8])
        else:
            ax.plot(itf_x*100, itf_y, '-', label='ITF', ms=ms, color=tableau20[6])
            ax.plot(otf_x*100, otf_y, '-', label='OTF', ms=ms, color=tableau20[8])

        ax.legend(fontsize=fontsize)
        ax.set_xlabel('Distance along probe (cm)', fontsize=fontsize)
        ax.set_ylabel('Deposition (arbitrary units)', fontsize=fontsize)
        ax.set_xlim([0, 10])
        ax.set_ylim([0,None])

        # Option to perform an exponential fit to the data.
        if fit_exp:
            def exp_fit(x, a, b):
                return a * np.exp(-b * x)

            popt_itf, pcov_itf = curve_fit(exp_fit, itf_x, itf_y, maxfev=5000)
            popt_otf, pcov_otf = curve_fit(exp_fit, otf_x, otf_y, maxfev=5000)

            fitx = np.linspace(0, 0.1, 100)
            fity_itf = exp_fit(fitx, *popt_itf)
            fity_otf = exp_fit(fitx, *popt_otf)

            if log:
                ax.semilogy(fitx*100, fity_itf, '--', ms=ms, color=tableau20[6])
                ax.semilogy(fitx*100, fity_otf, '--', ms=ms, color=tableau20[8])
            else:
                ax.plot(fitx*100, fity_itf, '--', ms=ms, color=tableau20[6])
                ax.plot(fitx*100, fity_otf, '--', ms=ms, color=tableau20[8])

            print("Lambdas")
            print("  ITF = {:.2f}".format(1/popt_itf[1]*100))
            print("  OTF = {:.2f}".format(1/popt_otf[1]*100))

        if plotnum ==0:
            fig.tight_layout()
            fig.show()

        print("Center ITF/OTF: {:.2f}".format(itf_y.sum()/otf_y.sum()))

    def deposition_contour(self, side, probe_width=0.015, rad_cutoff=0.1, plotnum=0):
        """
        Plot the 2D tungsten distribution across the face.

        side: Either 'ITF' or 'OTF'.
        probe_width: The half-width of the collector probe (the variable CPCO).
                     A = 0.015, B = 0.005, C = 0.0025
        rad_cutoff:  Only plot data from the tip down to rad_cutoff. Useful
                     if we want to compare to LAMS since those scans only go
                     down a certain length of the probe.

        *** To-Do ***
        - Instead of entering the width, pull out CPCO(?) from the netcdf file.
            Need to figure out the points being deposited outside the expected
            probe width first though.
        - Print out the ITF/OTF ratio from this analysis.
        """

        #The deposition array.
        dep_arr = self.get_dep_array()

        # Location of each P bin, and its width. Currently they all have the same width,
        # but it may end up such that there are custom widths so we leave it like this.
        ps     = np.array(self.nc.variables['PS'][:].data)
        pwids  = np.array(self.nc.variables['PWIDS'][:].data)

        # Array of poloidal locations (i.e. the center of each P bin).
        pol_locs = ps - pwids/2.0

        # Distance cell centers along surface (i.e. the radial locations).
        rad_locs = np.array(self.nc.variables['ODOUTS'][:].data)

        # Remove data beyond rad_cutoff.
        idx = np.where(np.abs(rad_locs)<rad_cutoff)[0]
        rad_locs = rad_locs[idx]
        dep_arr = dep_arr[:, idx]

        # Get only positive values of rad_locs for ITF...
        idx = np.where(rad_locs > 0.0)[0]
        X_itf, Y_itf = np.meshgrid(rad_locs[idx], pol_locs)
        Z_itf = dep_arr[:, idx]

        # ... negative for OTF.
        idx = np.where(rad_locs < 0.0)[0]
        X_otf, Y_otf = np.meshgrid(np.abs(rad_locs[idx][::-1]), pol_locs)
        Z_otf = dep_arr[:, idx][:, ::-1]

        # Make the levels for the contour plot out of whichever side has the max deposition.
        if Z_itf.max() > Z_otf.max():
            levels = np.linspace(0, Z_itf.max(), 15)
        else:
            levels = np.linspace(0, Z_otf.max(), 15)

        # Plotting commands.
        if side == 'ITF':
            X = X_itf; Y = Y_itf; Z = Z_itf
        else:
            X = X_otf; Y = Y_otf; Z = Z_otf

        if plotnum ==0:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            ax = self.master_fig.axes[plotnum-1]

        ax.contourf(X*100, Y*100, Z, levels=levels, cmap='Reds')
        ax.set_xlabel('Distance along probe (cm)', fontsize=fontsize)
        ax.set_ylabel('Z location (cm)', fontsize=fontsize)
        ax.set_ylim([-probe_width*100, probe_width*100])
        props = dict(facecolor='white')
        ax.text(0.75, 0.85, side, bbox=props, fontsize=fontsize*1.5, transform=ax.transAxes)

        if plotnum ==0:
            fig.tight_layout()
            fig.show()

        print('Total ITF/OTF (0-{} cm): {:.2f}'.format(rad_cutoff*100, Z_itf.sum()/Z_otf.sum()))

    def avg_pol_profiles(self, probe_width=0.015, rad_cutoff=0.5, plotnum=0):
        """
        Plot the average poloidal profiles for each side. Mainly to see if
        deposition peaks on the edges.

        probe_width: The half-width of the collector probe (the variable CPCO).
                     A = 0.015, B = 0.005, C = 0.0025
        rad_cutoff:  Only plot data from the tip down to rad_cutoff. Useful
                     if we want to compare to LAMS since those scans only go
                     down a certain length of the probe.
        """

        # Code copied from above function, deposition_contour. See for comments.
        dep_arr = np.array(self.nc.variables['NERODS3'][0] * -1)
        ps     = np.array(self.nc.variables['PS'][:].data)
        pwids  = np.array(self.nc.variables['PWIDS'][:].data)
        pol_locs = ps - pwids/2.0
        dep_arr = dep_arr[:-1, :]
        pol_locs = pol_locs[:-1]
        rad_locs = np.array(self.nc.variables['ODOUTS'][:].data)
        idx = np.where(np.abs(rad_locs)<rad_cutoff)[0]
        rad_locs = rad_locs[idx]
        dep_arr = dep_arr[:, idx]
        idx = np.where(rad_locs > 0.0)[0]
        X_itf, Y_itf = np.meshgrid(rad_locs[idx], pol_locs)
        Z_itf = dep_arr[:, idx]
        idx = np.where(rad_locs < 0.0)[0]
        X_otf, Y_otf = np.meshgrid(np.abs(rad_locs[idx][::-1]), pol_locs)
        Z_otf = dep_arr[:, idx][:, ::-1]

        # Average along the radial direction.
        avg_pol_itf = np.mean(Z_itf, 1)
        avg_pol_otf = np.mean(Z_otf, 1)

        # Get the centerline index (or closest to it).
        cline = np.abs(pol_locs).min()
        cline_idx = np.where(pol_locs == cline)[0][0]

        # Get average peaking factor for each side.
        peak1 = avg_pol_itf[:cline_idx].max() / avg_pol_itf[cline_idx]
        peak2 = avg_pol_itf[cline_idx:].max() / avg_pol_itf[cline_idx]
        itf_peak = (peak1 + peak2) / 2.0
        peak1 = avg_pol_otf[:cline_idx].max() / avg_pol_otf[cline_idx]
        peak2 = avg_pol_otf[cline_idx:].max() / avg_pol_otf[cline_idx]
        otf_peak = (peak1 + peak2) / 2.0

        # Plotting commands.
        if plotnum ==0:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            ax = self.master_fig.axes[plotnum-1]

        ax.plot(pol_locs, avg_pol_itf/avg_pol_itf.max(), label='ITF', color=tableau20[6])
        ax.plot(pol_locs, avg_pol_otf/avg_pol_otf.max(), label='OTF', color=tableau20[8])
        ax.legend(fontsize=fontsize)
        ax.set_xlabel('Poloidal (m)', fontsize=fontsize)
        ax.set_ylabel('Deposition (normalized)', fontsize=fontsize)
        ax.set_xlim([-probe_width, probe_width])

        if plotnum==0:
            fig.tight_layout()
            fig.show()

        # Print and then return the message for the GUI to use.
        message = "OTF/ITF Peaking Ratio: {:.2f}".format(otf_peak/itf_peak)
        print(message)
        return message

    def overviewplot(self):

        self.master_fig = plt.figure(figsize=(12,6))
        for x in range(1, 10):
            self.master_fig.add_subplot(3, 3, x)

        self.centerline(plotnum=1)
        self.deposition_contour(side='ITF', plotnum=2)
        self.deposition_contour(side='OTF', plotnum=3)

        self.master_fig.tight_layout()
        self.master_fig.show()

    def multiplot_start(self):

        self.master_fig = plt.figure(figsize=(12, 6))
        for x in range(1, 10):
            self.master_fig.add_subplot(3, 3, x)

    def multiplot_end(self):

        self.master_fig.tight_layout()
        self.master_fig.show()

