import numpy             as np
import matplotlib.pyplot as plt
import matplotlib        as mpl
import netCDF4
from matplotlib import colors
from collections import OrderedDict


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


class Readout:
        """
        Class object that holds a figure of a grid of plots from the netCDF
        output of 3DLIM. Example usage in controlling script lim_readout.py.
        """

        def __init__(self, netcdf_file=None, figsize=(15,10), grid_shape=(3,3)):
            """
            netcdf_file: Path to the 3DLIM netCDF file. If None is entered then
                         it will use 'colprobe-test-m2.nc' as a default, which
                         is convienent for testing.
            figsize:     Size of the figure to hold all the plots. The default of
                         (15, 10) is a good size.
            grid_shape:  The shape of the grid of plots. Change if you want to
                         add more plots or whatever.
            """

            # Create the master figure to hold all the plots.
            self.master_fig = plt.figure(figsize=figsize)

            # If no netCDF file given, just use this test one.
            if netcdf_file is None:
                self.netcdf = netCDF4.Dataset('colprobe-test-m2.nc')
            else:
                self.netcdf = netCDF4.Dataset(netcdf_file)


            # Create figure with array of empty plots.
            for plot_num in range(1, grid_shape[0] * grid_shape[1] + 1):
                self.master_fig.add_subplot(grid_shape[0], grid_shape[1], plot_num)


        def print_readout(self):
            """
            Output a table with relevant info from the netcdf file.
            """

            # Let's just put everything we want into a dict so printing is easy.
            output = OrderedDict()
            output['3DLIM Version'] = self.netcdf['VERSION'][:].data.tostring().decode()
            output['Title'] = self.netcdf['TITLE'][:].data.tostring().decode()
            output['File'] = self.netcdf['JOB'][:].data.tostring().decode().split(' ')[0]
            output['Particles'] = format(self.netcdf['MAXIMP'][:].data, ',')
            output['Conn. Length'] = self.netcdf['CL'][:].data

            # Find longest output for formatting.
            pad = 0
            for val in output.values():
                if len(str(val)) > pad:
                    pad = len(str(val))

            # Printing commands.
            num_stars = 2 + 15 + 2 + pad
            print("\n" + "*"*num_stars)
            for key, val in output.items():
                print("* {:15}{:<{pad}} *".format(key, val, pad=pad))
            print("*"*num_stars)


        def centerline(self, plot_num):
            """
            Plot the ITF and OTF deposition along the centerlines.

            plot_num: Location in grid to place this plot. I.e. if the grid_shape
                      is (3,3), then enter a number between 0-8, where the locations
                      are labelled left to right.
            """

            #The deposition array.
            dep_arr = np.array(self.netcdf.variables['NERODS3'][0] * -1)

            # Location of each P bin, and its width. Currently they all have the same width,
            # but it may end up such that there are custom widths so we leave it like this.
            ps     = np.array(self.netcdf.variables['PS'][:].data)
            pwids  = np.array(self.netcdf.variables['PWIDS'][:].data)

            # Array of poloidal locations (i.e. the center of each P bin).
            pol_locs = ps - pwids/2.0

            # Drop last row since it's garbage.
            dep_arr = dep_arr[:-1, :]
            pol_locs = pol_locs[:-1]

            # Distance cell centers along surface (i.e. the radial locations).
            rad_locs = np.array(self.netcdf.variables['ODOUTS'][:].data)

            # Get the centerline index (or closest to it).
            cline = np.abs(pol_locs).min()

            # Index the deposition array at the centerline for plotting.
            itf_x = rad_locs[np.where(rad_locs > 0.0)[0]]
            itf_y = dep_arr[np.where(pol_locs == cline)[0], np.where(rad_locs > 0.0)[0]]
            otf_x = rad_locs[np.where(rad_locs < 0.0)[0]] * -1
            otf_y = dep_arr[np.where(pol_locs == cline)[0], np.where(rad_locs < 0.0)[0]]

            # Plotting commands.
            ax = self.master_fig.axes[plot_num]
            ax.plot(itf_x, itf_y, '-', label='ITF', ms=ms, color=tableau20[6])
            ax.plot(otf_x, otf_y, '-', label='OTF', ms=ms, color=tableau20[8])
            ax.legend(fontsize=fontsize)
            ax.set_xlabel('Distance along probe (m)', fontsize=fontsize)
            ax.set_ylabel('Deposition (arbitrary units)', fontsize=fontsize)


        def deposition_contour(self, plot_num, side, probe_width=0.015, rad_cutoff=0.1):
            """
            Plot the 2D tungsten distribution across the face.

            plot_num:    Location in grid to place this plot. I.e. if the grid_shape
                         is (3,3), then enter a number between 0-8, where the locations
                         are labelled left to right.
            side: Either 'ITF' or 'OTF'.
            probe_width: The half-width of the collector probe (the variable CPCO).
                         A = 0.015, B = 0.005, C = 0.0025
            rad_cutoff:  Only plot data from the tip down to rad_cutoff. Useful
                         if we want to compare to LAMS since those scans only go
                         down a certain length of the probe.

            """

            #The deposition array.
            dep_arr = np.array(self.netcdf.variables['NERODS3'][0] * -1)

            # Location of each P bin, and its width. Currently they all have the same width,
            # but it may end up such that there are custom widths so we leave it like this.
            ps     = np.array(self.netcdf.variables['PS'][:].data)
            pwids  = np.array(self.netcdf.variables['PWIDS'][:].data)

            # Array of poloidal locations (i.e. the center of each P bin).
            pol_locs = ps - pwids/2.0

            # Drop last row since it's garbage.
            dep_arr = dep_arr[:-1, :]
            pol_locs = pol_locs[:-1]

            # Distance cell centers along surface (i.e. the radial locations).
            rad_locs = np.array(self.netcdf.variables['ODOUTS'][:].data)

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

            ax = self.master_fig.axes[plot_num]
            ax.contourf(X, Y, Z, levels=levels, cmap='Reds')
            ax.set_xlabel('Distance along probe (m)', fontsize=fontsize)
            ax.set_ylabel('Z location (m)', fontsize=fontsize)
            ax.set_ylim([-probe_width, probe_width])
            props = dict(facecolor='white')
            ax.text(0.75, 0.85, side, bbox=props, fontsize=fontsize*1.5, transform=ax.transAxes)


        def velocity_contour_pol(self, pol_slice=0):
            """
            Plot the 2D distribution of the (tungsten? plasma?) velocity at a
            poloidal slice.

            pol_slice: The poloidal coordinate to get a velocity plot in (R, B) space.
            """
            pass

        def velocity_contour_par(self, par_slice=0):
            """
            Plot the 2D distribution of the (tungsten? plasma?) velocity at a
            parallel (to B) slice.

            par_slice: The parallel coordinate to get a velocity plot in (R, P) space.
            """
            pass

        def te_plot(self):
            """
            Plot the input Te (which is at the midplane?).
            """
            pass

        def ne_plot(self):
            """
            Plot the input ne (which is at the midplane?).
            """
            pass

        def te_contour(self, plot_num):
            """
            Plot the 2D background electron plasma temperature.

            plot_num: Location in grid to place this plot. I.e. if the grid_shape
                      is (3,3), then enter a number between 0-8, where the locations
                      are labelled left to right.
            """

            # Get the connection length to restrict the plot between the two absorbing surfaces.
            cl = float(self.netcdf['CL'][:].data)

            # Same with the location of the plasma center (the top of the box).
            ca = float(self.netcdf['CA'][:].data)

            # Get the X and Y grid data.
            x = self.netcdf.variables['XOUTS'][:].data
            y = self.netcdf.variables['YOUTS'][:].data

            # 2D grid of the temperature data.
            Z = self.netcdf.variables['CTEMBS'][:].data

            # Trim the zeros from the edges of the x and y arrays, and the associated
            # data points as well. This is done to stop this data from messing up
            # the contours in the contour plot.
            xkeep_min = np.nonzero(x)[0].min()
            xkeep_max = np.nonzero(x)[0].max()
            ykeep_min = np.nonzero(y)[0].min()
            ykeep_max = np.nonzero(y)[0].max()
            x = x[xkeep_min:xkeep_max]
            y = y[ykeep_min:ykeep_max]
            Z = Z[ykeep_min:ykeep_max, xkeep_min:xkeep_max]

            # Furthermore, trim the data off that is beyond CL.
            ykeep_cl = np.where(np.abs(y) < cl)[0]
            y = y[ykeep_cl]
            Z = Z[ykeep_cl, :]

            # Replace zeros in Z with just the smallest density value. Again to
            # stop all these zeros from messing up the contour levels.
            Zmin = np.partition(np.unique(Z), 1)[1]
            Z = np.clip(Z, Zmin, None)

            # Create grid for plotting. Note we swap definitions for x and y since
            # we want the x-axis in the plot to be the parallel direction (it just
            # looks better that way).
            Y, X = np.meshgrid(x, y)

            # Plotting commands.
            ax = self.master_fig.axes[plot_num]
            cont = ax.contourf(X, Y, Z, cmap='magma', levels=10)
            ax.set_xlim([-cl, cl])
            #ax.set_ylim([None, ca])
            ax.set_ylim([None, 0.01])  # Contour weird near edge.
            ax.set_xlabel('Parallel (m)')
            ax.set_ylabel('Radial (m)')
            cbar = self.master_fig.colorbar(cont, ax=ax)
            cbar.set_label('Background Te (eV)')

        def ne_contour(self, plot_num):
            """
            Plot the 2D background plasma density.

            plot_num: Location in grid to place this plot. I.e. if the grid_shape
                      is (3,3), then enter a number between 0-8, where the locations
                      are labelled left to right.
            """
            # Get the connection length to restrict the plot between the two absorbing surfaces.
            cl = float(self.netcdf['CL'][:].data)
            # Same with the location of the plasma center (the top of the box)
            ca = float(self.netcdf['CA'][:].data)

            # Get the X and Y grid data.
            x = self.netcdf.variables['XOUTS'][:].data
            y = self.netcdf.variables['YOUTS'][:].data

            # 2D grid of the temperature data.
            Z = self.netcdf.variables['CRNBS'][:].data

            # Trim the zeros from the edges of the x and y arrays, and the associated
            # data points as well. This is done to stop this data from messing up
            # the contours in the contour plot.
            xkeep_min = np.nonzero(x)[0].min()
            xkeep_max = np.nonzero(x)[0].max()
            ykeep_min = np.nonzero(y)[0].min()
            ykeep_max = np.nonzero(y)[0].max()
            x = x[xkeep_min:xkeep_max]
            y = y[ykeep_min:ykeep_max]
            Z = Z[ykeep_min:ykeep_max, xkeep_min:xkeep_max]

            # Furthermore, trim the data off that is beyond CL.
            ykeep_cl = np.where(np.abs(y) < cl)[0]
            y = y[ykeep_cl]
            Z = Z[ykeep_cl, :]

            # Replace zeros in Z with just the smallest density value. Again to
            # stop all these zeros from messing up the contour levels.
            Zmin = np.partition(np.unique(Z), 1)[1]
            Z = np.clip(Z, Zmin, None)

            # Create grid for plotting. Note we swap definitions for x and y since
            # we want the x-axis in the plot to be the parallel direction (it just
            # looks better that way).
            Y, X = np.meshgrid(x, y)

            # Plotting commands.
            ax = self.master_fig.axes[plot_num]

            # Create our own levels since the automatic ones are bad.
            lev_exp = np.arange(np.floor(np.log10(Z.min())-1), np.ceil(np.log10(Z.max())+1), 0.25)
            levs = np.power(10, lev_exp)

            cont = ax.contourf(X, Y, Z, cmap='magma', levels=levs, norm=colors.LogNorm())
            ax.set_xlim([-cl, cl])
            #ax.set_ylim([None, ca])
            ax.set_ylim([None, 0.01])  # Contour weird near edge.
            ax.set_xlabel('Parallel (m)')
            ax.set_ylabel('Radial (m)')
            cbar = self.master_fig.colorbar(cont, ax=ax)
            cbar.set_label('Background ne (m-3)')

        def avg_imp_vely(self, plot_num):
            """
            SVYBAR: Average impurity velocity at X coordinates in QXS.

            plot_num: Location in grid to place this plot. I.e. if the grid_shape
                      is (3,3), then enter a number between 0-8, where the locations
                      are labelled left to right.
            """

            # Grab the data.
            x = self.netcdf.variables['QXS'][:].data
            y = self.netcdf.variables['SVYBAR'][:].data

            # Plotting commands.
            ax = self.master_fig.axes[plot_num]
            ax.plot(x, y, '.', ms=ms, color=tableau20[6])
            ax.set_xlabel('Radial coordinates (m)', fontsize=fontsize)
            ax.set_ylabel('Average Y imp. vel. (m/s)', fontsize=fontsize)

        def show_fig(self):
            """
            Accessor for showing the master_fig and cleaning up the layout.
            """
            self.master_fig.tight_layout()
            self.master_fig.show()
