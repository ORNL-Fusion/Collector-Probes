# Author: Shawn Zamperini
#
# This script was more or less copied from the oedge_omfit script by J. Nichols.
# It's been rewritten because it's a good exercise and has been tweaked to
# better fit in future plans for a GUI interface.

import netCDF4
import numpy             as np
import matplotlib        as mpl
import matplotlib.pyplot as plt
import pandas            as pd
from collections         import OrderedDict


# A nice looking font.
plt.rcParams['font.family'] = 'serif'

# Create class object to hold in netcdf data as well as plotting routines. Use
# Ordered Dict prototype so we can index the class object with the data name.
# Ex: self['BTS']
class OedgePlots:

    def __init__(self, netcdf_path):
        """
        Initialization for class object. Loads in the relevant netCDF variables
        needed for the plots.

        netcdf_path: Path location of .nc file.
        """

        # Load in the netCDF file.
        self.nc = netCDF4.Dataset(netcdf_path)

        # Initialize variable to hold dat file.
        self.dat_file = None

        # Load in some netCDF data that is used a lot.
        self.rs     = self.nc['RS'][:]
        self.zs     = self.nc['ZS'][:]
        self.nrs    = self.nc['NRS'][:]
        self.nks    = self.nc['NKS'][:]
        self.area   = self.nc['KAREAS'][:]
        self.korpg  = self.nc['KORPG'][:]
        self.rvertp = self.nc['RVERTP'][:]
        self.zvertp = self.nc['ZVERTP'][:]
        self.rvesm  = self.nc['RVESM'][:]
        self.zvesm  = self.nc['ZVESM'][:]
        self.irsep  = self.nc['IRSEP'][:]
        self.qtim   = self.nc['QTIM'][:]
        self.absfac = self.nc['ABSFAC'][:]
        self.kss    = self.nc['KSS'][:]
        self.kfizs  = self.nc['KFIZS'][:]

        # Create a mesh of of the corners of the each cell/polygon in the grid.
        #mesh  = np.array([])
        mesh = []
        num_cells = 0

        # Scan through the rings.
        for ir in range(self.nrs):

            # Scan through the knots.
            for ik in range(self.nks[ir]):

                # Get the cell index of this knot on this ring.
                index = self.korpg[ir,ik] - 1

                # Only if the area of this cell is not zero append the corners.
                if self.area[ir,ik] != 0.0:
                    vertices = list(zip(self.rvertp[index][0:4], self.zvertp[index][0:4]))
                    #mesh = np.append(mesh, vertices)
                    mesh.append(vertices)
                    num_cells = num_cells + 1

                    # Print out a warning is the cell center is not within the vertices.
                    cell = mpl.path.Path(list(vertices))
                    r = self.rs[ir, ik]
                    z = self.zs[ir, ik]
                    if not cell.contains_point([r, z]):
                        print("Error: Cell center not within vertices.")
                        print("  (ir, ik)    = ({}, {})".format(ir, ik))
                        print("  Vertices    = {}".format(vertices))
                        print("  Cell center = ({}, {})".format(r, z))

        # Save the results in the class.
        self.num_cells = num_cells
        self.mesh = mesh

    def add_dat_file(self, dat_path):
        """
        Quick function to read in the dat_file into the class. The file is not
        suitible for pandas or anything like that, so just read it in as a text
        with newlines at the end of each line.

        dat_path: Path location of .dat file.
        """

        with open(dat_path) as f:
            self.dat_file = f.read()

    def read_data_2d(self, dataname, charge=None, scaling=1.0):
        """
        Reads in 2D data into a 1D array, in a form that is then passed easily
        to PolyCollection for plotting.

        dataname : The 2D data as named in the netCDF file.
        charge   : The charge state to be plotted, if applicable.
        scaling  : Scaling factor to apply to the data, if applicable.
        """

        # Get the 2D data from the netCDF file.
        raw_data = self.nc[dataname][:]
        data = np.zeros(self.num_cells)

        # 'all' just sums up all charge states. So here instead of loop for speed.
        if charge == 'all':
            raw_data = raw_data.sum(axis=0)

        count = 0
        for ir in range(self.nrs):
            for ik in range(self.nks[ir]):
                if self.area[ir, ik] != 0.0:

                    # If charge is specifed, this will be the first dimension,
                    # and will need to charge + 1 to match index.
                    if charge in [None, 'all']:
                        data[count] = raw_data[ir][ik] * scaling
                    else:
                        data[count] = raw_data[charge + 1][ir][ik] * scaling

                    count = count + 1

        return data

    def read_data_2d_kvhsimp(self):
        """
        Special function for plotting the flow velocity of the impurities. This
        is because some DIVIMP options (T13, T31, T37?, T38?...) add
        additional flow values onto the impurities not reflected in KVHS. These
        additional values are in the .dat file, and thus it is required to run
        this function.
        """

        # Make sure .dat file has been loaded.
        if self.dat_file == None:
            print("Error: .dat file not loaded in. Run 'add_dat_file' first.")

        pol_opt = float(self.dat_file.split('POL DRIFT OPT')[1].split(':')[0])
        if pol_opt == 0.0:
            print("Error: Poloidal drift option T13 was not on for this run.")

        # Get the relevant table for the extra drifts out of the .dat file.
        add_data = self.dat_file.split('TABLE OF DRIFT REGION BY RING - RINGS ' + \
                                        'WITHOUT FLOW ARE NOT LISTED\n')[1]. \
                                        split('DRIFT')[0].split('\n')

        # Split the data between the spaces, put into DataFrame.
        add_data = [line.split() for line in add_data]
        add_df = pd.DataFrame(add_data[1:-1], columns=['IR', 'Vdrift (m/s)',
                              'S_START (m)', 'S_END (m)'], dtype=np.float64). \
                              set_index('IR')

        # Get the 2D data from the netCDF file.
        dataname = 'KVHS'
        scaling = 1.0 / self.qtim
        raw_data = self.nc[dataname][:]
        data = np.zeros(self.num_cells)

        # Convert the 2D data (ir, ik) into 1D for plotting in the PolyCollection
        # matplotlib function.
        count = 0
        for ir in range(self.nrs):
            for ik in range(self.nks[ir]):
                if self.area[ir, ik] != 0.0:

                    # Put the data from this [ring, knot] into a 1D array.
                    data[count] = raw_data[ir][ik] * scaling

                    # If this ring has additional drifts to be added.
                    if ir in add_df.index:

                        # Then add the drift along the appropriate s (or knot) range.
                        if self.kss[ir][ik] > add_df['S_START (m)'].loc[ir] and \
                           self.kss[ir][ik] < add_df['S_END (m)'].loc[ir]:

                           data[count] = data[count] + add_df['Vdrift (m/s)'].loc[ir]

                    count = count + 1

        return data

    def get_sep(self):
        """
        Return collection of lines to be plotted with LineCollection method
        of matplotlib for the separatrix.
        """

        # Get (R, Z) coordinates of separatrix.
        rsep = self.rvertp[self.korpg[self.irsep-1,:self.nks[self.irsep-1]]][:,0]
        zsep = self.zvertp[self.korpg[self.irsep-1,:self.nks[self.irsep-1]]][:,0]
        nsep = len(rsep)
        lines=[]

        # Construct separatrix as a series of pairs of coordinates (i.e. a line
        # between the coordinates), to be plotted. Don't connect final point to first.
        for i in range(nsep-2):
            lines.append([(rsep[i], zsep[i]), (rsep[i+1], zsep[i+1])])

        return lines

    def plot_contour_polygon(self, dataname, charge=None, scaling=1.0,
                             normtype='linear', cmap='plasma', xlim=[0.9, 2.5],
                             ylim = [-1.5, 1.5], plot_sep=True, levels=None,
                             cbar_label=None, fontsize=16, lut=21,
                             smooth_cmap=False, vmin=None, vmax=None,
                             show_cp=None, ptip=None, show_mr=False):

        """
        Create a standalone figure using the PolyCollection object of matplotlib.

        dataname:    The netCDF variable name to plot. Some special datanames
                       will perform extra data handling: KVHSimp, ...
        charge:      The charge state to be plotted, if applicable.
        scaling:     Scaling factor to apply to the data, if applicable.
        normtype:    One of 'linear', 'log', ... of how to normalize the data on
                       the plot.
        cmap:        The colormap to apply to the plot. Uses standard matplotlib
                       names.
        xlim:        X range of axes.
        ylim:        Y range of axes.
        plot_sep:    Include separatrix on plot or not.
        levels:      Number of levels for colorbar (needs work).
        cbar_label:  Label for the colorbar.
        fontsize:    Size of font for labels.
        lut:         Number of chunks to break the colormap into.
        smooth_cmap: Choose whether to break colormap up into chunks or not.
        vmin/vmax:   Option to choose own vmin, vmax for the colorbars. Useful
                       when plots need tweaking.
        show_cp:     Choose if collector probes are to be shown. Can have
                       multiple probes as a list. 1 = MiMES, 2 = Top, 3 = DiMES.
        ptip:        The tip (either R or Z, depends on which probe you picked)
                       of the probe. Can pass as list where each entry corresponds
                       to the one in show_cp.
        show_mr:     Option to show the metal rings or not from MRC-I.
        """

        # Make sure show_cp and ptip is in list form if not.
        #if show_cp != None:
        if type(show_cp) is not list:
            show_cp = [show_cp]
        if type(ptip) is not list:
            ptip = [ptip]

        # Read in the data into a form for PolyCollection.
        if dataname == 'KVHSimp':
            data = self.read_data_2d_kvhsimp()
        else:
            data = self.read_data_2d(dataname, charge, scaling)


        # Create a good sized figure with correct proportions.
        fig = plt.figure(figsize=(7, 9))
        ax  = fig.add_subplot(111)

        if normtype == 'linear':
            if vmin == None:
                vmin = data.min()
            if vmax == None:
                vmax = data.max()
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

        elif normtype == 'log':
            data[data == 0.0] = 1e-3
            if vmin == None:
                vmin = data.min()
            if vmax == None:
                vmax = data.max()
            norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)

        elif normtype == 'symlin':
            if vmin == None:
                vmin = -np.abs(data).max()
            if vmax == None:
                vmax = -vmin
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
            cmap = 'coolwarm'

        elif normtype == 'symlog':
            data[data == 0.0] = 1e-3
            if vmin == None:
                vmin = -np.abs(data).max()
            if vmax == None:
                vmax = -vmin
            norm = mpl.colors.SymLogNorm(linthresh=0.01 * vmax, vmin=vmin, vmax=vmax)
            cmap = 'coolwarm'

        # Create the PolyCollection object. Choose whether to discretize the
        # colormap or not first.
        if smooth_cmap:
            scalar_map = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        else:
            scalar_map = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.get_cmap(cmap, lut=lut))
        coll = mpl.collections.PolyCollection(self.mesh, array=data,
                                              cmap=scalar_map.cmap,
                                              norm=scalar_map.norm,
                                              edgecolors='none')

        # Add the PolyCollection to the Axes object.
        ax.add_collection(coll)
        ax.plot(self.rvesm, self.zvesm, color='k', linewidth=1)

        # Get the separatrix coordinates as a collection of lines and plot.
        if plot_sep:
            sep = self.get_sep()
            sc = mpl.collections.LineCollection(sep, color='k')
            ax.add_collection(sc)

        # Use correct amount of levels for colorbar, if specified.
        if levels is not None:
            cbar = fig.colorbar(coll, ax=ax, boundaries=levels, ticks=levels)
        else:
            cbar = fig.colorbar(coll, ax=ax, extend='both')

        # Add colorbar label.
        if cbar_label is None:
            cbar.ax.set_ylabel(dataname, fontsize=fontsize)
        else:
            cbar.ax.set_ylabel(cbar_label, fontsize=fontsize)

        # Option to add collector probes to plots.
        if show_cp:
            if ptip == None:
                print("Error: Location of tip of probe not given.")

        # MiMES probe.
        cp_width  = 0.03
        facecolor = 'grey'
        edgecolor = 'black'
        if 1 in show_cp:
            ptip_tmp = ptip[show_cp.index(1)]
            lower_xy = (ptip_tmp, -0.185 - cp_width/2.0)
            width    = 2.38 - ptip_tmp
            height   = cp_width
            rect = mpl.patches.Rectangle(lower_xy, width=width, height=height,
                                         facecolor=facecolor, edgecolor=edgecolor)
            ax.add_patch(rect)

        # Top crown probe. Point is to simulate a DiMES probe if it were USN.
        if 2 in show_cp:
            ptip_tmp = ptip[show_cp.index(2)]
            lower_xy = (1.485 - cp_width/2.0, ptip_tmp)
            width    = cp_width
            height   = ptip_tmp - 1.18
            rect = mpl.patches.Rectangle(lower_xy, width=width, height=-height,
                                         facecolor=facecolor, edgecolor=edgecolor)
            ax.add_patch(rect)

        # DiMES probe.
        if 3 in show_cp:
            ptip_tmp = ptip[show_cp.index(3)]
            lower_xy = (1.485 - cp_width/2.0, -1.25)
            width    = cp_width
            height   = ptip_tmp - (-1.25)
            rect = mpl.patches.Rectangle(lower_xy, width=width, height=height,
                                         facecolor=facecolor, edgecolor=edgecolor)
            ax.add_patch(rect)

        # Option to show metal rings. Floor (1.32-1.37) Shelf (1.404, 1.454).
        if show_mr:
            tile_height = 0.01
            facecolor   = 'red'
            edgecolor   = 'black'
            floor_xy    = (1.32, -1.363 - tile_height)
            floor_width = 1.37 - 1.32
            shelf_xy    = (1.404, -1.250 - tile_height)
            shelf_width = 1.454 - 1.404
            floor_rect  = mpl.patches.Rectangle(floor_xy, width = floor_width,
                                               height=tile_height,
                                               facecolor=facecolor,
                                               edgecolor=edgecolor)


            shelf_rect  = mpl.patches.Rectangle(shelf_xy, width = shelf_width,
                                               height=tile_height,
                                               facecolor=facecolor,
                                               edgecolor=edgecolor)
            ax.add_patch(floor_rect)
            ax.add_patch(shelf_rect)

        # Organize plot, add labels.
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel('R (m)', fontsize=fontsize)
        ax.set_ylabel('Z (m)', fontsize=fontsize)
        fig.tight_layout()
        fig.show()

        return fig

    def cp_plots(self, cp_path, xaxis='ROMP', yaxis='IMPFLUX', cp_num=1,
                 fontsize=16, log=False, print_to_file=False):
        """
        Plot some of the collector probe data. Note: This still requires work,
        namely with the crown (top) probe.

        cp_path:      Path to the .collector_probe file.
        xaxis:        Which column name to plot on the xaxis.
        yaxis:        Which column name to plot on the yaxis, excluding the _IDF, _ODF.
                         One of IMPFLUX or IMPDENS
        cp_num:        1 = midplane probe, 2 = crown probe.
        fontsize:      Font size for the plots.
        log:           Set to true to use a log y axis.
        print_to_file: Print output to file in directory of netCDF file for
                         external use.
        """

        # Open files and read all the lines in.
        with open(cp_path) as f:
            lines = np.array(f.readlines())

        # If there is another cp on this file it will start after two \n's.
        for i in range(len(lines) - 1):
            if lines[i] == '\n':
                if lines[i+1] == '\n':
                    next_cp_idx = i

        # Grab the ABSFAC for the scaling.
        absfac = np.float(lines[6].split('     ')[1])

        # Load first cp as DataFrame.
        df1 = pd.read_csv(cp_path, skiprows=7, nrows=next_cp_idx-8, sep=' ',
                         skipinitialspace=True)

        # The INDEX column on the file doesn't have numbers under it so it goofs
        # up he df with an extra NaN column. Fix real quick.
        correct_cols = df1.columns[1:]
        df1 = df1.drop(df1.columns[-1], axis=1)
        df1.columns = correct_cols

        # Read in the second probe if there is one.
        df2 = pd.read_csv(cp_path, skiprows=8+next_cp_idx, sep=' ',
                          skipinitialspace=True)
        correct_cols = df2.columns[1:]
        df2 = df2.drop(df2.columns[-1], axis=1)
        df2.columns = correct_cols

        # The midplane probe I think.
        if cp_num == 1:
            x = df1[xaxis].values
            y_itf = df1[yaxis + '_IDF'].values
            y_otf = df1[yaxis + '_ODF'].values

        # The crown probe I think.
        elif cp_num == 2:
            x = df2[xaxis].values
            y_itf = df2[yaxis + '_IDF'].values
            y_otf = df2[yaxis + '_ODF'].values

        if xaxis == 'ROMP':
            xlabel = 'R-Rsep OMP (cm)'
            x = x  * 100.0

        if yaxis == 'IMPFLUX':
            ylabel = 'Impurity Flux (m-2 s-1)'
            y_itf = y_itf * absfac
            y_otf = y_otf * absfac

        # Plotting commands.
        red  = (214/255, 39/255, 40/255)
        purp = (148/255, 103/255, 189/255)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if log:
            ax.semilogy(x, y_itf, label='ITF', lw=5, color=red)
            ax.semilogy(x, y_otf, label='OTF', lw=5, color=purp)
        else:
            ax.plot(x, y_itf, label='ITF', lw=5, color=red)
            ax.plot(x, y_otf, label='OTF', lw=5, color=purp)
        ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.legend(fontsize=fontsize)
        fig.tight_layout()
        fig.show()

        if print_to_file:
            import csv
            filename = cp_path.split('.')[0] + '.txt'
            with open(filename, 'w') as f:
                writer = csv.writer(f, delimiter='\t')
                f.write(xaxis+'\tITF\tOTF\n')
                writer.writerows(zip(x, y_itf, y_otf))
