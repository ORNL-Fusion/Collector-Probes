import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib as mpl
import sys
from collections import defaultdict


# Few constants.
e = 1.609e-19
md = 2 * 931.494 * 1e6 / (3e8)**2  # Mass of Dueterium in eV/m2/s.

# These are the "Tableau 20" colors as RGB. The last one is just black. I added it.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229),
             (0, 0, 0)]

# A nice looking font.
#plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.family'] = 'DejaVu Sans'

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

class SolpsPlots:
    """

    """

    def __init__(self, folder_path=None, verbal=False):

        # Shortcut to loading in the data on init.
        if folder_path != None:
            self.load_b2f_files(folder_path, verbal)
        self.folder_path = folder_path

    def __repr__(self):

        repr_str =  "SolpsPlots object\n"
        repr_str += "  Folder:   {}\n".format(self.folder_path)
        repr_str += "  Label:    {}\n".format(self.label)
        try:
            date = self.label.split()[4]
            user = self.label.split()[6]
        except:
            date = ""
            user = ""
        repr_str += "  Run Date: {}\n".format(date)
        repr_str += "  Run By:   {}\n".format(user)

        return repr_str

    def read_file(self, file_path, verbal=False):
        """
        Generic function for loading in data from the b2fgmtry, b2fplasmf and
        b2fstate output files.

        Inputs
        file_path (str): Path to the file we're reading in.
        verbal (bool):   If you want output to see what's being read.

        Outputs
        var_dict (dict): The output dictionary with all the read in data.
        """

        count = 0

        # Use defaultdict with lists since they're quick to append to.
        var_dict = defaultdict(list)
        with open(file_path) as f:
            while True:
                count += 1
                l = f.readline()

                # A *cf indicates this is a header telling us the vairable name
                #  and how many entries there are.
                if l[:3] == '*cf':

                    # If it's a char then there's just one line beneath it.
                    if l.split()[1] == 'char':
                        lines_to_read = 1
                    else:

                        # When it's a real, there are 6 entries per row. For ints
                        # there are 12.
                        if l.split()[1] == 'real':
                            num_per_line = 6
                        elif l.split()[1] == 'int':
                            num_per_line = 12
                        else:
                            print("Error: Do not know how to handle type '{}'.'".format(l.split()[1]))

                        # Calculate how many lines we will need to read from the file.
                        lines_to_read = int(np.ceil(int(l.split()[2]) / num_per_line))
                    var_name = l.split()[3]
                    if verbal:
                        print("{} {} {} {}".format(count, var_name, l.split()[1], lines_to_read))
                        #print("{}: ".format(var_name))

                    # Read in the specified number of lines, adding each set of
                    # values onto the list in out dictionary.
                    for i in range(0, lines_to_read):
                        count += 1
                        k = f.readline()
                        var_dict[var_name].extend(k.split())

                # An empty string means we're at EOF stuff, so we're done.
                if l == '':
                    break

            return var_dict

    def load_gmtry(self, gmtry_path, verbal=False):
        """
        Load the b2fgmtry file into the class object, storing each entry
        individually for easy access. See var_names.txt for descriptions
        of some of the output variables.

        Inputs
        gmtry_path (str): Path to the b2fgmtry file.
        verbal (bool):    If you want output to see what's being read.

        Outputs:
        None
        """

        # Use the read_file function to get everything from the gmtry file.
        gmtry_dict = self.read_file(gmtry_path, verbal)

        # Store all the goodies individually in the class for easy access. Do
        # some special handling on specific variables.
        for key, arr in gmtry_dict.items():
            if key == 'nx,ny':
                self.nx = int(arr[0])
                self.ny = int(arr[1])
            elif key == 'label':
                self.label = " ".join(arr)

            # Values that are just single ints store as such.
            elif key in ['isymm', 'nlreg', 'nlxlo', 'nlxhi', 'nlylo', 'nlyhi',
                         'nlloc', 'nncut', 'leftcut', 'rightcut', 'topcut',
                         'bottomcut', 'periodic_bc', 'reder_pbs']:
                setattr(self, key, int(arr[0]))

            # Normal lists, just store as numpy arrays.
            else:
                setattr(self, key, np.array(arr, dtype=np.float64))

    def load_plasmf(self, plasmf_path, verbal=False):
        """
        Load variables from the b2fplasmf file into the class. See var_names.txt
        for descriptions of some of the variables in this file.

        Inputs
        plasmf_path (str): Path to the b2fplasmf file.
        verbal (bool):    If you want output to see what's being read.

        Outputs:
        None
        """

        # Use the read_file function to get everything from the b2fplasmf file.
        plasmf_dict = self.read_file(plasmf_path, verbal)

        # Store the goodies in out class object for easy access.
        for key, arr in plasmf_dict.items():
            setattr(self, key, np.array(arr, dtype=np.float64))

    def load_state(self, state_path, verbal=False):
        """
        Load variables from the b2fstate file into the class. See var_names.txt
        for descriptions of some of the variables in this file.

        Inputs
        state_path (str): Path to the b2fstate file.
        verbal (bool):    If you want output to see what's being read.

        Outputs:
        None
        """

        # Use the read_file function to get everything from the b2fstate file.
        state_dict = self.read_file(state_path, verbal)

        # Store the goodies in the class object for easy access.
        for key, arr in state_dict.items():
            if key == 'nx,ny,ns':
                self.nx = int(arr[0])
                self.ny = int(arr[1])
                self.ns = int(arr[2])
            elif key == 'label':
                self.label = " ".join(arr)
            else:
                setattr(self, key, np.array(arr, dtype=np.float64))

    def load_b2f_files(self, folder_path, verbal=False):
        """
        Wrapper to just run the three functions for loading the needed files.

        Inputs
        folder_path (str): Path to the folder that holds all the files from the run.
        verbal (bool): Show readout of what's being read.

        Outputs
        None
        """
        try:
            self.load_gmtry(folder_path + 'b2fgmtry', verbal)
        except FileNotFoundError:
            print("Error: b2fgmtry not found.")
        try:
            self.load_plasmf(folder_path + 'b2fplasmf', verbal)
        except FileNotFoundError:
            print("Error: b2fplasmf not found.")
        try:
            self.load_state(folder_path + 'b2fstate', verbal)
        except FileNotFoundError:
            print("Error: b2fstate not found.")
        self.folder_path = folder_path

    def get_plot_params(self, var, nx, ny, charge_state=1, flatten=True):
        """
        Helper function to return a dictionary containing the data to be plotted
        formatted in a correct way, along with some associated labels for the
        plot.

        Inputs
        var (str): Which plot to get parameters/data for.
        nx (int): Number of knots in each ring.
        ny (int): Number of rings.
        charge_state (int): For variables that have values for each charge state.

        Outputs
        plot_dict (dict): Dictionary with things like the data, colorbar label, etc.
        """

        # Dictionary that will be returned with all the plot goodies. Put in
        # some defaults as well.
        plot_dict = defaultdict(None)
        if var.lower() in ['density', 'ne']:

            # Need to reshape to account for multiple charge states.
            data = self.na.reshape((nx, ny, self.ns), order='F').copy()
            plot_dict['data'] = data[:, :, charge_state].T.flatten()

            # Axes labels.
            if charge_state == 1:
                plot_dict['cbar_label'] = 'Density D+ (m' + r'$^{-3}$)'
            else:
                plot_dict['cbar_label'] = 'Density D{}+ (m'.format(charge_state) + r'$^{-3}$)'

            # Choose a good default scale for the colorbar.
            plot_dict['cbar_scale'] = 'log'

        elif var.lower() == 'te':
            data = self.te.reshape((nx, ny), order='F').copy() / e
            #plot_dict['data'] = data.T.flatten()
            plot_dict['data'] = data
            plot_dict['cbar_label'] = 'T' + r'$\mathrm{_e}$' + ' (eV)'
            plot_dict['cbar_scale'] = 'linear'

        elif var.lower() == 'ti':
            data = self.ti.reshape((nx, ny), order='F').copy() / e
            #plot_dict['data'] = data.T.flatten()
            plot_dict['data'] = data
            plot_dict['cbar_label'] = 'T' + r'$\mathrm{_i}$' + ' (eV)'
            plot_dict['cbar_scale'] = 'linear'

        elif var.lower() == 'bpol':
            data = self.bpol.reshape((nx, ny), order='F').copy()
            #plot_dict['data'] = data.T.flatten()
            plot_dict['data'] = data
            plot_dict['cbar_label'] = 'Bpol (T)'
            plot_dict['cbar_scale'] = 'linear'

        elif var.lower() == 'btor':
            data = self.btor.reshape((nx, ny), order='F').copy()
            #plot_dict['data'] = data.T.flatten()
            plot_dict['data'] = data
            plot_dict['cbar_label'] = 'Btor (T)'
            plot_dict['cbar_scale'] = 'linear'

        elif var.lower() == 'vpar':
            data = self.ua.reshape((nx, ny, self.ns), order='F').copy()
            #plot_dict['data'] = data[:, :, charge_state].T.flatten()
            plot_dict['data'] = data[:, :, charge_state]
            plot_dict['cbar_label'] = r'v$\mathrm{_{||}}$ (m/s)'
            plot_dict['cbar_scale'] = 'symlog'

        elif var.lower() == 'par_current':
            data = self.fchp.reshape((nx, ny, self.ns), order='F').copy()
            #plot_dict['data'] = data[:, :, charge_state].T.flatten()
            plot_dict['data'] = data[:, :, charge_state]
            plot_dict['cbar_label'] = r'$\mathrm{I_{par}}$ (A)'
            plot_dict['cbar_scale'] = 'linear'

        elif var.lower() == 'mach':
            print("Note: Assuming D+. Incorrect for anything else.")
            vpar = self.ua.reshape((nx, ny, self.ns), order='F').copy()
            te = self.te.reshape((nx, ny), order='F').copy() / e
            ti = self.ti.reshape((nx, ny), order='F').copy() / e
            cs = np.sqrt((te + ti) / md)
            mach = vpar[:, :, charge_state] / cs
            #plot_dict['data'] = mach.T.flatten()
            plot_dict['data'] = mach
            plot_dict['cbar_label'] = 'Mach Number'
            plot_dict['cbar_scale'] = 'symlog'

        elif var.lower() == 'psi':
            data = self.fpsi.reshape((nx, ny, 4), order='F').copy()
            plot_dict['data'] = data
            plot_dict['cbar_label'] = 'Psi'
            plot_dict['cbar_scale'] = 'linear'

        else:
            try:
                print("Warning: Variable not yet set up for plotting.")

                # Try and figure out the shape.
                data = getattr(self, var)
                if len(data) == nx * ny:
                    data = data.reshape((nx, ny), order='F').copy()
                elif len(data) == nx * ny * self.ns:
                    data = data.reshape((nx, ny, self.ns), order='F').copy()
                    data = data[:, :, charge_state]
                else:
                    print("Could not determine what to reshape {} into.".format(var))
                #plot_dict['data'] = data.T.flatten()
                plot_dict['data'] = data
                plot_dict['cbar_label'] = None
                plot_dict['cbar_scale'] = 'linear'
            except:
                print("Error: Could not find variable for plotting.")
                return None

        if flatten:
            plot_dict['data'] = plot_dict['data'].T.flatten()

        return plot_dict

    def plot_polygon(self, var='density', charge_state=1, fontsize=18,
                     normtype=None, cbar_lims=[None, None], lut=None,
                     cmap='PuRd', irsep=18):
        """
        The primary plotting fuction for plotting 2D data. The only required
        variable is "var", while all the others are for tweaking the plots. For
        each var a default set of plot parameters are chosen to act as a
        starting point for making pretty plots, but tweaking can be done using
        the other optional input parameters.

        Inputs
        var (str): Which plot type you want plotted. Not all variables from the
          output files are explicitly supported, but you can still enter those
          in. The reason is some variables require some extra formatting (like
          if it has repeat arrays, one for each charge state), so we handle all
          that behind the scenes here. But nonetheless, it is easy to add
          additional plot options since you can just use any of the other
          supported options as a template.
          Current supported options include, as strings (a * indicates
          charge_state is required): density*, te, bpol, btor, vpar, mach,
          par_current, mach, ...
        charge_state (int): Which charge state to plot, if required.
        fontsize (int): Font size for the plot labels.
        normtype (str): What type of colorbar scale to use. One of 'linear',
          'log', 'symlin', 'symlog'.
        cbar_lims (list, float): A 2-element list of the limits to use for the
          colorbar, i.e. [low, high].
        lut (int): How many levels to break the colobar into if you would like a
          discretized colorbar instead of a continuous one.
        cmap (str): What colormap to use. Check the matplotlib site for a
          listing of options.
        irsep (int): The ring index (i.e. y) of the separatrix. Typically 18,
          but no reason that can't be different.

        Outputs
        fig (matplotlib.pyplot.figure): The plot figure that is shown.
        """

        # First need to add 2 to the nx and ny to account for guard cells from
        # the grid making procedure in Carre.
        nx = self.nx + 2; ny = self.ny + 2

        # Reshape our X, Y (or knot, ring) data in a not at all obvious way.
        crx = self.crx.reshape((nx, ny, 4), order='F').copy()
        cry = self.cry.reshape((nx, ny, 4), order='F').copy()
        leftiy = self.leftiy.reshape((nx, ny), order='F').copy()

        # Create Patch objects of each cell out of the vertices for R, Z.
        patches = []
        for iy in np.arange(0, ny):
            for ix in np.arange(0, nx):
                rcol = crx[ix, iy, [0, 1, 3, 2]]
                zcol = cry[ix, iy, [0, 1, 3, 2]]
                rcol.shape=(4, 1)
                zcol.shape=(4, 1)
                polygon = Polygon(np.column_stack((rcol, zcol)), True, linewidth=0)
                patches.append(polygon)

        # Get our dictionary with some plotting values for this variable.
        plot_dict = self.get_plot_params(var, nx, ny, charge_state)
        self.plot_dict = plot_dict
        data = plot_dict['data']

        if normtype == None:
            normtype = plot_dict['cbar_scale']

        # First choose what kind of colorbar we will be using.
        if normtype == 'linear':
            if cbar_lims[0] == None:
                cbar_lims[0] = data.min()
            if cbar_lims[1] == None:
                cbar_lims[1] = data.max()
            norm = mpl.colors.Normalize(vmin=cbar_lims[0], vmax=cbar_lims[1])

        elif normtype == 'log':
            #data[data == 0.0] = 1e-3
            if cbar_lims[0] == None:
                cbar_lims[0] = data.min()
            if cbar_lims[1] == None:
                cbar_lims[1] = data.max()
            norm = mpl.colors.LogNorm(vmin=cbar_lims[0], vmax=cbar_lims[1])

        elif normtype == 'symlin':
            if cbar_lims[0] == None:
                cbar_lims[0] = -np.abs(data).max()
            if cbar_lims[1] == None:
                cbar_lims[1] = -cbar_lims[0]
            norm = mpl.colors.Normalize(vmin=cbar_lims[0], vmax=cbar_lims[1])
            cmap = 'coolwarm'

        elif normtype == 'symlog':
            data[data == 0.0] = 1e-3
            if cbar_lims[0] == None:
                cbar_lims[0] = -np.abs(data[~np.isnan(data)]).max()
            if cbar_lims[1] == None:
                cbar_lims[1] = -cbar_lims[0]
            norm = mpl.colors.SymLogNorm(linthresh=0.01 * cbar_lims[1], vmin=cbar_lims[0], vmax=cbar_lims[1])
            cmap = 'coolwarm'

        # Choose whether to discretize the colormap or not.
        if lut == None:
            scalar_map = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        else:
            scalar_map = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.get_cmap(cmap, lut=lut))

        # Create a PatchCollection object and add to our Axes.
        p = PatchCollection(patches, True, array=data, cmap=scalar_map.cmap,
                            norm=scalar_map.norm)

        fig, ax = plt.subplots(figsize=(6, 8))
        ax.add_collection(p)
        ax.plot()

        # Plot the separatrix.
        sep_x = crx[:, :, 0][leftiy==irsep]
        sep_y = cry[:, :, 0][leftiy==irsep]
        ax.plot(sep_x, sep_y, 'k')

        #cbar = plt.colorbar(p, ax=ax, pad=0.01)
        cbar = plt.colorbar(scalar_map, ax=ax, pad=0.01)
        cbar.ax.set_ylabel(plot_dict['cbar_label'], fontsize=fontsize)

        # Make axis the same relative scale (i.e. so a circle would be a circle).
        ax.axis('equal')
        ax.set_xlabel("R (m)", fontsize=fontsize)
        ax.set_ylabel("Z (m)", fontsize=fontsize)
        fig.tight_layout()
        fig.show()

        return fig

    def plot_along_ring(self, var, ring_num, charge_state=1, fontsize=18, lw=3,
                        logx=False, logy=False, ylims=[None, None],
                        xlims=[None, None]):
        """
        Plot a variable along a specific ring.

        Inputs:
        var (str): Variable name to plot. Check the plot_polygon documentation for a
          list of optional commands.
        ring_num (int): The ring number to plot along.
        charge_state (int): For vars that require it.
        fontsize (int): Font size for the plot.
        lw (int): Line width for the plot.
        logx, logy (bool): Log scales for axes.
        xlims, ylims (list, float): Limits for the x and y axes. Enter as
          [low, high].

        Outputs
        fig (matplotlib.pyplot.figure): The figure object so you can mess with it some more.
        """

        # Like above.
        nx = self.nx + 2; ny = self.ny + 2
        crx = self.crx.reshape((nx, ny, 4), order='F').copy()
        cry = self.cry.reshape((nx, ny, 4), order='F').copy()

        # Grab the data, reshaped.
        plot_dict = self.get_plot_params(var, nx, ny, charge_state)
        data = plot_dict['data'].reshape((nx, ny), order='F').copy()
        data = data[:, ring_num]
        crx_ring = crx[:, ring_num]; cry_ring = cry[:, ring_num]

        # The parallel s coordinate is not given, so let's just calculate it
        # here going from the centers of each edge of each cell. The indices of
        # each vertex is shown. S should increase from (1, 3) to (0, 2).
        #
        # 1 __________________________0
        #  |            s            |
        #  |-------------------------|
        #  |_________________________|
        # 3                          2
        #
        # Radial
        #   ^
        #   |
        #   --> Poloidal

        # To go from poloidal s to parallel s, we need to multiply by B/Bp. For
        # bb, (ix, iy, 0:3) are Bx, By and Bz. (ix, iy, 3) is total B.
        bb = self.bb.reshape((nx, ny, 4), order='F').copy()
        b = bb[:, ring_num, 3]
        bx = bb[:, ring_num, 0]; by = bb[:, ring_num, 1]

        # Bpol is sqrt(Bx^2 + By^2).
        bp = np.sqrt(np.square(bx) + np.square(by))

        spol = []; spar = []; ds = []
        for ik in range(0, len(data)):

            # Get coordinates of midpoint between the appropriate vertices.
            x1 = (crx_ring[ik][1] + crx_ring[ik][3]) / 2; y1 = (cry_ring[ik][1] + cry_ring[ik][3]) / 2
            x2 = (crx_ring[ik][0] + crx_ring[ik][2]) / 2; y2 = (cry_ring[ik][0] + cry_ring[ik][2]) / 2

            # Get distance between them.
            d = np.sqrt(np.square(x2 - x1) + np.square(y2 - y1))
            ds.append(d)

            # Append only half that distance (i.e. the cell center) to our s
            # list, adding onto the previous value.
            if ik == 0:
                spol.append(d / 2)
                spar.append((d / 2) * b[ik] / bp[ik])
            else:
                spol.append(spol[ik-1] + (d / 2) + (ds[ik-1] / 2))
                spar.append(spar[ik-1] + ((d / 2) + (ds[ik-1] / 2)) * b[ik] / bp[ik])

        # Plotting commands.
        fig, ax = plt.subplots()
        ax.plot(spar, data, color=tableau20[6], lw=lw)
        ax.set_xlabel("Parallel Distance (m)", fontsize=fontsize)
        ax.set_ylabel(plot_dict['cbar_label'], fontsize=fontsize)
        if logx:
            ax.set_xscale('log')
        if logy:
            if data.min() < 0:
                ax.set_yscale('symlog')
            else:
                ax.set_yscale('log')
        if xlims not in ([None, None], (None, None)):
            ax.set_xlim(xlims)
        if ylims != [None, None]:
            ax.set_ylim(ylims)
        fig.tight_layout()
        fig.show()

        return fig

    def identify_cells(self, rings=None, knots=None, irsep=18):
        """
        Simple helper function to help find where a particular, ring, knots or
        any combination of the two are.

        Inputs
        rings (list, int): Can either be a single ring, e.g. "23", or a list of
          ring indices, e.g. [23, 24, 25].
        knots (list, int): Likewise, but for knots. Can combine the two to
          identify a particular cell for example.

        Outputs
        None
        """

        # First need to add 2 to the nx and ny to account for guard cells from
        # the grid making procedure in Carre.
        nx = self.nx + 2; ny = self.ny + 2

        # Reshape our X, Y (or knot, ring) data in a not at all obvious way.
        crx = self.crx.reshape((nx, ny, 4), order='F').copy()
        cry = self.cry.reshape((nx, ny, 4), order='F').copy()
        leftiy = self.leftiy.reshape((nx, ny), order='F').copy()

        data = np.zeros((nx, ny))
        if type(knots) == type(None):
            data[:, rings] = 1
        elif type(rings) == type(None):
            data[knots, :] = 1
        else:
            data[knots, rings] = 1
        data = data.T.flatten()

        # Create Patch objects of each cell out of the vertices for R, Z.
        patches = []
        for iy in np.arange(0, ny):
            for ix in np.arange(0, nx):
                rcol = crx[ix, iy, [0, 1, 3, 2]]
                zcol = cry[ix, iy, [0, 1, 3, 2]]
                rcol.shape=(4, 1)
                zcol.shape=(4, 1)
                polygon = Polygon(np.column_stack((rcol, zcol)), True, linewidth=0.25, edgecolor='k')
                patches.append(polygon)

        scalar_map = mpl.cm.ScalarMappable(cmap=mpl.cm.get_cmap('Reds', lut=2))
        p = PatchCollection(patches, True, array=data, cmap=scalar_map.cmap,
                            norm=scalar_map.norm)

        fig, ax = plt.subplots(figsize=(6, 8))
        ax.add_collection(p)
        ax.plot()

        # Plot the separatrix.
        sep_x = crx[:, :, 0][leftiy==irsep]
        sep_y = cry[:, :, 0][leftiy==irsep]
        ax.plot(sep_x, sep_y, 'k')

        cbar = plt.colorbar(scalar_map, ax=ax, pad=0.01)
        #cbar.ax.set_ylabel(plot_dict['cbar_label'], fontsize=fontsize)

        # Make axis the same relative scale (i.e. so a circle would be a circle).
        ax.axis('equal')
        ax.set_xlabel("R (m)", fontsize=18)
        ax.set_ylabel("Z (m)", fontsize=18)
        fig.tight_layout()
        fig.show()

    def plot_radial(self, var, knots, charge_state=1, irsep=18, xaxis='rminrsep',
                    fontsize=18):
        """

        """
        # First need to add 2 to the nx and ny to account for guard cells from
        # the grid making procedure in Carre.
        nx = self.nx + 2; ny = self.ny + 2

        # Reshape our X, Y (or knot, ring) data in a not at all obvious way.
        crx = self.crx.reshape((nx, ny, 4), order='F').copy()
        cry = self.cry.reshape((nx, ny, 4), order='F').copy()
        leftiy = self.leftiy.reshape((nx, ny), order='F').copy()

        # Get our dictionary with some plotting values for this variable.
        plot_dict = self.get_plot_params(var, nx, ny, charge_state, flatten=False)
        self.plot_dict = plot_dict
        data = plot_dict['data']
        #data = data.reshape((nx, ny))
        data = data[knots, irsep:].T.flatten()

        # Need the center R value of each cell. This is just the mean of the x values.
        rs = crx[knots, irsep:].mean(axis=1)

        fig, ax = plt.subplots()
        xlabel = 'R (m)'
        if xaxis.lower() == 'rminrsep':
            sep_x0 = crx[:, :, 0][leftiy==irsep][knots]
            sep_x1 = crx[:, :, 1][leftiy==irsep][knots]
            rsep = (sep_x0 + sep_x1) / 2
            rs = (rs - rsep) * 100
            ax.axvline(0, color='k', linestyle='--')
            xlabel = 'R-Rsep (cm)'
        elif xaxis.lower() == 'psin':
            pass
        ax.plot(rs, data, color=tableau20[6], lw=3)
        ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.set_ylabel(plot_dict['cbar_label'], fontsize=fontsize)
        fig.tight_layout()
        fig.show()

        return fig
