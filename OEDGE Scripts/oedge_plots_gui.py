import tkinter as tk
from tkinter import filedialog
import oedge_plots as oedge


plot_opts = ['B Ratio', 'E Radial', 'E Poloidal', 'ExB Poloidal', 'ExB Radial',
             'Flow Velocity', 'Impurity Density', 'Impurity Flow Velocity',
             'Impurity Ionization', 'ne', 'ne - Divertor', 'Te', 'Te - Divertor']
plot_opts_cp = ['R-Rsep OMP vs. Flux - Midplane', 'R-Rsep OMP vs. Flux - Crown']

# Spacing constants for padding.
padx = 3
pady = 3

class Window(tk.Frame):

    def __init__(self, master=None):

        # This line does something so you can like use the Tk.Frame methods or
        # something. All the examples have it at least.
        super().__init__(master)
        self.master = master
        self.master.title('OEDGE Plotting GUI')
        self.netcdf_loaded = False
        self.cp_cb_mid_var = tk.IntVar()
        self.cp_cb_top_var = tk.IntVar()
        self.cp_cb_dim_var = tk.IntVar()
        self.create_widgets()


    def create_widgets(self):

        # Create an Entry box for location of netCDF file.
        tk.Label(self.master, text='NetCDF File: ').grid(row=0, column=0, sticky='WE', padx=padx, pady=pady)
        self.netcdf_entry = tk.Entry(self.master)
        self.netcdf_entry.grid(row=0, column=1, padx=padx, pady=pady, sticky='WE')

        # Add a Browse button next to it.
        self.netcdf_button = tk.Button(self.master, text='Browse...')
        self.netcdf_button.grid(row=0, column=2, padx=padx, pady=pady, sticky='WE')
        self.netcdf_button['command'] = self.browse_netcdf

        # Entry for dat file.
        tk.Label(self.master, text='Dat File: ').grid(row=1, column=0, sticky='E', padx=padx, pady=pady)
        self.dat_entry = tk.Entry(self.master)
        self.dat_entry.grid(row=1, column=1, padx=padx, pady=pady, sticky='WE')

        # Add a Browse button next to it as well.
        self.dat_button = tk.Button(self.master, text='Browse...')
        self.dat_button.grid(row=1, column=2, padx=padx, pady=pady, sticky='WE')
        self.dat_button['command'] = self.browse_dat

        # Add a drop down of what kind of plot to plot.
        tk.Label(self.master, text='Plot: ').grid(row=2, column=0, sticky='E', padx=padx, pady=pady)
        self.current_option = tk.StringVar(self.master)
        self.current_option.set(plot_opts[0])
        self.plot_option = tk.OptionMenu(self.master, self.current_option, *plot_opts)
        self.plot_option.grid(row=2, column=1, sticky='WE', padx=padx, pady=pady)

        # Put button to the right to do the plot.
        self.plot_button = tk.Button(self.master, text='Plot')
        self.plot_button.grid(row=2, column=2, padx=padx, pady=pady, sticky='WE')
        self.plot_button['command'] = self.plot_command

        # Add a quit button.
        self.quit_button = tk.Button(self.master, text='Quit')
        self.quit_button.grid(row=98, column=0, padx=padx, pady=pady*10, columnspan=3, sticky='S')
        self.quit_button['command'] = self.quit_command

        # Add a message box to act as readout for errors or such.
        self.message_box = tk.Text(self.master, height=15, width=47)
        self.message_box.grid(row=0, column=3, rowspan=99, padx=padx, pady=pady)
        self.message_box.insert(tk.END, "Click 'Browse...' to load path to netCDF file.\n")

        # Add scrollbar to message box.
        self.scroll = tk.Scrollbar(self.master)
        self.scroll.grid(row=0, column=4, rowspan=99, pady=pady, sticky='NS')
        self.scroll.config(command=self.message_box.yview)
        self.message_box.config(yscrollcommand=self.scroll.set)

        # Add button for plot options if the defaults aren't good enough.
        self.extra_plot_button = tk.Button(self.master, text='Plot Options...')
        self.extra_plot_button.grid(row=3, column=1, columnspan=2, padx=padx, pady=pady*4)
        self.extra_plot_button['command'] = self.extra_plot_opts

        # Entry for collectorprobe file.
        tk.Label(self.master, text='CP File: ').grid(row=4, column=0, sticky='E', padx=padx, pady=pady)
        self.cp_entry = tk.Entry(self.master)
        self.cp_entry.grid(row=4, column=1, padx=padx, pady=pady, sticky='WE')

        # Add a Browse button next to it as well.
        self.cp_button = tk.Button(self.master, text='Browse...')
        self.cp_button.grid(row=4, column=2, padx=padx, pady=pady, sticky='WE')
        self.cp_button['command'] = self.browse_cp

        # Add a drop down of what kind of plot to plot.
        tk.Label(self.master, text='Plot: ').grid(row=5, column=0, sticky='E', padx=padx, pady=pady)
        self.current_option_cp = tk.StringVar(self.master)
        self.current_option_cp.set(plot_opts_cp[0])
        self.plot_option_cp = tk.OptionMenu(self.master, self.current_option_cp, *plot_opts_cp)
        self.plot_option_cp.grid(row=5, column=1, sticky='WE', padx=padx, pady=pady)

        # Put button to the right to do the plot.
        self.plot_button_cp = tk.Button(self.master, text='Plot')
        self.plot_button_cp.grid(row=5, column=2, padx=padx, pady=pady, sticky='WE')
        self.plot_button_cp['command'] = self.plot_command_cp

    def browse_netcdf(self):
        """
        Function linked to 'Browse...' button to get location of netCDF file.
        """
        self.message_box.insert(tk.END, 'Loading...\n')
        root = tk.Tk(); root.withdraw()
        netcdf_path = tk.filedialog.askopenfilename(filetypes=(('NetCDF files', '*.nc'),))
        self.netcdf_entry.delete(0, tk.END)
        self.netcdf_entry.insert(0, netcdf_path)
        self.op = oedge.OedgePlots(self.netcdf_entry.get())
        self.message_box.insert(tk.END, 'Loaded file: {}\n'.format(netcdf_path))

        # Try and grab the collector_probe file while we're at it since it's probably
        # the same name. Dat file too.
        cp_path = netcdf_path.split('.nc')[0] + '.collector_probe'
        try:
            f = open(cp_path, 'r')
            self.cp_entry.delete(0, tk.END)
            self.cp_entry.insert(0, cp_path)
            self.message_box.insert(tk.END, 'Loaded file: {}\n'.format(cp_path))
            self.cp_path = cp_path
        except:
            pass
        dat_path = netcdf_path.split('.nc')[0] + '.dat'
        try:
            f = open(dat_path, 'r')
            self.dat_entry.delete(0, tk.END)
            self.dat_entry.insert(0, dat_path)
            self.op.add_dat_file(dat_path)
            self.message_box.insert(tk.END, 'Loaded file: {}\n'.format(dat_path))
            self.dat_path = dat_path
        except:
            pass


    def quit_command(self):
        """
        Quit commands to avoid hanging up the terminal or anything.
        """
        self.master.quit()
        self.master.destroy()

    def extra_plot_opts(self):
        """
        Additional plot options to turn on and off via checkbuttons.
        """

        # Some constants.
        # R of tip of MiMES probe.
        rtip = 2.26
        ztip_dim = -1.10
        ztip_top = 0.93

        # Open a new window.
        self.opt_window = tk.Toplevel(self.master)
        self.opt_window.title("Plot Options")

        # Create check boxes to turn on and off options.
        # Show MiMES probe.
        self.cp_cb_mid = tk.Checkbutton(self.opt_window, text='Show Collector Probe (MiMES) - Tip at R =', variable=self.cp_cb_mid_var, command=self.cp_cb_mid_com)
        self.cp_cb_mid.grid(row=0, column=0, padx=padx, pady=pady, sticky='W')
        self.cp_mid_entry = tk.Entry(self.opt_window)
        self.cp_mid_entry.insert(0, rtip)
        self.cp_mid_entry.grid(row=0, column=1, padx=padx, pady=pady)

        # Show DiMES probe.
        self.cp_cb_dim = tk.Checkbutton(self.opt_window, text='Show Collector Probe (DiMES) - Tip at Z =', variable=self.cp_cb_dim_var, command=self.cp_cb_dim_com)
        self.cp_cb_dim.grid(row=1, column=0, padx=padx, pady=pady, sticky='W')
        self.cp_dim_entry = tk.Entry(self.opt_window)
        self.cp_dim_entry.insert(0, ztip_dim)
        self.cp_dim_entry.grid(row=1, column=1, padx=padx, pady=pady)

        # Show top probe.
        self.cp_cb_top = tk.Checkbutton(self.opt_window, text='Show Collector Probe (top)      - Tip at Z =', variable=self.cp_cb_top_var, command=self.cp_cb_top_com)
        self.cp_cb_top.grid(row=2, column=0, padx=padx, pady=pady, sticky='W')
        self.cp_top_entry = tk.Entry(self.opt_window)
        self.cp_top_entry.insert(0, ztip_top)
        self.cp_top_entry.grid(row=2, column=1, padx=padx, pady=pady)

        # Show the metal rings.
        self.mr_cb_var = tk.IntVar()
        self.mr_cb = tk.Checkbutton(self.opt_window, text='Show Metal Rings', variable=self.mr_cb_var, command=self.mr_cb_com)
        self.mr_cb.grid(row=3, column=0, padx=padx, pady=pady, sticky='W')

    def cp_cb_mid_com(self):
        """
        Show the MiMES collector probe on the plot.
        """
        pass

    def cp_cb_dim_com(self):
        """
        Show the DiMES collector probe on the plot.
        """
        pass

    def cp_cb_top_com(self):
        """
        Show the top crown collector probe on the plot.
        """
        pass

    def mr_cb_com(self):
        pass

    def browse_cp(self):
        """
        Function linked to 'Browse...' button to get location of collectorprobe file.
        """
        self.message_box.insert(tk.END, 'Loading...\n')
        root = tk.Tk(); root.withdraw()
        self.cp_path = tk.filedialog.askopenfilename(filetypes=(('Collector Probe files', '*.collector_probe'),))
        self.cp_entry.delete(0, tk.END)
        self.cp_entry.insert(0, self.cp_path)
        #self.op = oedge.OedgePlots(self.netcdf_entry.get())
        self.message_box.insert(tk.END, 'Loaded file: {}\n'.format(self.cp_path))

    def browse_dat(self):
        """
        Function linked to 'Browse...' button to get location of dat file.
        """
        self.message_box.insert(tk.END, 'Loading...\n')
        root = tk.Tk(); root.withdraw()
        self.dat_path = tk.filedialog.askopenfilename(filetypes=(('Dat files', '*.dat'),))
        self.dat_entry.delete(0, tk.END)
        self.dat_entry.insert(0, self.dat_path)
        self.op.add_dat_file(self.dat_path)
        self.message_box.insert(tk.END, 'Loaded file: {}\n'.format(self.dat_path))


    def plot_command_cp(self):
        """
        Plotting commands for the collector probe plots.
        """

        if self.current_option_cp.get() == 'R-Rsep OMP vs. Flux - Midplane':
            plot_args = {'cp_path':self.cp_path, 'xaxis':'ROMP', 'yaxis':'IMPFLUX', 'cp_num':1, 'log':True}

        elif self.current_option_cp.get() == 'R-Rsep OMP vs. Flux - Crown':
            plot_args = {'cp_path':self.cp_path, 'xaxis':'ROMP', 'yaxis':'IMPFLUX', 'cp_num':2, 'log':True}

        else:
            self.message_box.insert(tk.END, 'Plot option not found.')

        self.op.cp_plots(**plot_args)

    def plot_command(self):
        """
        For the selected DIVIMP plot, pass the correct parameters to the plotting
        function.
        """

        if self.current_option.get() == 'Te':
            plot_args = {'dataname'  :'KTEBS',
                         'cmap'      :'inferno',
                         'cbar_label':'Te (eV)',
                         'normtype'  :'log'}

        elif self.current_option.get() == 'ne':
            plot_args = {'dataname'  :'KNBS',
                         'cmap'      :'viridis',
                         'cbar_label':'ne (m-3)',
                         'normtype'  :'log'}

        if self.current_option.get() == 'Te - Divertor':
            plot_args = {'dataname'  :'KTEBS',
                         'cmap'      :'inferno',
                         'cbar_label':'Te (eV)',
                         'normtype'  :'log',
                         'xlim'      :[1.3, 1.75],
                         'ylim'      :[-1.4, -0.9]}

        elif self.current_option.get() == 'ne - Divertor':
            plot_args = {'dataname'  :'KNBS',
                         'cmap'      :'viridis',
                         'cbar_label':'ne (m-3)',
                         'normtype'  :'log',
                         'xlim'      :[1.3, 1.75],
                         'ylim'      :[-1.4, -0.9]}

        elif self.current_option.get() == 'ExB Radial':
            plot_args = {'dataname'  :'EXB_R',
                         'cbar_label':'ExB Radial Drift (m/s?)',
                         'normtype'  :'symlog',
                         'scaling':1.0 / self.op.qtim}

        elif self.current_option.get() == 'ExB Poloidal':
            plot_args = {'dataname'  :'EXB_P',
                         'cbar_label':'ExB Poloidal Drift (m/s?)',
                         'normtype'  :'symlog',
                         'scaling':1.0 / self.op.qtim}

        elif self.current_option.get() == 'B Ratio':
            plot_args = {'dataname'  :'BRATIO',
                         'cbar_label':'B Ratio'}

        elif self.current_option.get() == 'Flow Velocity':

            # Flow velocity should be scaled by 1/QTIM.
            scaling = 1.0 / self.op.qtim
            plot_args = {'dataname'  :'KVHS',
                         'cbar_label':'Flow Velocity (m/s)',
                         'normtype'  :'symlog',
                         'scaling'   :scaling}

        elif self.current_option.get() == 'Impurity Flow Velocity':
            plot_args = {'dataname'  :'KVHSimp',
                         'cbar_label':'Impurity Flow Velocity (m/s)',
                         'normtype'  :'symlog'}

        elif self.current_option.get() == 'E Radial':
            plot_args = {'dataname'  :'E_RAD',
                         'cbar_label':'E Radial (V/m)',
                         'normtype'  :'symlin'}

        elif self.current_option.get() == 'E Poloidal':
            plot_args = {'dataname'  :'E_POL',
                         'cbar_label':'E Poloidal (V/m)',
                         'normtype'  :'symlin'}

        elif self.current_option.get() == 'Impurity Density':

            # Impurity density is scaled by ABSFAC.
            scaling = self.op.absfac
            plot_args = {'dataname'  :'DDLIMS',
                         'cbar_label':'Tungsten Density (m-3)',
                         'normtype'  :'log',
                         'charge'    :'all',
                         'vmin'      :1e13,
                         'vmax'      :1e17,
                         'scaling'   :scaling,
                         'cmap'      :'nipy_spectral'}

        elif self.current_option.get() == 'Impurity Ionization':
            plot_args = {'dataname'  :'KFIZS',
                         'cbar_label':'Impurity Ionization Rate',
                         'charge'    :30,
                         'normtype'  :'log'}

        else:
            self.message_box.insert(tk.END, 'Plot option not found.')

        # Options for including collector probes or metal rings. Add into
        # dictionary if we want them.
        if self.cp_cb_mid_var.get() == 1:

            # Append the probe and tip value on if a list has already been made.
            if 'show_cp' in plot_args.keys():
                plot_args['show_cp'].append(1)
                plot_args['ptip'].append(float(self.cp_mid_entry.get()))
            else:
                plot_args['show_cp'] = [1]
                plot_args['ptip'] = [float(self.cp_mid_entry.get())]

        if self.cp_cb_top_var.get() == 1:

            # Append the probe and tip value on if a list has already been made.
            if 'show_cp' in plot_args.keys():
                plot_args['show_cp'].append(2)
                plot_args['ptip'].append(float(self.cp_top_entry.get()))
            else:
                plot_args['show_cp'] = [2]
                plot_args['ptip'] = [float(self.cp_top_entry.get())]

        if self.cp_cb_dim_var.get() == 1:

            # Append the probe and tip value on if a list has already been made.
            if 'show_cp' in plot_args.keys():
                plot_args['show_cp'].append(3)
                plot_args['ptip'].append(float(self.cp_dim_entry.get()))
            else:
                plot_args['show_cp'] = [3]
                plot_args['ptip'] = [float(self.cp_dim_entry.get())]

        # Option to show metal rings.
        if self.mr_cb_var.get() == 1:
            plot_args['show_mr'] = True

        message = "Plotting " + self.current_option.get() + ".\n"
        self.message_box.insert(tk.END, message)
        fig = self.op.plot_contour_polygon(**plot_args)



def main():

    root = tk.Tk()
    window = Window(root)
    window.mainloop()

if __name__ == '__main__':
    main()
