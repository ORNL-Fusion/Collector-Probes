import tkinter as tk
from tkinter import filedialog
import oedge_plots as oedge
import numpy as np


plot_opts = ['B Ratio', 'E Radial', 'E Poloidal', 'ExB Poloidal', 'ExB Radial',
             'Flow Velocity', 'Flow Velocity (with T13)', 'Impurity Density',
             'Impurity Ionization', 'ne', 'ne - Divertor', 'Rings', 'S Coordinate', 'Te',
             'Te - Divertor']
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
        self.mr_cb_var     = tk.IntVar()
        self.charge_entry  = tk.Entry(self.master)
        self.charge_entry.insert(0, 30)
        self.create_widgets()


    def create_widgets(self):

        row = 0

        # Add a message box to act as readout for errors or such.
        self.message_box = tk.Text(self.master, height=25, width=55)
        self.message_box.grid(row=0, column=3, rowspan=99, padx=padx, pady=pady)
        self.message_box.insert(tk.END, "Click 'Browse...' to load path to netCDF file.\n")

        # Add scrollbar to message box.
        self.scroll = tk.Scrollbar(self.master)
        self.scroll.grid(row=0, column=4, rowspan=99, pady=pady, sticky='NS')
        self.scroll.config(command=self.message_box.yview)
        self.message_box.config(yscrollcommand=self.scroll.set)

        # Create an Entry box for location of netCDF file.
        tk.Label(self.master, text='NetCDF File: ').grid(row=row, column=0, sticky='WE', padx=padx, pady=pady)
        self.netcdf_entry = tk.Entry(self.master)
        self.netcdf_entry.grid(row=row, column=1, padx=padx, pady=pady, sticky='WE')

        # Add a Browse button next to it.
        self.netcdf_button = tk.Button(self.master, text='Browse...')
        self.netcdf_button.grid(row=row, column=2, padx=padx, pady=pady, sticky='WE')
        self.netcdf_button['command'] = self.browse_netcdf
        row += 1

        # Entry for dat file.
        tk.Label(self.master, text='Dat File: ').grid(row=row, column=0, sticky='E', padx=padx, pady=pady)
        self.dat_entry = tk.Entry(self.master)
        self.dat_entry.grid(row=row, column=1, padx=padx, pady=pady, sticky='WE')

        # Add a Browse button next to it as well.
        self.dat_button = tk.Button(self.master, text='Browse...')
        self.dat_button.grid(row=row, column=2, padx=padx, pady=pady, sticky='WE')
        self.dat_button['command'] = self.browse_dat
        row += 1

        # Add a drop down of what kind of plot to plot.
        tk.Label(self.master, text='Plot: ').grid(row=row, column=0, sticky='E', padx=padx, pady=pady)
        self.current_option = tk.StringVar(self.master)
        self.current_option.set(plot_opts[0])
        self.plot_option = tk.OptionMenu(self.master, self.current_option, *plot_opts)
        self.plot_option.grid(row=row, column=1, sticky='WE', padx=padx, pady=pady)

        # Put button to the right to do the plot.
        self.plot_button = tk.Button(self.master, text='Plot')
        self.plot_button.grid(row=row, column=2, padx=padx, pady=pady, sticky='WE')
        self.plot_button['command'] = self.plot_command
        row += 1

        # Add button for plot options if the defaults aren't good enough.
        self.extra_plot_button = tk.Button(self.master, text='Plot Options...')
        self.extra_plot_button.grid(row=row, column=1, columnspan=2, padx=padx, pady=pady*4)
        self.extra_plot_button['command'] = self.extra_plot_opts
        row += 1

        # Some spacing between sections.
        tk.Label(self.master, text=' ').grid(row=row, column=0, padx=padx, pady=pady)
        row +=1

        # Entry for collectorprobe file.
        tk.Label(self.master, text='CP File: ').grid(row=row, column=0, sticky='E', padx=padx, pady=pady)
        self.cp_entry = tk.Entry(self.master)
        self.cp_entry.grid(row=row, column=1, padx=padx, pady=pady, sticky='WE')

        # Add a Browse button next to it as well.
        self.cp_button = tk.Button(self.master, text='Browse...')
        self.cp_button.grid(row=row, column=2, padx=padx, pady=pady, sticky='WE')
        self.cp_button['command'] = self.browse_cp
        row += 1

        # Add a drop down of what kind of plot to plot.
        tk.Label(self.master, text='Plot: ').grid(row=row, column=0, sticky='E', padx=padx, pady=pady)
        self.current_option_cp = tk.StringVar(self.master)
        self.current_option_cp.set(plot_opts_cp[0])
        self.plot_option_cp = tk.OptionMenu(self.master, self.current_option_cp, *plot_opts_cp)
        self.plot_option_cp.grid(row=row, column=1, sticky='WE', padx=padx, pady=pady)

        # Put button to the right to do the plot.
        self.plot_button_cp = tk.Button(self.master, text='Plot')
        self.plot_button_cp.grid(row=row, column=2, padx=padx, pady=pady, sticky='WE')
        self.plot_button_cp['command'] = self.plot_command_cp
        row += 1

        # Some spacing between sections.
        tk.Label(self.master, text=' ').grid(row=row, column=0, padx=padx, pady=pady)
        row += 1

        # Entry for Thomson input file.
        tk.Label(self.master, text='Thomson\n Input File: ').grid(row=row, column=0, sticky='E', padx=padx, pady=pady, rowspan=2)
        self.ts_entry = tk.Entry(self.master)
        self.ts_entry.grid(row=row, column=1, padx=padx, pady=pady, sticky='WE', rowspan=2)

        # Add a Browse button next to it as well.
        self.ts_button = tk.Button(self.master, text='Browse...')
        self.ts_button.grid(row=row, column=2, padx=padx, pady=pady, sticky='WE')
        self.ts_button['command'] = self.browse_ts
        row += 1

        # Add a Browse button next to it as well.
        self.ts_button = tk.Button(self.master, text='Create...')
        self.ts_button.grid(row=row, column=2, padx=padx, pady=pady, sticky='WE')
        self.ts_button['command'] = self.create_ts
        row += 1

        # Entry for Thomson file.
        tk.Label(self.master, text='Thomson\n Output File: ').grid(row=row, column=0, sticky='E', padx=padx, pady=pady, rowspan=2)
        self.ts_out_entry = tk.Entry(self.master)
        self.ts_out_entry.grid(row=row, column=1, padx=padx, pady=pady, sticky='WE', rowspan=2)

        # Add button to make PDF of Thomson comparisons.
        self.compare_button = tk.Button(self.master, text='Plot')
        self.compare_button.grid(row=row, column=2, padx=padx, pady=pady*4, sticky='WE')
        self.compare_button['command'] = self.compare_ts_command
        row += 1

        # Add a quit button.
        self.quit_button = tk.Button(self.master, text='Quit')
        self.quit_button.grid(row=98, column=0, padx=padx, pady=pady*10, columnspan=3, sticky='S')
        self.quit_button['command'] = self.quit_command

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
        self.message_box.insert(tk.END, 'Loaded file: {}\n'.format(netcdf_path.split('/')[-1]))

        # Try and grab the collector_probe file while we're at it since it's probably
        # the same name. Dat file too.
        cp_path = netcdf_path.split('.nc')[0] + '.collector_probe'
        try:
            f = open(cp_path, 'r')
            self.cp_entry.delete(0, tk.END)
            self.cp_entry.insert(0, cp_path)
            self.message_box.insert(tk.END, 'Loaded file: {}\n'.format(cp_path.split('/')[-1]))
            self.cp_path = cp_path
        except:
            pass
        dat_path = netcdf_path.split('.nc')[0] + '.dat'
        try:
            f = open(dat_path, 'r')
            self.dat_entry.delete(0, tk.END)
            self.dat_entry.insert(0, dat_path)
            self.op.add_dat_file(dat_path)
            self.message_box.insert(tk.END, 'Loaded file: {}\n'.format(dat_path.split('/')[-1]))
            self.dat_path = dat_path
        except:
            pass

        # Add a generic name for the Thomson output file.
        ts_out = netcdf_path.split('.nc')[0] + '_ts.pdf'
        self.ts_out_entry.insert(0, ts_out)

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
        self.cp_cb_mid = tk.Checkbutton(self.opt_window, text='Show Collector Probe (MiMES) - Tip at R =', variable=self.cp_cb_mid_var)
        self.cp_cb_mid.grid(row=0, column=0, padx=padx, pady=pady, sticky='W')
        self.cp_mid_entry = tk.Entry(self.opt_window)
        self.cp_mid_entry.insert(0, rtip)
        self.cp_mid_entry.grid(row=0, column=1, padx=padx, pady=pady)

        # Show DiMES probe.
        self.cp_cb_dim = tk.Checkbutton(self.opt_window, text='Show Collector Probe (DiMES) - Tip at Z =', variable=self.cp_cb_dim_var)
        self.cp_cb_dim.grid(row=1, column=0, padx=padx, pady=pady, sticky='W')
        self.cp_dim_entry = tk.Entry(self.opt_window)
        self.cp_dim_entry.insert(0, ztip_dim)
        self.cp_dim_entry.grid(row=1, column=1, padx=padx, pady=pady)

        # Show top probe.
        self.cp_cb_top = tk.Checkbutton(self.opt_window, text='Show Collector Probe (top)      - Tip at Z =', variable=self.cp_cb_top_var)
        self.cp_cb_top.grid(row=2, column=0, padx=padx, pady=pady, sticky='W')
        self.cp_top_entry = tk.Entry(self.opt_window)
        self.cp_top_entry.insert(0, ztip_top)
        self.cp_top_entry.grid(row=2, column=1, padx=padx, pady=pady)

        # Show the metal rings.
        self.mr_cb_var = tk.IntVar()
        self.mr_cb = tk.Checkbutton(self.opt_window, text='Show Metal Rings', variable=self.mr_cb_var)
        self.mr_cb.grid(row=3, column=0, padx=padx, pady=pady, sticky='W')

        # Option to choose which charge state to plot.
        tk.Label(self.opt_window, text=' ').grid(row=4, column=0)
        tk.Label(self.opt_window, text='Charge State to Plot:').grid(row=5, column=0, padx=padx, pady=pady, sticky='W')
        self.charge_entry = tk.Entry(self.opt_window)
        self.charge_entry.insert(0, 30)
        self.charge_entry.grid(row=5, column=1, padx=padx, pady=pady)

        # Option to change vmin/vmax for plot colorbars.
        tk.Label(self.opt_window, text='Colorbar Scale Min:').grid(row=6, column=0, padx=padx, pady=pady, sticky='W')
        tk.Label(self.opt_window, text='Colorbar Scale Max:').grid(row=7, column=0, padx=padx, pady=pady, sticky='W')
        self.vmin_entry = tk.Entry(self.opt_window)
        self.vmax_entry = tk.Entry(self.opt_window)
        self.vmin_entry.insert(0, 'auto')
        self.vmax_entry.insert(0, 'auto')
        self.vmin_entry.grid(row=6, column=1, padx=padx, pady=pady)
        self.vmax_entry.grid(row=7, column=1, padx=padx, pady=pady)

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
        self.message_box.insert(tk.END, 'Loaded file: {}\n'.format(self.dat_path.split('/')[-1]))

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
        self.message_box.insert(tk.END, 'Loaded file: {}\n'.format(self.cp_path.split('/')[-1]))

    def browse_ts(self):
        """
        Function linked to 'Browse...' button to get location of Thomson file.
        """
        self.message_box.insert(tk.END, 'Loading...\n')
        root = tk.Tk(); root.withdraw()
        self.ts_path = tk.filedialog.askopenfilename(filetypes=(('Excel files', '*.xlsx'),))
        self.ts_entry.delete(0, tk.END)
        self.ts_entry.insert(0, self.ts_path)
        #self.op = oedge.OedgePlots(self.netcdf_entry.get())
        self.message_box.insert(tk.END, 'Loaded file: {}\n'.format(self.ts_path.split('/')[-1]))

    def create_ts(self):

        # Open a new window.
        self.create_window = tk.Toplevel(self.master)
        self.create_window.title("Create TS file")

        # Create Entry for shots to get TS data for.
        tk.Label(self.create_window, text='Shots (separated by commas):').grid(row=0, column=0, sticky='E', padx=padx, pady=pady)
        self.shot_entry = tk.Entry(self.create_window)
        self.shot_entry.grid(row=0, column=1, padx=padx, pady=pady, sticky='WE', columnspan=3)

        # Entry for times.
        tk.Label(self.create_window, text='Times (start, stop, step):').grid(row=1, column=0, sticky='E', padx=padx, pady=pady)
        self.time_start = tk.Entry(self.create_window, width=10)
        self.time_start.grid(row=1, column=1, padx=padx, pady=pady, sticky='W')
        self.time_end = tk.Entry(self.create_window, width=10)
        self.time_end.grid(row=1, column=2, padx=padx, pady=pady, sticky='W')
        self.time_step = tk.Entry(self.create_window, width=10)
        self.time_step.grid(row=1, column=3, padx=padx, pady=pady, sticky='W')

        # Entry for reference time to map all the TS measurements to.
        tk.Label(self.create_window, text='Reference time:').grid(row=2, column=0, sticky='E', padx=padx, pady=pady)
        self.time_ref = tk.Entry(self.create_window)
        self.time_ref.grid(row=2, column=1, padx=padx, pady=pady, sticky='WE', columnspan=2)

        # Entry for reference time to map all the TS measurements to.
        tk.Label(self.create_window, text='Filename:').grid(row=3, column=0, sticky='E', padx=padx, pady=pady)
        self.ts_filename = tk.Entry(self.create_window)
        self.ts_filename.grid(row=3, column=1, padx=padx, pady=pady, sticky='WE', columnspan=2)

        # Fill in some default values.
        self.shot_entry.insert(0, '167192, 167193, 167194, 167195')
        self.time_start.insert(0, 2500)
        self.time_end.insert(0, 5000)
        self.time_step.insert(0, 500)
        self.time_ref.insert(0, 3000)
        self.ts_filename.insert(0, 'my_ts_filename.xlsx')

        # Create a TS file with these parameters.
        self.ts_create_button = tk.Button(self.create_window, text='Create')
        self.ts_create_button.grid(row=4, column=1, padx=padx, pady=pady, sticky='WE')
        self.ts_create_button['command'] = self.create_ts_command

        # Note that you can use zero for step to get all TS times.
        tk.Label(self.create_window, text='Note: Use 0 in step to load all available times between start and stop.').grid(row=5, column=0, columnspan=4, padx=padx, pady=pady, sticky='W')

    def create_ts_command(self):

        # Printout to message box.
        self.message_box.insert(tk.END, 'Creating Thomson scattering data file... ')

        # Convert the shots input to a list.
        shots = self.shot_entry.get()
        shots = np.array(shots.split(','), dtype=np.int)

        if float(self.time_step.get()) == 0.0:
            load_all_ts = True
            times=np.arange(float(self.time_start.get()), float(self.time_end.get()), 10)
        else:
            load_all_ts = False
            times=np.arange(float(self.time_start.get()), float(self.time_end.get()), float(self.time_step.get()))


        fname = self.op.create_ts(shots=shots,
                                  times=times,
                                  ref_time=float(self.time_ref.get()),
                                  filename=self.ts_filename.get(),
                                  load_all_ts=load_all_ts)

        self.ts_entry.insert(0, fname)
        self.message_box.insert(tk.END, 'Done.\n')

    def compare_ts_command(self):
        """
        Generate the pdf of Thomson comparisons.
        """

        self.message_box.insert(tk.END, 'Generating PDF...')
        ts_filename = self.ts_entry.get()
        rings = np.append(np.arange(19, 49), np.arange(178, 190))
        self.op.compare_ts(ts_filename, rings, show_legend='short', output_file=self.ts_out_entry.get())
        self.message_box.insert(tk.END, ' Done.\n')

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

        elif self.current_option.get() == 'Flow Velocity (with T13)':
            plot_args = {'dataname'  :'KVHSimp',
                         'cbar_label':'Flow Velocity (m/s)',
                         'normtype'  :'symlog'}

        elif self.current_option.get() == 'E Radial':
            plot_args = {'dataname'  :'E_RAD',
                         'cbar_label':'E Radial (V/m)',
                         'normtype'  :'symlog'}

        elif self.current_option.get() == 'E Poloidal':
            plot_args = {'dataname'  :'E_POL',
                         'cbar_label':'E Poloidal (V/m)',
                         'normtype'  :'symlog'}

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
            charge = int(self.charge_entry.get())
            plot_args = {'dataname'  :'KFIZS',
                         'cbar_label':'Impurity Ionization Rate W{}+'.format(charge),
                         'charge'    :charge,
                         'normtype'  :'log',
                         'fix_fill'  :True}

        elif self.current_option.get() == 'S Coordinate':
            plot_args = {'dataname'  :'Snorm',
                         'cbar_label':'S Parallel (Normalized)',
                         'normtype'  :'linear',
                         'lut'       :10,
                         'cmap'      :'nipy_spectral'}

        elif self.current_option.get() == 'Rings':
            plot_args = {'dataname'  :'Ring',
                         'cbar_label':'Ring'}

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

        # Option to supply a vmin or vmax. Put in a try since this may not even
        # be defined yet if the extra plot option window isn't open.
        try:
            if self.vmin_entry.get() != 'auto':
                plot_args['vmin'] = float(self.vmin_entry.get())
            if self.vmax_entry.get() != 'auto':
                plot_args['vmax'] = float(self.vmax_entry.get())
        except:
            pass

        message = "Plotting " + self.current_option.get() + ".\n"
        self.message_box.insert(tk.END, message)
        fig = self.op.plot_contour_polygon(**plot_args)



def main():

    root = tk.Tk()
    window = Window(root)
    window.mainloop()

if __name__ == '__main__':
    main()
