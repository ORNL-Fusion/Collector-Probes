import tkinter as tk
from tkinter import filedialog
from tkinter import font
import numpy as np
import dlim_plots as dlim
import lim_plots as limpt

# Color name to indicate background color of section.
cname = 'gray75'
cname2 = 'gray90'

#spacing constants
padx = 3
pady = 4

#Width of columns
col1_width = 10
col2_width = 30
col3_width = 7

#plot options
plot_op = ['Center Line', 'Contour', 'Poloidal Profiles', 'Temperature']

class Window(tk.Frame):

    def __init__(self, master=None):

        # This line does something so you can like use the Tk.Frame methods or
        # something. All the examples have it at least.
        super().__init__(master)
        self.master = master
        self.master.title('3dlim plotting GUI')
        self.netcdf_loaded = False
        self.cp_cb_mid_var = tk.IntVar()
        self.cp_cb_top_var = tk.IntVar()
        self.cp_cb_dim_var = tk.IntVar()
        self.mr_cb_var     = tk.IntVar()
        self.core_cb_var   = tk.IntVar()
        self.charge_entry  = tk.Entry(self.master)
        self.charge_entry.insert(0, 30)
        self.vzmult_entry = tk.Entry(self.master)
        self.vzmult_entry.insert(0, 0)
        self.create_widgets()

    def create_widgets(self):
        """
        Create all the buttons, message box, etc. and lay them all out.
        """

        # Variable to keep track of row number so we don't have to!
        row = 0

        # Color name to indicate background color of section.


        #Width of columns
        col1_width = 10
        col2_width = 30
        col3_width = 7

        #Add a message box
        self.message_box = tk.Text(self.master, height = 7, width=65)
        self.message_box.grid(row=0, column=4, rowspan=5, padx=padx, sticky='NS', )
        self.message_box.insert(tk.END, "Click 'Browse...' to load path to netCDF file.\n")

        #Add scrollbar to message box
        self.scroll = tk.Scrollbar(self.master)
        self.scroll.grid(row=0, column=5, rowspan=5, pady=pady, sticky='NS')
        self.scroll.config(command=self.message_box.yview)
        self.message_box.config(yscrollcommand=self.scroll.set)

        #Make new frame for file stuff
        self.netcdf_frame = tk.Frame(self.master, bg=cname)
        self.netcdf_frame.grid(row=row, column=0, columnspan=4, sticky='WE')
        tk.Label(self.netcdf_frame, text='Input Files:', bg=cname, width=col1_width).grid(row=row, column=1, sticky='E')

        #Browse Button to browse for netcdf files
        self.netcdf_button = tk.Button(self.netcdf_frame, text='Browse...',width=col3_width)
        self.netcdf_button.grid(row=row, column=2, padx=padx, pady=pady, sticky='WE')
        self.netcdf_button['command'] = self.browse_netcdf

        row += 1

        #Place for netcdf file
        self.netcdf_entry = tk.Entry(self.netcdf_frame, width=col2_width)
        self.netcdf_entry.grid(row=row, column=0, columnspan=4, padx=padx, pady=pady, sticky='WE')

        tk.Label(self.netcdf_frame, text='.nc', bg=cname, width=col1_width, anchor='w').grid(row=row, column=4)

        row += 1

        #Place for dat file
        self.dat_entry = tk.Entry(self.netcdf_frame, width=col2_width)
        self.dat_entry.grid(row=row, column=0, columnspan=4, padx=padx, pady=pady, sticky='WE')

        tk.Label(self.netcdf_frame, text='.dat', bg=cname, width=col1_width, anchor='w').grid(row=row, column=4)

        row += 1

        #Place for lim file
        self.lim_entry = tk.Entry(self.netcdf_frame, width=col2_width)
        self.lim_entry.grid(row=row, column=0, columnspan=4, padx=padx, pady=pady, sticky='WE')

        tk.Label(self.netcdf_frame, text='.lim', bg=cname, width=col1_width, anchor='w').grid(row=row, column=4)

        row += 1

        #Plot overview button to plot simple graphs
        self.overview_button = tk.Button(self.netcdf_frame, text='Plot Overview')
        self.overview_button.grid(row=row, column=2, padx=padx, pady=pady, sticky='WE')

        row += 1

        #Make a new frame for plot selection
        self.sel_frame = tk.Frame(self.master, bg=cname2)
        self.sel_frame.grid(row=row, column=0, columnspan=4, sticky='WE')
        tk.Label(self.sel_frame, text='Plot Selection:', bg=cname2).grid(row=row, column=1)

        #Make a new frame for plot options
        self.opt_frame = tk.Frame(self.master, bg=cname)
        self.opt_frame.grid(row=row, column=4, sticky='WE')
        tk.Label(self.opt_frame, text='Additional Plot Options', bg=cname).grid(row=row, column=5)

        row += 1

        #Option Menu for different plots
        self.current_option = tk.StringVar(self.sel_frame)
        self.current_option.set(plot_op[0])
        self.plot_options = tk.OptionMenu(self.sel_frame, self.current_option, *plot_op)
        self.plot_options.grid(row=row, column=1, padx=padx, pady=pady)
        self.current_option.trace('w', self.option_selection)

        '''
        #Option Selection Button
        self.option_button = tk.Button(self.sel_frame, text='Extra Options')
        self.option_button.grid(row=row, column=2, padx=padx, pady=pady)
        self.option_button['command'] = self.option_selection
        '''

        row += 1

        #Plot button to plot the selected plot
        self.plot_button = tk.Button(self.sel_frame, text='Plot', width=col3_width)
        self.plot_button.grid(row=row, column=2, padx=padx, pady=pady)
        self.plot_button['command'] = self.plot_action

        row += 1

        tk.Label(self.master, text='').grid(row=row, column=1)

        row += 1

        self.quit_button = tk.Button(self.master, text='Quit')
        self.quit_button.grid(row=row, column=1)
        self.quit_button['command'] = self.quit_command


    def plot_action(self):
        """
        Function to take current options and plot the selected graph.
        """

        #checks to see if there is a file selected
        if self.netcdf_entry.get() == '':
            self.message_box.insert(tk.END, 'No netCDF File Selected\n')

        elif self.current_option.get() == 'Center Line':

            try:
                if self.centline_log == 1:
                    plot_args = {'log': True}
                else:
                    plot_args = {'log': False}

                if self.centline_exp ==1:
                    plot_args['fit_exp'] = True
                else:
                    plot_args['fit_exp'] = False

                self.dl.centerline(**plot_args)
            except:
                self.dl.centerline()

        elif self.current_option.get() == 'Contour':

            if self.side_option.get() == 2:
                plot_args = {'side': 'ITF'}
            else:
                plot_args = {'side': 'OTF'}

            if self.width_option.get() == 2:
                plot_args['probe_width'] = 0.005
            elif self.width_option.get() == 3:
                plot_args['probe_width'] = 0.0025
            elif self.width_option.get() == 4:
                try:
                    plot_args['probe_width'] = float(self.cont_widEnt.get())
                except:
                    self.message_box.insert(tk.END, 'Custom width failed.\n')
                    plot_args['probe_width'] = 0.015
            else:
                plot_args['probe_width'] = 0.015

            self.dl.deposition_contour(**plot_args)

        elif self.current_option.get() == 'Poloidal Profiles':

            if self.width_option.get() == 2:
                plot_args = {'probe_width': 0.005}
            elif self.width_option.get() == 3:
                plot_args = {'probe_width': 0.0025}
            elif self.width_option.get() == 4:
                try:
                    plot_args = {'probe_width': float(self.cont_widEnt.get())}
                except:
                    self.message_box.insert(tk.END, 'Custom width failed.\n')
                    plot_args = {'probe_width': 0.015}
            else:
                plot_args = {'probe_width': 0.015}

            self.dl.avg_pol_profiles(**plot_args)

        else:
            #Graphs temperature graph
            if self.current_option.get() == 'Temperature':

                #gets the ylim
                try:
                    ylim = int(self.temp_ylim.get())
                except:
                    ylim=500

                self.dl.plot_2d('CTEMBS', ylim=ylim)
            self.message_box.insert(tk.END, 'Plotted')



    def option_selection(self, var, ind, mode):
        """
        Function to make selection of plot and show additional plot options.
        """

        #The Center Line is selected this is what happens
        if self.current_option.get() == 'Center Line':

            #creates new frame for the center line selection
            self.delete_frames()
            self.opt_centline = tk.Frame(self.master, bg=cname)
            self.opt_centline.grid(row=5, column=4, columnspan=4, sticky='WE')

            #Adds log options checkbox
            tk.Label(self.opt_centline, text='Log plot: ', bg=cname, width=10, anchor='e').grid(row=5, column=4)
            self.log_option = tk.IntVar()
            self.centline_log = tk.Checkbutton(self.opt_centline, variable=self.log_option, onvalue=1, offvalue=0, bg=cname)
            self.centline_log.grid(row=5, column=5)

            #Adds exp options checkbox
            tk.Label(self.opt_centline, text='Exponential fit: ', bg=cname, anchor='e').grid(row=5, column=6)
            self.exp_option = tk.IntVar()
            self.centline_exp = tk.Checkbutton(self.opt_centline, variable=self.exp_option, onvalue=1, offvalue=0, bg=cname)
            self.centline_exp.grid(row=5, column=7)

        #The contour is selected this is what happens
        elif self.current_option.get() == 'Contour':

            #creates new frame for the contour selection
            self.delete_frames()
            self.opt_cont = tk.Frame(self.master, bg=cname)
            self.opt_cont.grid(row=5, column=4, columnspan=4, sticky='WE')

            tk.Label(self.opt_cont, text='Side:', bg=cname).grid(row=5, column=4, columnspan=2)

            tk.Label(self.opt_cont, text='Probe Width:', bg=cname).grid(row=5, column=6, columnspan=2)

            #next row

            #Adds side selection OTF
            tk.Label(self.opt_cont, text='OTF: ', width=10, bg=cname, anchor='e').grid(row=6, column=4)
            self.side_option = tk.IntVar()
            self.side_option.set(1)
            self.cont_OTF = tk.Radiobutton(self.opt_cont, variable=self.side_option, value=1, bg=cname)
            self.cont_OTF.grid(row=6, column=5)

            tk.Label(self.opt_cont, text='0.015: ', width=10, bg=cname, anchor='e').grid(row=6, column=6)
            self.width_option = tk.IntVar()
            self.width_option.set(1)
            self.cont_widA = tk.Radiobutton(self.opt_cont, variable=self.width_option, value=1, bg=cname)
            self.cont_widA.grid(row=6, column=7)

            #next row

            #Adds side selection ITF
            tk.Label(self.opt_cont, text='ITF: ', width=10, bg=cname, anchor='e').grid(row=7, column=4)
            self.cont_ITF = tk.Radiobutton(self.opt_cont, variable=self.side_option, value=2, bg=cname)
            self.cont_ITF.grid(row=7, column=5)

            tk.Label(self.opt_cont, text='0.005: ', width=10, bg=cname, anchor='e').grid(row=7, column=6)
            self.cont_widB = tk.Radiobutton(self.opt_cont, variable=self.width_option, value=2, bg=cname)
            self.cont_widB.grid(row=7, column=7)

            #next row

            tk.Label(self.opt_cont, text='0.0025: ', width=10, bg=cname, anchor='e').grid(row=8, column=6)
            self.cont_widC = tk.Radiobutton(self.opt_cont, variable=self.width_option, value=3, bg=cname)
            self.cont_widC.grid(row=8, column=7)

            #next row

            tk.Label(self.opt_cont, text='Other: ', width=10, bg=cname, anchor='e').grid(row=9, column=6)
            self.cont_widO = tk.Radiobutton(self.opt_cont, variable=self.width_option, value=4, bg=cname)
            self.cont_widO.grid(row=9, column=7)

            self.cont_widEnt = tk.Entry(self.opt_cont)
            self.cont_widEnt.grid(row=9, column=8)

        elif self.current_option.get() == 'Poloidal Profiles':

            self.delete_frames()
            self.opt_polo = tk.Frame(self.master, bg=cname)
            self.opt_polo.grid(row=5, column=4, columnspan=4, sticky='WE')

            tk.Label(self.opt_polo, text='',bg=cname, width=10).grid(row=5, column=4, columnspan=2)

            tk.Label(self.opt_polo, text='Probe Width:', bg=cname).grid(row=5, column=6, columnspan=2)

            tk.Label(self.opt_polo, text='0.015: ', width=10, bg=cname, anchor='e').grid(row=6, column=6)
            self.width_option = tk.IntVar()
            self.width_option.set(1)
            self.cont_widA = tk.Radiobutton(self.opt_polo, variable=self.width_option, value=1, bg=cname)
            self.cont_widA.grid(row=6, column=7)

            tk.Label(self.opt_polo, text='0.005: ', width=10, bg=cname, anchor='e').grid(row=7, column=6)
            self.cont_widB = tk.Radiobutton(self.opt_polo, variable=self.width_option, value=2, bg=cname)
            self.cont_widB.grid(row=7, column=7)

            tk.Label(self.opt_polo, text='0.0025: ', width=10, bg=cname, anchor='e').grid(row=8, column=6)
            self.cont_widC = tk.Radiobutton(self.opt_polo, variable=self.width_option, value=3, bg=cname)
            self.cont_widC.grid(row=8, column=7)

            tk.Label(self.opt_polo, text='Other: ', width=10, bg=cname, anchor='e').grid(row=9, column=6)
            self.cont_widO = tk.Radiobutton(self.opt_polo, variable=self.width_option, value=4, bg=cname)
            self.cont_widO.grid(row=9, column=7)

            self.cont_widEnt = tk.Entry(self.opt_polo)
            self.cont_widEnt.grid(row=9, column=8)

        #The temperature is selected this is what happens
        elif self.current_option.get() == 'Temperature':

            #creates new frame for the temperature selection
            self.delete_frames()
            self.opt_Temp = tk.Frame(self.master, bg=cname)
            self.opt_Temp.grid(row=5, column=4, columnspan=4, sticky='WE')

            #Adds Entry box for the ylim
            tk.Label(self.opt_Temp, text='y Position: ', bg=cname).grid(row=5, column=4)
            self.temp_ylim = tk.Entry(self.opt_Temp, width=col2_width)
            self.temp_ylim.grid(row=5, column=5, padx=padx, pady=pady)

    def delete_frames(self):
        try:
            self.opt_frame.grid_forget()
        except:
            pass

        try:
            self.opt_Temp.grid_forget()
        except:
            pass

        try:
            self.opt_polo.grid_forget()
        except:
            pass

        try:
            self.opt_cont.grid_forget()
        except:
            pass

        try:
            self.opt_centline.grid_forget()
        except:
            pass

    def browse_netcdf(self):
        """
        Function for Browse button to get netcdf file.
        """

        #adds the netcdf file
        self.message_box.insert(tk.END, 'Loading...\n')
        root = tk.Tk(); root.withdraw()
        netcdf_path = tk.filedialog.askopenfilename(filetypes=(('NetCDF files', '*.nc'),))
        self.netcdf_entry.delete(0, tk.END)
        self.netcdf_entry.insert(0, netcdf_path)
        self.dl = limpt.LimPlots(self.netcdf_entry.get())
        self.message_box.insert(tk.END, 'Loaded file: {}\n'.format(netcdf_path.split('/')[-1]))

        cp_path = netcdf_path.split('.nc')[0] + '.collector_probe'
        try:
            f = open(cp_path, 'r')
            self.cp_entry.delete(0, tk.END)
            self.cp_entry.insert(0, cp_path)
            self.message_box.insert(tk.END, 'Loaded file: {}\n'.format(cp_path.split('/')[-1]))
            self.cp_path = cp_path
        except:
            pass

        #adds the dat file
        dat_path = netcdf_path.split('.nc')[0] + '.dat'
        try:
            f = open(dat_path, 'r')
            self.dat_entry.delete(0, tk.END)
            self.dat_entry.insert(0, dat_path)
            self.message_box.insert(tk.END, 'Loaded file: {}\n'.format(dat_path.split('/')[-1]))
            self.dat_path = dat_path
        except:
            pass

        #adds the lim file
        lim_path = netcdf_path.split('.nc')[0] + '.lim'
        try:
            f = open(lim_path, 'r')
            self.lim_entry.delete(0,tk.END)
            self.lim_entry.insert(0, lim_path)
            self.message_box.insert(tk.END, 'Loaded file: {}\n'.format(lim_path.split('/')[-1]))
            self.lim_path = lim_path
        except:
            pass

    def quit_command(self):

        self.master.quit()
        self.master.destroy()


def main():

    root = tk.Tk()
    window = Window(root)
    window.mainloop()

if __name__ == '__main__':
    main()