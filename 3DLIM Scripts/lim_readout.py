from Readout import Readout
from tkinter import filedialog, Tk


# Ask user for location to netCDF file.
root = Tk(); root.withdraw()
netcdf_path = filedialog.askopenfilename(filetypes=(('NetCDF files', '*.nc'),))

if netcdf_path == ():
    print("Using testfile...")
    netcdf_path = None
    dat_path = None
else:
    # Ask for .dat file as well.
    dat_path = filedialog.askopenfilename(filetypes=(('.dat files', '*.dat'),))

# Option if dat file not included.
if dat_path == ():
    dat_path = None

# Variables for plots.
probe_width = 0.015
rad_cutoff = 0.05

# Controlling routine for Readout.
grid = Readout(netcdf_file=netcdf_path, dat_file=dat_path)
grid.print_readout()
grid.centerline(0)
grid.avg_imp_vely(1)
grid.te_contour(2)
grid.deposition_contour(3, probe_width=probe_width, rad_cutoff=rad_cutoff, side='ITF')
grid.deposition_contour(4, probe_width=probe_width, rad_cutoff=rad_cutoff, side='OTF')
grid.ne_contour(5)
grid.avg_pol_profiles(6, probe_width=probe_width, rad_cutoff=rad_cutoff)
grid.show_fig()
