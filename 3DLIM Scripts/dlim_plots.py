import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

class dlim_plots:

    def __init__(self, netcdf_path):

        self.nc = netCDF4.Dataset(netcdf_path)
        self.ncpath = netcdf_path

        self.datfile = None

    def plot_2d(self, dataname, ylim=500):
        x = self.nc.variables['XOUTS'][:]

        if dataname == 'CTEMBS':
            y = self.nc.variables['CTEMBS'][ylim]

        plt.plot(x,y)
        plt.show()




