# OEDGE GUI

This is a GUI written by Shawn Zamperini to plot the results from a DIVIMP run. To run the GUI, simply run oedge_plots_gui from the command line:

```
python3 oedge_plots_gui.py
```

The interface is fairly self-explanatory. Just click browse to load in the path to the NetCDF file from the OEDGE output. The GUI will automatically look in the same folder as the NetCDF file for the dat and collectorprobe files of the same run name. 

Required packages: tkinter, netCDF4, pandas, numpy, matplotlib.

Report bugs or feature requests to zamp@utk.edu.