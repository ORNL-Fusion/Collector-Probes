# Collector-Probes

These are files related to the creation and accessing of the deposition probe MDSplus tree. Everything is written in python2. 

**ProbeClass.py** - Highest level program that incorporates the others. Most users will use this script to access the data.  
**meas\_locations.py** - Accounts for the 13 degree angle the probes are inserted in by. Feeds into get\_Rsep.py.   
**pull\_data\_dp.py** - Pulls RBS data from the MDS+ tree on the R2D2 server.  
**get\_Rsep.py** - Uses meas\_locations and load\_gfile\_d3d from the EFIT github repository to find the average R-Rsep and
R\_omp - Rsep\_omp using EFIT.   

The highest level program is ProbeClass.py, and the main function used is "get\_multiple". This function will return a list of up three
probes (an A, B, and/or C), ONLY IF RBS data is available. Each object in the list will be of "Probe" class. The Probe class includes 
two dictionaries: one with data from the R2D2 server, and another with data analyzed using EFIT on atlas. The parameters passed to
get\_multiple are:
  
```python
ProbeClass.get_multiple(  
      aNumber=None       # A probe number  
      bNumber=None       # B probe number  
      cNumber=None       # C probe number  
      MDStunnel=False    # Use True is accessing remotely outside the DIII-D network.  
      startTime=2500     # The start time of the range that R - Rsep will be averaged for.  
      endTime=5000       # The end time for the range.  
      step=500           # Time step for the time range above.  
      efitTree='EFIT01'  # Which EFIT tree is used for the R - Rsep calculations.  
      toHDF5=False       # Save data to HDF5 file format.
      )
```  
  
To use these functions import ProbeClass into your python2 session, then call the function. An example session is as follows:   
```
$ import ProbeClass as Probe  
$ probeList = Probe.get_multiple(aNumber=2, bNumber=2, cNumber=2, MDStunnel=True)  
SSH into R2D2. Press enter when ready...  
```

_In a new terminal, ssh linking r2d2 to your localhost:_  
```
$ ssh -Y -p 2039 -L 8000:r2d2.gat.com:8000 username@cybele.gat.com
```
_Then login._  

_Back in the python terminal press enter:_ 
```
AU Run: 1  
AU Run: 2  
...  
CD Run: 21  
CD Run: 22  
Data for run 22 not available.
SSH into Atlas. Press enter when ready...  
```
  
_Back to the other terminal logout and then ssh into atlas:_ 
```
$ ssh -Y -p 2039 -L 8000:atlas.gat.com:8000 username@cybele.gat.com  
```
  
_Back in the python terminal press enter:_  
```
Analyzing AU2 Data...  
Shot: 167196  
Time: 2500  
  
Shot: 167196  
Time: 3000  
  
...  
``` 
_Now all the dictionaries in the Probe class objects are filled._  
```
$ probeList[0].letter  
'A'  
$ probeList[1].letter  
'B'  
$ probeList[2].letter  
'C'  
$ probeList[0].number  
2  
$ probeList[0].plot_norm()
```

_Plots the W areal density vs. R-Rsep._  
  
```
$ probeList[0].plot_omp  
```
_Plots the W areal density vs. R\_omp-Rsep\_omp._  
```
$ dict1 = probeList[0].r2d2DICT  
```
_Return a dictionary of the data pulled from R2D2. Check the r2d2 function docstring._  
```
$ dict2 = probeList[0].atlas  
```
_Return a dictionary of the data analyzed using EFIT. Check the atlas function docstring._  
```
$ w_areal = probeList[0].r2d2['w_areal_U']  
$ rminrsep = probeList[0].atlas['rminrsep_U']  
```
_You can use these two lists to plot W areal density vs. R-Rsep for example._  
```
$ probeList[0].to_matlab()
```
_Save the data to a matlab file. Useful for the curve fitting tool._
