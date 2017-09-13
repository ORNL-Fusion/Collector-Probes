# Collector-Probes

These are files related to the creation and accessing of the deposition probe MDSplus tree. The highest level program is 
ProbeClass.py. The main function used is "get\_multiple". This function will return a list of up three probes (an A, B, 
and/or C), ONLY IF RBS data is available. Each object in the list will be of "Probe" class. The Probe class includes two dictionaries: one with data from
the R2D2 server, and another with data analyzed using EFIT on atlas. The parameters passed to get\_multiple are:

ProbeClass.get\_multiple(aNumber=None       --> A probe number 
                         bNumber=None       --> B probe number
                         cNumber=None       --> C probe number
                         MDStunnel=False    --> Use True is accessing remotely outside the DIII-D network.
                         startTime=2500     --> The start time of the range that R - Rsep will be averaged for.
                         endTime=5000       --> The end time for the range.
                         step=500           --> Time step for the time range above.
                         efitTree='EFIT01'  --> Which EFIT tree is used for the R - Rsep calculations.
                         toHDF5=False       --> Save data to HDF5 file format.)
                       
To use these functions import ProbeClass into your python2 session, then call the function. An example session is as follows:

$ import ProbeClass as Probe
$ probeList = Probe.get\_multiple(aNumber=2, bNumber=2, cNumber=2, MDStunnel=True)
SSH into R2D2. Press enter when ready...

_In a new terminal ssh linking r2d2 to your localhost:_
$ ssh -Y -p 2039 -L 8000:r2d2.gat.com:8000 username@cybele.gat.com
_Then login_

_Back in the python terminal press enter:_
AU Run: 1
AU Run: 2
...
CD Run: 21
CD Run: 22
Data for run 22 not available. _Note: This is OK._
SSH into Atlas. Press enter when ready...

_Back to the other terminal logout and then ssh into atlas:_
$ ssh -Y -p 2039 -L 8000:atlas.gat.com:8000 username@cybele.gat.com

_Back in the python terminal press enter:_
Analyzing AU2 Data...
Shot: 167196
Time: 2500

Shot: 167196
Time: 3000

...

_Now all the dictionaries in the Probe class objects are filled._

$ probeList[0].letter
'A'
$ probeList[1].letter
'B'
$ probeList[2].letter
'C'
$ probeList[0].number
2
$ probeList[0].plot\_norm()
_Plots the W areal density vs. R-Rsep._
$ probeList[0].plot\_omp
