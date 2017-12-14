import get_lp as lp
import shawn_model as model
import compareLPtoTS_psiN as comp

import numpy as np
import matplotlib.pyplot as plt
import EFIT.load_gfile_d3d as loadg
from scipy import interpolate
import MDSplus as mds


lp_dict = lp.get_dict_of_lps(167192)
plunge_dict = comp.runScript()

# Average between 1900 and 2000 for the divertor LPs since this is when the
# plunging LP was in for.
psin_list = []
te_list = []
for lp_name in lp_dict.keys():
    print "Averaging " + lp_name + "..."
    lp = lp_dict[lp_name]
    times = lp["time"]
    temp = lp["temp"]
    psin = lp["psin"]

    # Find what index 1900 and 2000 is at.
    time1 = 1900
    while(True):
        try:
            index1, = np.where(times==time1)[0]
            break
        except:
            time1 = time1 + 1

    time2 = 2000
    while(True):
        try:
            index2, = np.where(times==time2)[0]
            break
        except:
            time2 = time2 + 1


    # Find the average temp in this range.
    avgTemp = np.mean(temp[index1:index2])
    avgPsin = np.mean(psin[index1:index2])

    # Then add onto lists for (R,Te).
    psin_list.append(avgPsin)
    te_list.append(avgTemp)

plungeLPpsin = plunge_dict["LPpsiNs"]
plungeLPtemps = plunge_dict["LPtemps"][:61]

plt.plot(psin_list, te_list, ".", label="Divertor")
plt.plot(plungeLPpsin, plungeLPtemps, label="Plunging")
plt.title("Divertor vs. Plunging LPs Shot #167192")
plt.xlabel(r"$\psi_N$")
plt.ylabel(r"$\mathrm{T_e\ (eV)}$")
plt.legend()
plt.text(1,1,"Divertor avg. from 1900-2000 ms")

plt.show()
