import sys
sys.path.append("/home/shawn/DIII-D/ORNL-Fusion/Collector-Probes")
import get_ts as ts
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d


ts192 = ts.run_script(167192, "core", tmin=1900, tmax=3000, tstep=150)
ts193 = ts.run_script(167193, "core", tmin=1900, tmax=3000, tstep=150)
ts194 = ts.run_script(167194, "core", tmin=1900, tmax=3000, tstep=150)
ts195 = ts.run_script(167195, "core", tmin=1900, tmax=3000, tstep=150)

all_psin = np.array([])
all_Te = np.array([])
all_Te_err = np.array([])
for ts_dict in [ts192, ts193, ts194, ts195]:
    psin = ts_dict["psins"]["avg_psins"][:-1]
    all_psin = np.append(all_psin, psin)
    Te = ts_dict["psins"]["avg_Tes"][:-1]
    Te_err = ts_dict["psins"]["avg_Tes_err"][:-1]
    all_Te = np.append(all_Te, Te)
    all_Te_err = np.append(all_Te_err, Te_err)
    shot = ts_dict["shot"]
    #plt.plot(psin, Te, label=str(shot))

# Sort all_psin and all_Te.
idx = np.argsort(all_psin)
all_psin = all_psin[idx]
all_Te = all_Te[idx]
all_Te_err = all_Te_err[idx]

# Smooth the Te data.
#smooth_Te = savgol_filter(all_Te, 3, 2)
#smooth_Te = interp1d(all_psin, all_Te, kind='cubic')

plt.errorbar(all_psin, all_Te, yerr=all_Te_err, label="167192-167195 (1900-3000ms)")
#plt.axis([0.9, 1.5, 0, 50])
plt.xlabel("Psin")
plt.ylabel("Te (eV)")
plt.legend()
plt.show()
