import numpy as np
import matplotlib.pyplot as plt
from ThomsonClass import ThomsonClass
from scipy import interpolate
from matplotlib.lines import Line2D


# Only things you need to change are here. This script has only been tested on 178351.
shot = 178351
time = 1700  # 1580, 1740, 1880, 3000
cm_flux = 4  # Which cm flux surface to plot.

# Some constants. Probe tip.
r_probe =  1.4895
z_probe = -1.205

# Probe bottom.
#r_probe =  1.4895
#z_probe = -1.245

# Load gfile, scavenging the code from ThomsonClass.
ts = ThomsonClass(shot, "core")
gfile = ts.load_gfile_mds(shot, time, tree="EFIT01")

# Get the relevant arrays from the gfile.
psin = gfile["psiRZn"]
R, Z = np.meshgrid(gfile["R"], gfile["Z"])
Z_axis = gfile['ZmAxis']
R_axis = gfile['RmAxis']
Zes = np.copy(gfile['lcfs'][:, 1][13:-17])
Res = np.copy(gfile['lcfs'][:, 0][13:-17])

# Find the nearest psin to the probe coordinates.
r_idx = np.where(np.abs(gfile["R"]-r_probe)==np.abs(gfile["R"]-r_probe).min())[0][0]
z_idx = np.where(np.abs(gfile["Z"]-z_probe)==np.abs(gfile["Z"]-z_probe).min())[0][0]

# The below is used to find out what psin the probe is on. Use interpolation
# functions on the grid supplied from the gfile. Limited to the right half of the
# plasma otherwise the interpolations get wonky from trying to interpolate
# too much.
#Rs_trunc = R > R_axis
Rs_trunc = R > 1.30
Rs_right_side = R > 1.2
f_psin = interpolate.Rbf(R[Rs_right_side], Z[Rs_right_side], psin[Rs_right_side])
f_Rs   = interpolate.interp1d(Zes, Res, assume_sorted=False)
f_Zs   = interpolate.interp1d(Res[Zes<Z_axis], Zes[Zes<Z_axis], assume_sorted=False)
psin_probe = f_psin(r_probe, z_probe)
print("Probe on psin = {:.3f}.".format(psin_probe))
f_R_right = interpolate.Rbf(gfile['psiRZn'][Rs_trunc], Z[Rs_trunc], R[Rs_trunc], epsilon=0.00001)
rsep_omp = f_Rs(Z_axis)
zsep_bot = f_Zs(r_probe)
cm_flux_surf = f_psin(rsep_omp + cm_flux / 100, Z_axis)
rsep_mimes = f_Rs(-0.188)
r_mimes = f_R_right(psin_probe, -0.188)
rminrsep_omp_dimes = f_R_right(psin_probe, Z_axis) - rsep_omp

# Print some info.
print("MiMES inserted to (R, Z) = ({:.3f}, {:.3f})".format(-0.188, r_mimes))
print("First MiMES contact at is {:.3f} m, or {:.3f} cm from separatrix".format(r_mimes, (r_mimes-rsep_mimes)*100))
print("DiMES tip is on the {:.2f} cm flux surface.".format(rminrsep_omp_dimes * 100))
print("GAPBOT = {:.2f} cm".format((zsep_bot - (-1.250))*100))

# Plotting commands.
fig, ax = plt.subplots(figsize=(4,7))
cont = ax.contour(R, Z, psin, levels=[psin_probe], colors='r')
cont2 = ax.contour(R, Z, psin, levels=[cm_flux_surf], colors='b')

# A line from the floor to show GAPBOT.
ax.plot((r_probe, r_probe), (-1.250, zsep_bot), 'r-')

ax.plot(r_probe, z_probe, 'k.', ms=5)
ax.plot([r_probe, r_probe], [z_probe, -1.25], 'k')
ax.plot(r_mimes, -0.118, 'k.', ms=5)
ax.plot([r_mimes, 2.35], [-0.188, -0.188], 'k')
ax.plot(gfile['lcfs'][:,0], gfile['lcfs'][:,1], 'k')
ax.plot(gfile['wall'][:,0], gfile['wall'][:,1], 'k')
ax.text(0.5, 0.5, "DiMES on {:.2f} cm\n flux surface".format(rminrsep_omp_dimes * 100), transform=ax.transAxes, horizontalalignment="center", fontsize=12)
ax.axis("equal")
ax.set_title("{} {}ms".format(shot, time))
ax.set_ylim([-1.5, 1.5])
custom_lines = [Line2D([0], [0], color='r', lw=2),
                Line2D([0], [0], color='b', lw=2)]
ax.legend(custom_lines, ["DiMES Tip", "{}cm flux surface".format(cm_flux)])
fig.tight_layout()
fig.show()
