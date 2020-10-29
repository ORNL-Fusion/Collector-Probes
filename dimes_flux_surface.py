import numpy as np
import matplotlib.pyplot as plt
from ThomsonClass import ThomsonClass
from scipy import interpolate
from matplotlib.lines import Line2D


# Some options we've been looking at.
# 178351: 1580, 1740, 1880, 3000
# 178341: 1700
# 969002:

# Only things you need to change are here. This script has only been tested on 178351.
shot        = 969002
time        = 2900
tree        = "EFITRT1"  # EFIT01 or EFITRT1
cm_flux     = np.arange(1, 8) # Which cm flux surface to plot.
extra_mimes = 2  # How many extra cm to insert MiMES from the tube where it intercepts the tip of DiMES.

# Some constants. DiMES probe tip.
r_probe =  1.4895
z_probe = -1.205

# Probe bottom.
r_probe_bot =  1.4895
z_probe_bot = -1.245

# Length of MiMES.
mimes_len = 10  # cm

# Load gfile, scavenging the code from ThomsonClass.
ts = ThomsonClass(shot, "core")
gfile = ts.load_gfile_mds(shot, time, tree=tree)

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
psin_probe_bot = f_psin(r_probe_bot, z_probe_bot)
f_R_right = interpolate.Rbf(gfile['psiRZn'][Rs_trunc], Z[Rs_trunc], R[Rs_trunc], epsilon=0.00001)
rsep_omp = f_Rs(Z_axis)
zsep_bot = f_Zs(r_probe)
cm_flux_surf = f_psin(rsep_omp + cm_flux / 100, np.full(len(cm_flux), Z_axis))
rsep_mimes = f_Rs(-0.188)
r_mimes = f_R_right(psin_probe, -0.188) - extra_mimes / 100
rminrsep_omp_dimes = f_R_right(psin_probe, Z_axis) - rsep_omp
rminrsep_omp_dimes_bot = f_R_right(psin_probe_bot, Z_axis) - rsep_omp
psin_mimes_top = f_psin(r_mimes, -0.188)
psin_mimes_bot = f_psin(r_mimes + mimes_len / 100, -0.188)

# Printout of information.
print("\nDiMES")
print("  Psin:     {:.3f} - {:.3f}".format(psin_probe, psin_probe_bot))
print("  Surfaces: {:.2f}  - {:.2f}".format(rminrsep_omp_dimes*100, rminrsep_omp_dimes_bot*100))
print("\nMiMES")
print("  Centimeters inserted past surface that hits DiMES tip: {:.2f}".format(extra_mimes))
print("  Insert MiMES to (R, Z): ({:.3f}, {:.3f})".format(r_mimes, -0.188))
print("  Psin:     {:.3f} - {:.3f}".format(psin_mimes_top, psin_mimes_bot))
print("  Surfaces: {:.2f}  - {:.2f}".format((r_mimes-rsep_mimes)*100, (r_mimes-rsep_mimes)*100+mimes_len))
print("\nGAPBOT = {:.2f} cm".format((zsep_bot - (-1.250))*100))


# This will set up the coordinates needed to do a fill between a vertical line that
# extends from the inner wall to some arbitrary distance outward (towards lower R)
# so that we can cover flux tubes that go past the inner wall. Covering the rest
# of the flux tubes that go outside the vessel can be done with a normal fill between.

#    inner_wall_R = gfile['wall'][:,0][0]
#    outer_wall_R = gfile['wall'][:,0][90]
#    right_bound_Z = np.append(gfile['wall'][:,1][10:14][::-1], [-2])
#    right_bound_R = np.append(gfile['wall'][:,0][10:14][::-1], [inner_wall_R])
#    right_bound_Z = np.insert(right_bound_Z, 0, 2)
#    right_bound_R = np.insert(right_bound_R, 0, right_bound_R[-1])
inner_wall_R = 1.0160000324249268
outer_wall_R = 2.3515799045562744
right_bound_Z = np.array([2.0, 1.21700001, 1.21700001, 1.16499996, 1.14699996,-2.0])
right_bound_R = np.array([1.01600003, 1.02900004, 1.00100005, 1.01199996, 1.01600003, 1.01600003])

# Plotting commands.
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,7), sharex=True)
ax1.contourf(R, Z, psin, levels=[psin_probe, psin_probe_bot], colors='r')
ax1.contour(R, Z, psin, levels=cm_flux_surf, colors='k', alpha=0.25)
ax2.contourf(R, Z, psin, levels=[psin_mimes_top, psin_mimes_bot], colors='b')
ax2.contour(R, Z, psin, levels=cm_flux_surf, colors='k', alpha=0.25)

# A line from the floor to show GAPBOT.
ax1.plot((r_probe, r_probe), (-1.250, zsep_bot), 'r-')
ax2.plot((r_probe, r_probe), (-1.250, zsep_bot), 'r-')

ax1.plot(r_probe, z_probe, 'k.', ms=5)
ax1.plot([r_probe, r_probe], [z_probe, -1.25], 'k')
ax1.plot(r_mimes, -0.188, 'k.', ms=5)
ax1.plot([r_mimes, 2.35], [-0.188, -0.188], 'k')
ax1.plot(gfile['lcfs'][:,0], gfile['lcfs'][:,1], 'k')
ax1.plot(gfile['wall'][:,0], gfile['wall'][:,1], 'k', zorder=100)
ax1.axis("equal")
ax1.set_title("DiMES", color='r', fontsize=16)
ax1.fill_between(gfile['wall'][:,0], gfile['wall'][:,1], gfile['wall'][:,1]*100, color='w', zorder=99)
ax1.fill_betweenx(right_bound_Z, right_bound_R, right_bound_R - 0.2, color='w', zorder=99)
ax1.fill_betweenx([-2, 2], [outer_wall_R, outer_wall_R], [outer_wall_R+2, outer_wall_R+2], color='w', zorder=99)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.get_xaxis().set_visible(False)
ax1.get_yaxis().set_visible(False)
ax1.set_ylim(-1, 1)

ax2.plot(r_probe, z_probe, 'k.', ms=5)
ax2.plot([r_probe, r_probe], [z_probe, -1.25], 'k')
ax2.plot(r_mimes, -0.188, 'k.', ms=5)
ax2.plot([r_mimes, 2.35], [-0.188, -0.188], 'k')
ax2.plot(gfile['lcfs'][:,0], gfile['lcfs'][:,1], 'k')
ax2.plot(gfile['wall'][:,0], gfile['wall'][:,1], 'k', zorder=100)
ax2.axis("equal")
ax2.set_title("MiMES", color='b', fontsize=16)
ax2.set_ylim([-1.5, 1.5])
ax2.fill_between(gfile['wall'][:,0], gfile['wall'][:,1], gfile['wall'][:,1]*100, color='w', zorder=99)
ax2.fill_betweenx(right_bound_Z, right_bound_R, right_bound_R - 0.2, color='w', zorder=99)
ax2.fill_betweenx([-2, 2], [outer_wall_R, outer_wall_R], [outer_wall_R+2, outer_wall_R+2], color='w', zorder=99)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.get_xaxis().set_visible(False)
ax2.get_yaxis().set_visible(False)

plt.suptitle("{} {}ms".format(shot, time), fontsize=16)
fig.tight_layout()
fig.show()
plt.savefig("fig.png")
