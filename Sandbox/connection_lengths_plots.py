import openpyxl as xl
import numpy as np
import matplotlib.pyplot as plt
import EFIT.load_gfile_d3d as loadg
import scipy.interpolate as scinter
import MDSplus as mds

conn_wb = xl.load_workbook("conn_lengths.xlsx", data_only=True)
samp_wb = xl.load_workbook("link_to_lp_with_fits.xlsx", data_only=True)

conn_sheet = conn_wb.get_sheet_by_name("Sheet1")
samp_sheet = samp_wb.get_sheet_by_name("Fit Te's")

def returnArray(sheet, lowRange, highRange):
    cells = sheet[lowRange:highRange]
    cells = np.transpose(cells)
    cells = np.reshape(cells, cells.size)
    values = np.array([cell.value for cell in cells])
    return values

# Get the sampling length data. This is in R-Rsep vs. L.
samp_R        = returnArray(samp_sheet, "M2", "M348")
samp_L_A        = returnArray(samp_sheet, "Q2", "Q348")
samp_L_B        = returnArray(samp_sheet, "R2", "R348")
samp_L_C        = returnArray(samp_sheet, "S2", "S348")

# Get the connection length data. This is in R-Rsep_omp vs. L.
conn_RminRsep_omp_A = returnArray(conn_sheet, "A17", "A113")
conn_odf_A          = returnArray(conn_sheet, "B17", "B113")
conn_idf_A          = returnArray(conn_sheet, "C17", "C113")
conn_RminRsep_omp_B = returnArray(conn_sheet, "E17", "E113")
conn_odf_B          = returnArray(conn_sheet, "F17", "F113")
conn_idf_B          = returnArray(conn_sheet, "G17", "G113")
conn_RminRsep_omp_C = returnArray(conn_sheet, "I17", "I113")
conn_odf_C          = returnArray(conn_sheet, "J17", "J113")
conn_idf_C          = returnArray(conn_sheet, "K17", "K113")



# Since the data from D. Elder is from the omp, I'll have to translate the
# sampling length data to the omp. This just means converting an R value to
# an R-Rsep_omp value. Easy peasy.
conn = mds.Connection("localhost")
print "Loading gfile..."
parmDICT = loadg.read_g_file_mds(167196, 3000, connection=conn, write2file=False, tree="EFIT01")
Rs, Zs = np.meshgrid(parmDICT['R'], parmDICT['Z'])
Z_axis = parmDICT['ZmAxis']
R_axis = parmDICT['RmAxis']
Zes = np.copy(parmDICT['lcfs'][:, 1][13:-12])
Res = np.copy(parmDICT['lcfs'][:, 0][13:-12])
f_Rs = scinter.interp1d(Zes, Res, assume_sorted=False)

# Functions for inteprolation.
Rs_trunc = Rs > R_axis
f_psiN = scinter.Rbf(Rs[Rs_trunc],
                     Zs[Rs_trunc],
                     parmDICT['psiRZn'][Rs_trunc], function='linear')
f_Romp = scinter.interp2d(parmDICT['psiRZn'][Rs_trunc], Zs[Rs_trunc], Rs[Rs_trunc])

psin_A = np.array([])
psin_B = np.array([])
psin_C = np.array([])
for R in samp_R:
    # Get the psin value. Use m instead of cm.
    tmp_psin_A = f_psiN(R / 100.0, -0.18)
    tmp_psin_B = f_psiN(R / 100.0, -0.1546)
    tmp_psin_C = f_psiN(R / 100.0, -0.2054)
    psin_A = np.append(psin_A, tmp_psin_A)
    psin_B = np.append(psin_B, tmp_psin_B)
    psin_C = np.append(psin_C, tmp_psin_C)

# Now go from psin and the known Z of the omp (Zaxis) to R_omp.
R_omp_A = np.array([])
R_omp_B = np.array([])
R_omp_C = np.array([])
for psin in psin_A:
    tmp_R_omp = f_Romp(psin, Z_axis)
    R_omp_A = np.append(R_omp_A, tmp_R_omp)
for psin in psin_B:
    tmp_R_omp = f_Romp(psin, Z_axis)
    R_omp_B = np.append(R_omp_B, tmp_R_omp)
for psin in psin_C:
    tmp_R_omp = f_Romp(psin, Z_axis)
    R_omp_C = np.append(R_omp_C, tmp_R_omp)

# Get whatever the Rsep at the omp is.
R_omp = f_Rs(Z_axis)

# Finally get R-Rsep_omp.
RminRsep_omp_A = (R_omp_A - R_omp) * 100.0
RminRsep_omp_B = (R_omp_B - R_omp) * 100.0
RminRsep_omp_C = (R_omp_C - R_omp) * 100.0

plt.rcParams.update({'font.size': 22})
plt.semilogy(RminRsep_omp_A, samp_L_A*2.0, label="Sampling Length A x 2")
plt.semilogy(RminRsep_omp_B, samp_L_B*2.0, label="Sampling Length B x 2")
plt.semilogy(RminRsep_omp_C, samp_L_C*2.0, label="Sampling Length C x 2")
#plt.semilogy(conn_RminRsep_omp_A, conn_odf_A, label="Outer Divertor Facing A")
#plt.semilogy(conn_RminRsep_omp_A, conn_idf_A, label="Inner Divertor Facing A")
plt.semilogy(conn_RminRsep_omp_B, conn_odf_B, label="Outer Divertor Connection Length")
plt.semilogy(conn_RminRsep_omp_B, conn_idf_B, label="Inner Divertor Connection Length")
#plt.semilogy(conn_RminRsep_omp_C, conn_odf_C, label="Outer Divertor Facing C")
#plt.semilogy(conn_RminRsep_omp_C, conn_idf_C, label="Inner Divertor Facing C")
plt.legend(prop={"size":18})
plt.xlabel("R-Rsep omp (cm)")
plt.ylabel("Length (cm)")
plt.title("Connection vs. Sampling Lengths")
plt.show()
