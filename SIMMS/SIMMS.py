import numpy as np
import scipy.interpolate as sinter
import pandas as pd


def doALL(EFfileNAM="AU02_dict.csv", RBSfileNAM="A02_RminRsep_data.csv", FOLDpath='/.',
          simmBG=0.0435, HDFdump=0):

    efFOLDpath = FOLDpath+'CP_LAMS_Dictionary/'
    efDF = getCPenrichment(fileNAM=EFfileNAM, datFOLDpath=efFOLDpath)

    rbsFOLDpath = FOLDpath+'RminRsep_RBS_Dictionary/'
    rbsDIC = makRBS_DIC(RBSfileNAM, datFOLDpath=rbsFOLDpath)

    simmDF = mak2sourceSIMMS(efDF, BGfrac=simmBG, HDFdump=0)

    # make some names
    rbsDICnam = efDF.DFname[1]+"df"
    rminrsepNAM = "rminrsep_"+efDF.DFname[1]
    arealNAM = "w_areal_"+efDF.DFname[1]
    err_arealNAM = "w_areal_err_"+efDF.DFname[1]

    # RBS probe distance
    rbs_dprobe_mm = np.array(rbsDIC[rbsDICnam].index[::-1])*10.

    # individual SIMM interpolation
    # -----------------------------
    simm_dprobe = np.around(np.array(simmDF.index), decimals=2)

    ftemp = sinter.interp1d(simm_dprobe, simmDF['Ffloor_8082'], kind='cubic',
                            fill_value=0.0, assume_sorted=True)
    Ffloor_8082_interp = ftemp(rbs_dprobe_mm)
    ftemp = sinter.interp1d(simm_dprobe, simmDF['err_Ffloor_8082'], kind='cubic',
                            fill_value=0.0, assume_sorted=True)
    err_Ffloor_8082_interp = ftemp(rbs_dprobe_mm)
    ftemp = sinter.interp1d(simm_dprobe, simmDF['Ffloor_8382'], kind='cubic',
                            fill_value=0.0, assume_sorted=True)
    Ffloor_8382_interp = ftemp(rbs_dprobe_mm)
    ftemp = sinter.interp1d(simm_dprobe, simmDF['err_Ffloor_8382'], kind='cubic',
                            fill_value=0.0, assume_sorted=True)
    err_Ffloor_8382_interp = ftemp(rbs_dprobe_mm)
    ftemp = sinter.interp1d(simm_dprobe, simmDF['Ffloor_8482'], kind='cubic',
                            fill_value=0.0, assume_sorted=True)
    Ffloor_8482_interp = ftemp(rbs_dprobe_mm)
    ftemp = sinter.interp1d(simm_dprobe, simmDF['err_Ffloor_8482'], kind='cubic',
                            fill_value=0.0, assume_sorted=True)
    err_Ffloor_8482_interp = ftemp(rbs_dprobe_mm)
    ftemp = sinter.interp1d(simm_dprobe, simmDF['Ffloor_8682'], kind='cubic',
                            fill_value=0.0, assume_sorted=True)
    Ffloor_8682_interp = ftemp(rbs_dprobe_mm)
    ftemp = sinter.interp1d(simm_dprobe, simmDF['err_Ffloor_8682'], kind='cubic',
                            fill_value=0.0, assume_sorted=True)
    err_Ffloor_8682_interp = ftemp(rbs_dprobe_mm)

    ftemp = sinter.interp1d(simm_dprobe, simmDF['Fshelf_8082'], kind='cubic',
                            fill_value=0.0, assume_sorted=True)
    Fshelf_8082_interp = ftemp(rbs_dprobe_mm)
    ftemp = sinter.interp1d(simm_dprobe, simmDF['err_Fshelf_8082'], kind='cubic',
                            fill_value=0.0, assume_sorted=True)
    err_Fshelf_8082_interp = ftemp(rbs_dprobe_mm)
    ftemp = sinter.interp1d(simm_dprobe, simmDF['Fshelf_8382'], kind='cubic',
                            fill_value=0.0, assume_sorted=True)
    Fshelf_8382_interp = ftemp(rbs_dprobe_mm)
    ftemp = sinter.interp1d(simm_dprobe, simmDF['err_Fshelf_8382'], kind='cubic',
                            fill_value=0.0, assume_sorted=True)
    err_Fshelf_8382_interp = ftemp(rbs_dprobe_mm)
    ftemp = sinter.interp1d(simm_dprobe, simmDF['Fshelf_8482'], kind='cubic',
                            fill_value=0.0, assume_sorted=True)
    Fshelf_8482_interp = ftemp(rbs_dprobe_mm)
    ftemp = sinter.interp1d(simm_dprobe, simmDF['err_Fshelf_8482'], kind='cubic',
                            fill_value=0.0, assume_sorted=True)
    err_Fshelf_8482_interp = ftemp(rbs_dprobe_mm)
    ftemp = sinter.interp1d(simm_dprobe, simmDF['Fshelf_8682'], kind='cubic',
                            fill_value=0.0, assume_sorted=True)
    Fshelf_8682_interp = ftemp(rbs_dprobe_mm)
    ftemp = sinter.interp1d(simm_dprobe, simmDF['err_Fshelf_8682'], kind='cubic',
                            fill_value=0.0, assume_sorted=True)
    err_Fshelf_8682_interp = ftemp(rbs_dprobe_mm)
    # -----------------------------------------------

    frac_areal_rminrsep = np.array(rbsDIC[rbsDICnam][rminrsepNAM][::-1])
    rbsTEMP = rbsDIC[rbsDICnam][arealNAM].values[::-1]*1e15
    err_rbsTEMP = rbsDIC[rbsDICnam][err_arealNAM].values[::-1]*1e15

    ftemp = sinter.interp1d(simm_dprobe, simmDF['simm_floor_AVG'], kind='cubic',
                            fill_value=0.0, assume_sorted=True)
    Ffloor_AVG_interp = ftemp(rbs_dprobe_mm)
    ftemp = sinter.interp1d(simm_dprobe, simmDF['err_simm_floor_AVG'], kind='slinear',
                            fill_value=0.0, assume_sorted=True)
    err_Ffloor_AVG_interp = ftemp(rbs_dprobe_mm)
    ftemp = sinter.interp1d(simm_dprobe, simmDF['simm_shelf_AVG'], kind='cubic',
                            fill_value=0.0, assume_sorted=True)
    Fshelf_AVG_interp = ftemp(rbs_dprobe_mm)
    ftemp = sinter.interp1d(simm_dprobe, simmDF['err_simm_shelf_AVG'], kind='slinear',
                            fill_value=0.0, assume_sorted=True)
    err_Fshelf_AVG_interp = ftemp(rbs_dprobe_mm)

    arealFRAC_floor = Ffloor_AVG_interp * rbsTEMP
    err_arealFRAC_floor = np.abs(arealFRAC_floor * np.sqrt((err_Ffloor_AVG_interp / Ffloor_AVG_interp)**2 +
                                                           (err_rbsTEMP / rbsTEMP)**2))
    arealFRAC_shelf = Fshelf_AVG_interp * rbsTEMP
    err_arealFRAC_shelf = np.abs(arealFRAC_shelf * np.sqrt((err_Fshelf_AVG_interp / Fshelf_AVG_interp)**2 +
                                                           (err_rbsTEMP / rbsTEMP)**2))
    arealFRACdic = {
                   'arealFRAC_dprobe_mm': pd.Series(rbs_dprobe_mm, index=frac_areal_rminrsep),
                   'arealFRAC_floor': pd.Series(arealFRAC_floor, index=frac_areal_rminrsep),
                   'err_arealFRAC_floor': pd.Series(err_arealFRAC_floor, index=frac_areal_rminrsep),
                   'arealFRAC_shelf': pd.Series(arealFRAC_shelf, index=frac_areal_rminrsep),
                   'err_arealFRAC_shelf': pd.Series(err_arealFRAC_shelf, index=frac_areal_rminrsep)
                   }
    arealfrDF = pd.DataFrame(arealFRACdic)

    if HDFdump:
        HDFfilename = str(efDF.DFname)+'_arealFRAC_full_2igor.h5'

        hdf = pd.HDFStore(HDFfilename)
        hdf.put(str(simmDF.DFname), simmDF, format='table', data_columns=True)
        hdf.put(str(efDF.DFname)+"_ef", efDF, format='table', data_columns=True, append=True)
        hdf.put(str(efDF.DFname)+"_rbs", rbsDIC[rbsDICnam], format='table', data_columns=True,
                append=True)
        hdf.put(str(efDF.DFname)+"_arealF", arealfrDF, format='table', data_columns=True,
                append=True)
        hdf.close()

    return simmDF, rbsDIC, efDF, arealfrDF


def mak2sourceSIMMS(enrichDF, BGfrac=0.0435, HDFdump=0):
    # First get the source standard numbers
    STANDdic = makSTANDARDS()

    # Now what through the analysis for each ratio
    Rprob_8082 = enrichDF['EF180']/enrichDF['EF182']
    Rprob_8082 = np.nan_to_num(Rprob_8082)
    Rprob_8082[Rprob_8082 > 2.0] = 2.0
    Rprob_8082[Rprob_8082 < -2.0] = 0.0
    Rprob_8082_err = Rprob_8082 * np.sqrt((enrichDF['EF180err'] / enrichDF['EF180'])**2 +
                                          (enrichDF['EF182err'] / enrichDF['EF182'])**2)
    Rprob_8082_err = np.nan_to_num(Rprob_8082_err)
    Rprob_8082_err[Rprob_8082_err > 1.0] = 1.0
    Rprob_8082_err[Rprob_8082_err < -1.0] = 1.0
    delW_probe_8082 = (Rprob_8082/STANDdic['Rnist_8082']) - 1.
    err_delW_probe_8082 = delW_probe_8082*np.sqrt((Rprob_8082_err / Rprob_8082)**2 +
                                                  (STANDdic['err_Rnist_8082'] /
                                                  STANDdic['Rnist_8082'])**2)
    Fshelf_8082 = (delW_probe_8082 - (STANDdic['delW_floor_8082'] * (1-BGfrac))) / \
                  (STANDdic['delW_shelf_8082'] - STANDdic['delW_floor_8082'])
    err_Fshelf_8082 = Fshelf_8082*np.sqrt((err_delW_probe_8082/delW_probe_8082)**2 +
                                          (STANDdic['err_delW_floor_8082'] /
                                          STANDdic['delW_floor_8082'])**2 +
                                          (STANDdic['err_delW_shelf_8082'] /
                                          STANDdic['delW_shelf_8082'])**2)
    Ffloor_8082 = 1 - Fshelf_8082 - BGfrac

    Rprob_8382 = enrichDF['EF183']/enrichDF['EF182']
    Rprob_8382 = np.nan_to_num(Rprob_8382)
    Rprob_8382[Rprob_8382 > 2.0] = 2.0
    Rprob_8382[Rprob_8382 < -2.0] = -2.0
    Rprob_8382_err = Rprob_8382 * np.sqrt((enrichDF['EF183err'] / enrichDF['EF183'])**2 +
                                          (enrichDF['EF182err'] / enrichDF['EF182'])**2)
    Rprob_8382_err = np.nan_to_num(Rprob_8382_err)
    Rprob_8382_err[Rprob_8382_err > 1.0] = 1.0
    Rprob_8382_err[Rprob_8382_err < -1.0] = 1.0
    delW_probe_8382 = (Rprob_8382/STANDdic['Rnist_8382']) - 1.
    err_delW_probe_8382 = delW_probe_8382*np.sqrt((Rprob_8382_err / Rprob_8382)**2 +
                                                  (STANDdic['err_Rnist_8382'] /
                                                  STANDdic['Rnist_8382'])**2)
    Fshelf_8382 = (delW_probe_8382 - (STANDdic['delW_floor_8382'] * (1-BGfrac))) / \
                  (STANDdic['delW_shelf_8382'] - STANDdic['delW_floor_8382'])
    err_Fshelf_8382 = Fshelf_8382*np.sqrt((err_delW_probe_8382/delW_probe_8382)**2 +
                                          (STANDdic['err_delW_floor_8382'] /
                                          STANDdic['delW_floor_8382'])**2 +
                                          (STANDdic['err_delW_shelf_8382'] /
                                          STANDdic['delW_shelf_8382'])**2)
    Ffloor_8382 = 1 - Fshelf_8382 - BGfrac

    Rprob_8482 = enrichDF['EF184']/enrichDF['EF182']
    Rprob_8482 = np.nan_to_num(Rprob_8482)
    Rprob_8482[Rprob_8482 > 2.0] = 2.0
    Rprob_8482[Rprob_8482 < -2.0] = -2.0
    Rprob_8482_err = Rprob_8482 * np.sqrt((enrichDF['EF184err'] / enrichDF['EF184'])**2 +
                                          (enrichDF['EF182err'] / enrichDF['EF182'])**2)
    Rprob_8482_err = np.nan_to_num(Rprob_8482_err)
    Rprob_8482_err[Rprob_8482_err > 1.0] = 1.0
    Rprob_8482_err[Rprob_8482_err < -1.0] = 1.0
    delW_probe_8482 = (Rprob_8482/STANDdic['Rnist_8482']) - 1.
    err_delW_probe_8482 = delW_probe_8482*np.sqrt((Rprob_8482_err / Rprob_8482)**2 +
                                                  (STANDdic['err_Rnist_8482'] /
                                                  STANDdic['Rnist_8482'])**2)
    Fshelf_8482 = (delW_probe_8482 - (STANDdic['delW_floor_8482'] * (1-BGfrac))) / \
                  (STANDdic['delW_shelf_8482'] - STANDdic['delW_floor_8482'])
    err_Fshelf_8482 = Fshelf_8482*np.sqrt((err_delW_probe_8482/delW_probe_8482)**2 +
                                          (STANDdic['err_delW_floor_8482'] /
                                          STANDdic['delW_floor_8482'])**2 +
                                          (STANDdic['err_delW_shelf_8482'] /
                                          STANDdic['delW_shelf_8482'])**2)
    Ffloor_8482 = 1 - Fshelf_8482 - BGfrac

    Rprob_8682 = enrichDF['EF186']/enrichDF['EF182']
    Rprob_8682 = np.nan_to_num(Rprob_8682)
    Rprob_8682[Rprob_8682 > 2.0] = 2.0
    Rprob_8682[Rprob_8682 < -2.0] = -2.0
    Rprob_8682_err = Rprob_8682 * np.sqrt((enrichDF['EF186err'] / enrichDF['EF186'])**2 +
                                          (enrichDF['EF182err'] / enrichDF['EF182'])**2)
    Rprob_8682_err = np.nan_to_num(Rprob_8682_err)
    Rprob_8682_err[Rprob_8682_err > 1.0] = 1.0
    Rprob_8682_err[Rprob_8682_err < -1.0] = 1.0
    delW_probe_8682 = (Rprob_8682/STANDdic['Rnist_8682']) - 1.
    err_delW_probe_8682 = delW_probe_8682*np.sqrt((Rprob_8682_err / Rprob_8682)**2 +
                                                  (STANDdic['err_Rnist_8682'] /
                                                  STANDdic['Rnist_8682'])**2)
    Fshelf_8682 = (delW_probe_8682 - (STANDdic['delW_floor_8682'] * (1-BGfrac))) / \
                  (STANDdic['delW_shelf_8682'] - STANDdic['delW_floor_8682'])
    Fshelf_8682 = np.nan_to_num(Fshelf_8682)
    err_Fshelf_8682 = Fshelf_8682*np.sqrt((err_delW_probe_8682/delW_probe_8682)**2 +
                                          (STANDdic['err_delW_floor_8682'] /
                                          STANDdic['delW_floor_8682'])**2 +
                                          (STANDdic['err_delW_shelf_8682'] /
                                          STANDdic['delW_shelf_8682'])**2)
    Ffloor_8682 = 1 - Fshelf_8682 - BGfrac

    # Full averaging
    # simm_floor_AVG = (Ffloor_8082 + Ffloor_8382 + Ffloor_8482 + Ffloor_8682) / 4.0
    # err_simm_floor_AVG = np.sqrt(err_Ffloor_8082**2 + err_Ffloor_8382**2 +
    # err_Ffloor_8482**2 + err_Ffloor_8682**2)

    simm_floor_AVG = (Ffloor_8382 + Ffloor_8482 + Ffloor_8682) / 3.0
    simm_shelf_AVG = (Fshelf_8382 + Fshelf_8482 + Fshelf_8682) / 3.0
    err_simm_floor_AVG = err_simm_shelf_AVG = np.sqrt(err_Fshelf_8382**2 + err_Fshelf_8482**2 +
                                                      err_Fshelf_8682**2)

    outDIC = {'f_BG': BGfrac,
              'Fshelf_8082': pd.Series(Fshelf_8082, index=enrichDF.index),
              'Ffloor_8082': pd.Series(Ffloor_8082, index=enrichDF.index),
              'err_Fshelf_8082': pd.Series(err_Fshelf_8082, index=enrichDF.index),
              'err_Ffloor_8082': pd.Series(err_Fshelf_8082, index=enrichDF.index),
              'Fshelf_8382': pd.Series(Fshelf_8382, index=enrichDF.index),
              'Ffloor_8382': pd.Series(Ffloor_8382, index=enrichDF.index),
              'err_Fshelf_8382': pd.Series(err_Fshelf_8382, index=enrichDF.index),
              'err_Ffloor_8382': pd.Series(err_Fshelf_8382, index=enrichDF.index),
              'Fshelf_8482': pd.Series(Fshelf_8482, index=enrichDF.index),
              'Ffloor_8482': pd.Series(Ffloor_8482, index=enrichDF.index),
              'err_Fshelf_8482': pd.Series(err_Fshelf_8482, index=enrichDF.index),
              'err_Ffloor_8482': pd.Series(err_Fshelf_8482, index=enrichDF.index),
              'Fshelf_8682': pd.Series(Fshelf_8682, index=enrichDF.index),
              'Ffloor_8682': pd.Series(Ffloor_8682, index=enrichDF.index),
              'err_Fshelf_8682': pd.Series(err_Fshelf_8682, index=enrichDF.index),
              'err_Ffloor_8682': pd.Series(err_Fshelf_8682, index=enrichDF.index),
              'simm_floor_AVG': pd.Series(simm_floor_AVG, index=enrichDF.index),
              'err_simm_floor_AVG': pd.Series(err_simm_floor_AVG, index=enrichDF.index),
              'simm_shelf_AVG': pd.Series(simm_shelf_AVG, index=enrichDF.index),
              'err_simm_shelf_AVG': pd.Series(err_simm_shelf_AVG, index=enrichDF.index)
              }

    simmsDF = pd.DataFrame(outDIC)
    # simmsDF = simmsDF.fillna(value=0.0)
    simmsDF.DFname = str(enrichDF.DFname)+'_simms'

    if HDFdump:
        HDFfilename = str(simmsDF.DFname)+'.h5'
        hdf = pd.HDFStore(HDFfilename)
        hdf.put(str(simmsDF.DFname), simmsDF, format='table', data_columns=True)
    return simmsDF


def getCPenrichment(fileNAM="AU17_dict.csv", datFOLDpath='/.'):

    # Obtain the data
    useCOLS = [1, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
    enrichDF = pd.read_csv(datFOLDpath+fileNAM, index_col=0, usecols=useCOLS)

    enrichDF.DFname = fileNAM[0:4]
    return enrichDF


def makRBS_DIC(fileNAM, datFOLDpath='/.'):
    rbs_U_COLS = [1, 3, 4, 7, 8, 12, 13]
    rbsDF_U = pd.read_csv(datFOLDpath+fileNAM, index_col=2, usecols=rbs_U_COLS)

    rbs_D_COLS = [0, 2, 5, 6, 9, 10, 11]
    rbsDF_D = pd.read_csv(datFOLDpath+fileNAM, index_col=5, usecols=rbs_D_COLS)

    rbsDIC = {'Udf': rbsDF_U, 'Ddf': rbsDF_D}
    return rbsDIC


def pltCPenrichment(DF):
    nams = ['EF180', 'EF182', 'EF183', 'EF184', 'EF186']
    namsERR = ['EF180err', 'EF182err', 'EF183err', 'EF184err', 'EF186err']

    ax = DF.plot(y=nams[0], yerr=namsERR[0], logy=True, ylim=(0.0001, 1.0), marker='o',
                 linestyle='None', ms=3, mfc='none', alpha=0.5)

    for i, j in enumerate(nams[1::]):
        DF.plot(y=nams[i+1], yerr=namsERR[i+1], ylim=(0.0001, 1.0), marker='o', linestyle='None',
                ms=3, mfc='none', alpha=0.5, ax=ax, logy=True, title='CP enrichments')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    return


def pltArealFRAC(arealfrDF):  # arealfrDF
    nams = ['arealFRAC_floor', 'arealFRAC_shelf']
    namsERR = ['err_arealFRAC_floor', 'err_arealFRAC_shelf']

    ax = arealfrDF.plot(y=nams[0], yerr=namsERR[0], logy=False, marker='o',
                        linestyle='None', ms=3, mfc='none', alpha=0.5)

    for i, j in enumerate(nams[1::]):
        arealfrDF.plot(y=nams[i+1], yerr=namsERR[i+1], marker='o', linestyle='None',
                       ms=3, mfc='none', alpha=0.5, ax=ax, logy=False,
                       title='Fractional Areal Density')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    return


def makSTANDARDS():
    # Standards go like W-180, W-182, W-183, W-184, W-186
    STAND_NIST = [0.0012, 0.2650, 0.1431, 0.3064, 0.2843]
    err_STAND_NIST = [0.0002, 0.0032, 0.0008, 0.0004, 0.0038]

    # Old estimates should update with new LAMS -- EAU 1/10/2018
    # STAND_ULTRAMET = [0.0180,0.2652,0.1402,0.3055,0.2883]
    # err_STAND_ULTRAMET = [0.1,0.0065,0.0016,0.0061,0.0033]
    # try #2 from LAMS
    STAND_ULTRAMET = [0.0012, 0.2652, 0.1397, 0.3056, 0.2860]
    err_STAND_ULTRAMET = [0.0100, 0.0065, 0.0016, 0.0061, 0.0033]

    # ICPMS-TOF-MS
    # STAND_ORNL182 = [0.0005,0.9299,0.0284,0.0279,0.0138]
    # err_STAND_ORNL182 = [0.0001,0.0036,0.0031,0.0024,0.0006]

    # ICPMS-LA-MS
    STAND_ORNL182 = [6.2169e-5, 0.9307, 0.0284, 0.0279, 0.0140]
    err_STAND_ORNL182 = [2e-5, 0.005, 0.005, 0.005, 0.005]

    R_NIST_8082 = STAND_NIST[0]/STAND_NIST[1]
    R_NIST_8382 = STAND_NIST[2]/STAND_NIST[1]
    R_NIST_8482 = STAND_NIST[3]/STAND_NIST[1]
    R_NIST_8682 = STAND_NIST[4]/STAND_NIST[1]
    err_NIST_8082 = STAND_NIST[0]*np.sqrt((err_STAND_NIST[0]/STAND_NIST[0])**2 +
                                          (err_STAND_NIST[1]/STAND_NIST[1])**2)
    err_NIST_8382 = STAND_NIST[2]*np.sqrt((err_STAND_NIST[2]/STAND_NIST[2])**2 +
                                          (err_STAND_NIST[1]/STAND_NIST[1])**2)
    err_NIST_8482 = STAND_NIST[3]*np.sqrt((err_STAND_NIST[3]/STAND_NIST[3])**2 +
                                          (err_STAND_NIST[1]/STAND_NIST[1])**2)
    err_NIST_8682 = STAND_NIST[4]*np.sqrt((err_STAND_NIST[4]/STAND_NIST[4])**2 +
                                          (err_STAND_NIST[1]/STAND_NIST[1])**2)

    R_ULTRAMET_8082 = STAND_ULTRAMET[0]/STAND_ULTRAMET[1]
    R_ULTRAMET_8382 = STAND_ULTRAMET[2]/STAND_ULTRAMET[1]
    R_ULTRAMET_8482 = STAND_ULTRAMET[3]/STAND_ULTRAMET[1]
    R_ULTRAMET_8682 = STAND_ULTRAMET[4]/STAND_ULTRAMET[1]
    err_ULTRAMET_8082 = STAND_ULTRAMET[0]*np.sqrt((err_STAND_ULTRAMET[0]/STAND_ULTRAMET[0])**2 +
                                                  (err_STAND_ULTRAMET[1]/STAND_ULTRAMET[1])**2)
    err_ULTRAMET_8382 = STAND_ULTRAMET[2]*np.sqrt((err_STAND_ULTRAMET[2]/STAND_ULTRAMET[2])**2 +
                                                  (err_STAND_ULTRAMET[1]/STAND_ULTRAMET[1])**2)
    err_ULTRAMET_8482 = STAND_ULTRAMET[3]*np.sqrt((err_STAND_ULTRAMET[3]/STAND_ULTRAMET[3])**2 +
                                                  (err_STAND_ULTRAMET[1]/STAND_ULTRAMET[1])**2)
    err_ULTRAMET_8682 = STAND_ULTRAMET[4]*np.sqrt((err_STAND_ULTRAMET[4]/STAND_ULTRAMET[4])**2 +
                                                  (err_STAND_ULTRAMET[1]/STAND_ULTRAMET[1])**2)

    R_ORNL_8082 = STAND_ORNL182[0]/STAND_ORNL182[1]
    R_ORNL_8382 = STAND_ORNL182[2]/STAND_ORNL182[1]
    R_ORNL_8482 = STAND_ORNL182[3]/STAND_ORNL182[1]
    R_ORNL_8682 = STAND_ORNL182[4]/STAND_ORNL182[1]
    err_ORNL_8082 = STAND_ORNL182[0]*np.sqrt((err_STAND_ORNL182[0]/STAND_ORNL182[0])**2 +
                                             (err_STAND_ORNL182[1]/STAND_ORNL182[1])**2)
    err_ORNL_8382 = STAND_ORNL182[2]*np.sqrt((err_STAND_ORNL182[2]/STAND_ORNL182[2])**2 +
                                             (err_STAND_ORNL182[1]/STAND_ORNL182[1])**2)
    err_ORNL_8482 = STAND_ORNL182[3]*np.sqrt((err_STAND_ORNL182[3]/STAND_ORNL182[3])**2 +
                                             (err_STAND_ORNL182[1]/STAND_ORNL182[1])**2)
    err_ORNL_8682 = STAND_ORNL182[4]*np.sqrt((err_STAND_ORNL182[4]/STAND_ORNL182[4])**2 +
                                             (err_STAND_ORNL182[1]/STAND_ORNL182[1])**2)

    delW_floor_8082 = (R_ULTRAMET_8082/R_NIST_8082) - 1.
    delW_floor_8382 = (R_ULTRAMET_8382/R_NIST_8382) - 1.
    delW_floor_8482 = (R_ULTRAMET_8482/R_NIST_8482) - 1.
    delW_floor_8682 = (R_ULTRAMET_8682/R_NIST_8682) - 1.
    err_delW_floor_8082 = delW_floor_8082*np.sqrt((err_ULTRAMET_8082/R_ULTRAMET_8082)**2 +
                                                  (err_NIST_8082/R_NIST_8082)**2)
    err_delW_floor_8382 = delW_floor_8382*np.sqrt((err_ULTRAMET_8382/R_ULTRAMET_8382)**2 +
                                                  (err_NIST_8382/R_NIST_8382)**2)
    err_delW_floor_8482 = delW_floor_8482*np.sqrt((err_ULTRAMET_8482/R_ULTRAMET_8482)**2 +
                                                  (err_NIST_8482/R_NIST_8482)**2)
    err_delW_floor_8682 = delW_floor_8682*np.sqrt((err_ULTRAMET_8682/R_ULTRAMET_8682)**2 +
                                                  (err_NIST_8682/R_NIST_8682)**2)

    delW_shelf_8082 = (R_ORNL_8082/R_NIST_8082) - 1.
    delW_shelf_8382 = (R_ORNL_8382/R_NIST_8382) - 1.
    delW_shelf_8482 = (R_ORNL_8482/R_NIST_8482) - 1.
    delW_shelf_8682 = (R_ORNL_8682/R_NIST_8682) - 1.
    err_delW_shelf_8082 = delW_floor_8082*np.sqrt((err_ORNL_8082/R_ORNL_8082)**2 +
                                                  (err_NIST_8082/R_NIST_8082)**2)
    err_delW_shelf_8382 = delW_floor_8382*np.sqrt((err_ORNL_8382/R_ORNL_8382)**2 +
                                                  (err_NIST_8382/R_NIST_8382)**2)
    err_delW_shelf_8482 = delW_floor_8482*np.sqrt((err_ORNL_8482/R_ORNL_8482)**2 +
                                                  (err_NIST_8482/R_NIST_8482)**2)
    err_delW_shelf_8682 = delW_floor_8682*np.sqrt((err_ORNL_8682/R_ORNL_8682)**2 +
                                                  (err_NIST_8682/R_NIST_8682)**2)

    outDIC = {'delW_floor_8082': delW_floor_8082,
              'delW_floor_8382': delW_floor_8382,
              'delW_floor_8482': delW_floor_8482,
              'delW_floor_8682': delW_floor_8682,
              'err_delW_floor_8082': np.abs(err_delW_floor_8082),
              'err_delW_floor_8382': np.abs(err_delW_floor_8382),
              'err_delW_floor_8482': np.abs(err_delW_floor_8482),
              'err_delW_floor_8682': np.abs(err_delW_floor_8682),
              'delW_shelf_8082': delW_shelf_8082,
              'delW_shelf_8382': delW_shelf_8382,
              'delW_shelf_8482': delW_shelf_8482,
              'delW_shelf_8682': delW_shelf_8682,
              'err_delW_shelf_8082': np.abs(err_delW_shelf_8082),
              'err_delW_shelf_8382': np.abs(err_delW_shelf_8382),
              'err_delW_shelf_8482': np.abs(err_delW_shelf_8482),
              'err_delW_shelf_8682': np.abs(err_delW_shelf_8682),
              'Rnist_8082': R_NIST_8082,
              'Rnist_8382': R_NIST_8382,
              'Rnist_8482': R_NIST_8482,
              'Rnist_8682': R_NIST_8682,
              'err_Rnist_8082': err_NIST_8082,
              'err_Rnist_8382': err_NIST_8382,
              'err_Rnist_8482': err_NIST_8482,
              'err_Rnist_8682': err_NIST_8682}
    return outDIC
