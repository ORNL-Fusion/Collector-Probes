import numpy as np
import pandas as pd


def makSIMMS(enrichDF, BGfrac=0.0435, HDFdump=0):
    # First get the source standard numbers
    STANDdic = makSTANDARDS()

    # Now what through the analysis for each ratio
    Rprob_8082 = pd.Series(data=enrichDF['W180R']/enrichDF['W182R'],
                           index=enrichDF.index)
    Rprob_8082_err = Rprob_8082 * np.sqrt((enrichDF['W180 error'] / enrichDF['W180R'])**2 +
                                          (enrichDF['W182 error'] / enrichDF['W182R'])**2)
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

    Rprob_8382 = pd.Series(data=enrichDF['W183R']/enrichDF['W182R'],
                           index=enrichDF.index)
    Rprob_8382_err = Rprob_8382 * np.sqrt((enrichDF['W183 error'] / enrichDF['W183R'])**2 +
                                          (enrichDF['W182 error'] / enrichDF['W182R'])**2)
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

    Rprob_8482 = enrichDF['W184R']/enrichDF['W182R']
    Rprob_8482_err = Rprob_8482 * np.sqrt((enrichDF['W184 error'] / enrichDF['W184R'])**2 +
                                          (enrichDF['W182 error'] / enrichDF['W182R'])**2)
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

    Rprob_8682 = enrichDF['W186R']/enrichDF['W182R']
    Rprob_8682_err = Rprob_8682 * np.sqrt((enrichDF['W186 error'] / enrichDF['W186R'])**2 +
                                          (enrichDF['W182 error'] / enrichDF['W182R'])**2)
    delW_probe_8682 = (Rprob_8682/STANDdic['Rnist_8682']) - 1.
    err_delW_probe_8682 = delW_probe_8682*np.sqrt((Rprob_8682_err / Rprob_8682)**2 +
                                                  (STANDdic['err_Rnist_8682'] /
                                                  STANDdic['Rnist_8682'])**2)
    Fshelf_8682 = (delW_probe_8682 - (STANDdic['delW_floor_8682'] * (1-BGfrac))) / \
                  (STANDdic['delW_shelf_8682'] - STANDdic['delW_floor_8682'])
    err_Fshelf_8682 = Fshelf_8682*np.sqrt((err_delW_probe_8682/delW_probe_8682)**2 +
                                          (STANDdic['err_delW_floor_8682'] /
                                          STANDdic['delW_floor_8682'])**2 +
                                          (STANDdic['err_delW_shelf_8682'] /
                                          STANDdic['delW_shelf_8682'])**2)
    Ffloor_8682 = 1 - Fshelf_8682 - BGfrac

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
              'err_Ffloor_8682': pd.Series(err_Fshelf_8682, index=enrichDF.index)
              }

    simmsDF = pd.DataFrame(outDIC)
    simmsDF.DFname = str(enrichDF.DFname)+'_simms'

    if HDFdump:
        HDFfilename = str(simmsDF.DFname)+'.h5'
        hdf = pd.HDFStore(HDFfilename)
        hdf.put(str(simmsDF.DFname), simmsDF, format='table', data_columns=True)
    return simmsDF


def getCPenrichment(fileNAM="AU17.xlsx", shNAM='TriPlot', rowSTART=2,
                    enrichCOLS=[1, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18],
                    rbsCOLS=[20, 21, 22, 24, 25, 26]):

    # Obtain the data
    enrichDF = pd.read_excel(fileNAM, skiprows=rowSTART, sheetname=shNAM, usecols=enrichCOLS,
                             index_col=[0], engine='xlrd')
    rbsDF = pd.read_excel(fileNAM, skiprows=rowSTART, sheetname=shNAM, usecols=rbsCOLS,
                          skip_footer=360, engine='xlrd')
    enrichDF.DFname = fileNAM[0:4]
    return rbsDF, enrichDF


def pltCPenrichment(DF):
    nams = ['W180R', 'W182R', 'W183R', 'W184R', 'W186R']
    namsERR = ['W180 error', 'W182 error', 'W183 error', 'W184 error', 'W186 error']

    ax = DF.plot(y=nams[0], yerr=namsERR[0], logy=True, ylim=(0.0001, 1.0), marker='o',
                 linestyle='None', ms=3, mfc='none', alpha=0.5)

    for i, j in enumerate(nams[1::]):
        DF.plot(y=nams[i+1], yerr=namsERR[i+1], ylim=(0.0001, 1.0), marker='o', linestyle='None',
                ms=3, mfc='none', alpha=0.5, ax=ax, logy=True, title='CP enrichments')
    ax.legend()
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
