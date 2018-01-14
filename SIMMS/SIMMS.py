

import numpy as np

def makSTANDARDS():
    # Standards go like W-180, W-182, W-183, W-184, W-186
    STAND_NIST = [0.0012,0.2650,0.1431,0.3064,0.2843]
    err_STAND_NIST = [0.0002,0.0032,0.0008,0.0004,0.0038]

    # Old estimates should update with new LAMS -- EAU 1/10/2018
    # STAND_ULTRAMET = [0.0180,0.2652,0.1402,0.3055,0.2883]
    # err_STAND_ULTRAMET = [0.1,0.0065,0.0016,0.0061,0.0033]
    # try #2 from LAMS
    STAND_ULTRAMET = [0.0012,0.2652,0.1397,0.3056,0.2860]
    err_STAND_ULTRAMET = [0.0100,0.0065,0.0016,0.0061,0.0033]

    # ICPMS-TOF-MS
    # STAND_ORNL182 = [0.0005,0.9299,0.0284,0.0279,0.0138]
    # err_STAND_ORNL182 = [0.0001,0.0036,0.0031,0.0024,0.0006]

    # ICPMS-LA-MS
    STAND_ORNL182 = [6.2169e-5,0.9307,0.0284,0.0279,0.0140]
    err_STAND_ORNL182 = [2e-5,0.005,0.005,0.005,0.005]

    R_NIST_8082 = STAND_NIST[0]/STAND_NIST[1]
    R_NIST_8382 = STAND_NIST[2]/STAND_NIST[1]
    R_NIST_8482 = STAND_NIST[3]/STAND_NIST[1]
    R_NIST_8682 = STAND_NIST[4]/STAND_NIST[1]
    err_NIST_8082 = STAND_NIST[0]*np.sqrt((err_STAND_NIST[0]/STAND_NIST[0])**2+
                                           (err_STAND_NIST[1]/STAND_NIST[1])**2)
    err_NIST_8382 = STAND_NIST[2]*np.sqrt((err_STAND_NIST[2]/STAND_NIST[2])**2+
                                           (err_STAND_NIST[1]/STAND_NIST[1])**2)
    err_NIST_8482 = STAND_NIST[3]*np.sqrt((err_STAND_NIST[3]/STAND_NIST[3])**2+
                                           (err_STAND_NIST[1]/STAND_NIST[1])**2)
    err_NIST_8682 = STAND_NIST[4]*np.sqrt((err_STAND_NIST[4]/STAND_NIST[4])**2+
                                           (err_STAND_NIST[1]/STAND_NIST[1])**2)

    R_ULTRAMET_8082 = STAND_ULTRAMET[0]/STAND_ULTRAMET[1]
    R_ULTRAMET_8382 = STAND_ULTRAMET[2]/STAND_ULTRAMET[1]
    R_ULTRAMET_8482 = STAND_ULTRAMET[3]/STAND_ULTRAMET[1]
    R_ULTRAMET_8682 = STAND_ULTRAMET[4]/STAND_ULTRAMET[1]
    err_ULTRAMET_8082 = STAND_ULTRAMET[0]*np.sqrt((err_STAND_ULTRAMET[0]/STAND_ULTRAMET[0])**2+
                                           (err_STAND_ULTRAMET[1]/STAND_ULTRAMET[1])**2)
    err_ULTRAMET_8382 = STAND_ULTRAMET[2]*np.sqrt((err_STAND_ULTRAMET[2]/STAND_ULTRAMET[2])**2+
                                           (err_STAND_ULTRAMET[1]/STAND_ULTRAMET[1])**2)
    err_ULTRAMET_8482 = STAND_ULTRAMET[3]*np.sqrt((err_STAND_ULTRAMET[3]/STAND_ULTRAMET[3])**2+
                                           (err_STAND_ULTRAMET[1]/STAND_ULTRAMET[1])**2)
    err_ULTRAMET_8682 = STAND_ULTRAMET[4]*np.sqrt((err_STAND_ULTRAMET[4]/STAND_ULTRAMET[4])**2+
                                           (err_STAND_ULTRAMET[1]/STAND_ULTRAMET[1])**2)

    R_ORNL_8082 = STAND_ORNL182[0]/STAND_ORNL182[1]
    R_ORNL_8382 = STAND_ORNL182[2]/STAND_ORNL182[1]
    R_ORNL_8482 = STAND_ORNL182[3]/STAND_ORNL182[1]
    R_ORNL_8682 = STAND_ORNL182[4]/STAND_ORNL182[1]
    err_ORNL_8082 = STAND_ORNL182[0]*np.sqrt((err_STAND_ORNL182[0]/STAND_ORNL182[0])**2+
                                           (err_STAND_ORNL182[1]/STAND_ORNL182[1])**2)
    err_ORNL_8382 = STAND_ORNL182[2]*np.sqrt((err_STAND_ORNL182[2]/STAND_ORNL182[2])**2+
                                           (err_STAND_ORNL182[1]/STAND_ORNL182[1])**2)
    err_ORNL_8482 = STAND_ORNL182[3]*np.sqrt((err_STAND_ORNL182[3]/STAND_ORNL182[3])**2+
                                           (err_STAND_ORNL182[1]/STAND_ORNL182[1])**2)
    err_ORNL_8682 = STAND_ORNL182[4]*np.sqrt((err_STAND_ORNL182[4]/STAND_ORNL182[4])**2+
                                           (err_STAND_ORNL182[1]/STAND_ORNL182[1])**2)

    delW_floor_8082 = (R_ULTRAMET_8082/R_NIST_8082) - 1.
    delW_floor_8382 = (R_ULTRAMET_8382/R_NIST_8382) - 1.
    delW_floor_8482 = (R_ULTRAMET_8482/R_NIST_8482) - 1.
    delW_floor_8682 = (R_ULTRAMET_8682/R_NIST_8682) - 1.
    err_delW_floor_8082 = delW_floor_8082*np.sqrt((err_ULTRAMET_8082/R_ULTRAMET_8082)**2+
                                           (err_NIST_8082/R_NIST_8082)**2)
    err_delW_floor_8382 = delW_floor_8382*np.sqrt((err_ULTRAMET_8382/R_ULTRAMET_8382)**2+
                                           (err_NIST_8382/R_NIST_8382)**2)
    err_delW_floor_8482 = delW_floor_8482*np.sqrt((err_ULTRAMET_8482/R_ULTRAMET_8482)**2+
                                           (err_NIST_8482/R_NIST_8482)**2)
    err_delW_floor_8682 = delW_floor_8682*np.sqrt((err_ULTRAMET_8682/R_ULTRAMET_8682)**2+
                                           (err_NIST_8682/R_NIST_8682)**2)

    delW_shelf_8082 = (R_ORNL_8082/R_NIST_8082) - 1.
    delW_shelf_8382 = (R_ORNL_8382/R_NIST_8382) - 1.
    delW_shelf_8482 = (R_ORNL_8482/R_NIST_8482) - 1.
    delW_shelf_8682 = (R_ORNL_8682/R_NIST_8682) - 1.
    err_delW_shelf_8082 = delW_floor_8082*np.sqrt((err_ORNL_8082/R_ORNL_8082)**2+
                                                  (err_NIST_8082/R_NIST_8082)**2)
    err_delW_shelf_8382 = delW_floor_8382*np.sqrt((err_ORNL_8382/R_ORNL_8382)**2+
                                                  (err_NIST_8382/R_NIST_8382)**2)
    err_delW_shelf_8482 = delW_floor_8482*np.sqrt((err_ORNL_8482/R_ORNL_8482)**2+
                                                  (err_NIST_8482/R_NIST_8482)**2)
    err_delW_shelf_8682 = delW_floor_8682*np.sqrt((err_ORNL_8682/R_ORNL_8682)**2+
                                                  (err_NIST_8682/R_NIST_8682)**2)

    outDIC = {'delW_floor_8082':delW_floor_8082,
              'delW_floor_8382':delW_floor_8382,
              'delW_floor_8482':delW_floor_8482,
              'delW_floor_8682':delW_floor_8682,
              'err_delW_floor_8082':np.abs(err_delW_floor_8082),
              'err_delW_floor_8382':np.abs(err_delW_floor_8382),
              'err_delW_floor_8482':np.abs(err_delW_floor_8482),
              'err_delW_floor_8682':np.abs(err_delW_floor_8682),
              'delW_shelf_8082':delW_shelf_8082,
              'delW_shelf_8382':delW_shelf_8382,
              'delW_shelf_8482':delW_shelf_8482,
              'delW_shelf_8682':delW_shelf_8682,
              'err_delW_shelf_8082':np.abs(err_delW_shelf_8082),
              'err_delW_shelf_8382':np.abs(err_delW_shelf_8382),
              'err_delW_shelf_8482':np.abs(err_delW_shelf_8482),
              'err_delW_shelf_8682':np.abs(err_delW_shelf_8682),
              'R_NIST_8082':R_NIST_8082,
              'R_NIST_8382':R_NIST_8382,
              'R_NIST_8482':R_NIST_8482,
              'R_NIST_8682':R_NIST_8682
              }

    return outDIC


def makSIMMS(xLOCs, arr180='none',arr182='none',arr183='none',
             arr184='none', arr186='none', BGfrac=0.0435):
    # First get the source standard numbers
    STANDdic = makSTANDARDS()

    # Now what through the analysis for each ratio
    if arr180 is not 'none':
        R_prob_8082 = arr180/arr182
        delW_probe_8082 = (R_prob_8082/STANDdic['R_NIST_8082']) - 1.
        f_shelf_8082 = (delW_probe_8082 - (STANDdic['delW_floor_8082']*(1-BGfrac)))/(STANDdic['delW_shelf_8082']-STANDdic['delW_floor_8082'])
        f_floor_8082 = 1 - f_shelf_8082 - BGfrac
    else:
        print("no W-180 data")
    if arr183 is not 'none':
        R_prob_8382 = arr183/arr182
        delW_probe_8382 = (R_prob_8382/STANDdic['R_NIST_8382']) - 1.
        f_shelf_8382 = (delW_probe_8382 - (STANDdic['delW_floor_8382']*(1-BGfrac)))/(STANDdic['delW_shelf_8382']-STANDdic['delW_floor_8382'])
        f_floor_8382 = 1 - f_shelf_8382 - BGfrac
    else:
        print("no W-183 data")
    if arr184 is not 'none':
        R_prob_8482 = arr184/arr182
        delW_probe_8482 = (R_prob_8482/STANDdic['R_NIST_8482']) - 1.
        f_shelf_8482 = (delW_probe_8482 - (STANDdic['delW_floor_8482']*(1-BGfrac)))/(STANDdic['delW_shelf_8482']-STANDdic['delW_floor_8482'])
        f_floor_8482 = 1 - f_shelf_8482 - BGfrac
    if arr186 is not 'none':
        R_prob_8682 = arr186/arr182
        delW_probe_8682 = (R_prob_8682/STANDdic['R_NIST_8682']) - 1.
        f_shelf_8682 = (delW_probe_8682 - (STANDdic['delW_floor_8682']*(1-BGfrac)))/(STANDdic['delW_shelf_8682']-STANDdic['delW_floor_8682'])
        f_floor_8682 = 1 - f_shelf_8682 - BGfrac
    else:
        print("no W-186 data")

    if arr180 is not 'none' and arr183 is 'none' and arr184 is 'none' and arr186 is 'none':
        outDIC = {'f_BG':BGfrac,
                  'f_shelf_8082':f_shelf_8082,
                  'f_floor_8082':f_floor_8082,
                  'f_xlocs':xLOCs}
    if arr180 is 'none' and arr183 is not 'none' and arr184 is 'none' and arr186 is 'none':
        outDIC = {'f_BG':BGfrac,
                  'f_shelf_8382':f_shelf_8382,
                  'f_floor_8382':f_floor_8382,
                  'f_xlocs':xLOCs}
    if arr180 is not 'none' and arr183 is not 'none' and arr184 is 'none' and arr186 is 'none':
        outDIC = {'f_BG':BGfrac,
                  'f_shelf_8082':f_shelf_8082,
                  'f_floor_8082':f_floor_8082,
                  'f_shelf_8382':f_shelf_8382,
                  'f_floor_8382':f_floor_8382,
                  'f_xlocs':xLOCs}
    if arr180 is not 'none' and arr183 is not 'none' and arr184 is not 'none' and arr186 is 'none':
        outDIC = {'f_BG':BGfrac,
                  'f_shelf_8082':f_shelf_8082,
                  'f_floor_8082':f_floor_8082,
                  'f_shelf_8382':f_shelf_8382,
                  'f_floor_8382':f_floor_8382,
                  'f_shelf_8482':f_shelf_8482,
                  'f_floor_8482':f_floor_8482,
                  'f_xlocs':xLOCs}
    if arr180 is not 'none' and arr183 is not 'none' and arr184 is not 'none' and arr186 is not 'none':
        outDIC = {'f_BG':BGfrac,
                  'f_shelf_8082':f_shelf_8082,
                  'f_floor_8082':f_floor_8082,
                  'f_shelf_8382':f_shelf_8382,
                  'f_floor_8382':f_floor_8382,
                  'f_shelf_8482':f_shelf_8482,
                  'f_floor_8482':f_floor_8482,
                   'f_shelf_8682':f_shelf_8682,
                   'f_floor_8682':f_floor_8682,
                  'f_xlocs':xLOCs}

    return outDIC


def getCPenrichment(fileNAM="AU17.xlsx"):
    from xlrd import open_workbook
    #Open the sheet
    book = open_workbook(fileNAM)
    sheet = book.sheet_by_index(3) #If your data is on sheet 4

    #Define the lists
    scanTime = []
    scanDist = []
    W180R = []
    W180_err = []
    W182R = []
    W182_err = []
    W183R = []
    W183_err = []
    W184R = []
    W184_err = []
    W186R = []
    W186_err = []
    W_RBS_x = []
    W_RBS = []
    W_RBS_err = []
    RminRsep_x = []

    #Obtain the data
    for row in range(4, 400): #start from 1, to leave out row 0
        scanTime.append(float(sheet.cell(row, 0).value)) #extract from first col
        scanDist.append(float(sheet.cell(row, 1).value))
        W180R.append(float(sheet.cell(row, 9).value))
        W180_err.append(float(sheet.cell(row, 10).value))
        W182R.append(float(sheet.cell(row, 11).value))
        W182_err.append(float(sheet.cell(row, 12).value))
        W183R.append(float(sheet.cell(row, 13).value))
        W183_err.append(float(sheet.cell(row, 14).value))
        W184R.append(float(sheet.cell(row, 15).value))
        W184_err.append(float(sheet.cell(row, 16).value))
        W186R.append(float(sheet.cell(row, 17).value))
        W186_err.append(float(sheet.cell(row, 18).value))


    for row in range(3,23):
        W_RBS_x.append(float(sheet.cell(row, 20).value))
        W_RBS.append(float(sheet.cell(row, 21).value))
        W_RBS_err.append(float(sheet.cell(row, 22).value))
        RminRsep_x.append(float(sheet.cell(row,26).value))

    #Generate error band values
    RBS_perr = [i - j for i, j in zip(W_RBS, W_RBS_err)]
    RBS_merr = [i + j for i, j in zip(W_RBS, W_RBS_err)]
    W180merr = [i - j for i, j in zip(W180R, W180_err)]
    W180perr = [i + j for i, j in zip(W180R, W180_err)]
    W182merr = [i - j for i, j in zip(W182R, W182_err)]
    W182perr = [i + j for i, j in zip(W182R, W182_err)]
    W183merr = [i - j for i, j in zip(W183R, W183_err)]
    W183perr = [i + j for i, j in zip(W183R, W183_err)]
    W184merr = [i - j for i, j in zip(W184R, W184_err)]
    W184perr = [i + j for i, j in zip(W184R, W184_err)]
    W186merr = [i - j for i, j in zip(W186R, W186_err)]
    W186perr = [i + j for i, j in zip(W186R, W186_err)]

    outDIC ={'Xprobe':np.array(scanDist),
             'W180ric':np.array(W180R),
             'W180_err':np.array(W180_err),
             'W182ric':np.array(W182R),
             'W182_err':np.array(W182_err),
             'W183ric':np.array(W183R),
             'W183_err':np.array(W183_err),
             'W184ric':np.array(W183R),
             'W184_err':np.array(W183_err),
             'W186ric':np.array(W183R),
             'W186_err':np.array(W183_err)
             }

    return outDIC
