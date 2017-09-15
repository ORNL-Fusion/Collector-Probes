# In separate terminal use
# ssh -Y -p 2039 -L 8000:atlas.gat.com:8000 username@cybele.gat.com
# where you substitute in your username. The MDSplus connection is
# then MDSplus.Connection('localhost').
from __future__ import print_function
import numpy as np
import scipy.interpolate as scinter

# Normal load_gfile_d3d except needed to cast 'g['lcfs']' as an int.
import MDSplus as mds
import EFIT.load_gfile_d3d as loadg
import meas_locations as geo


# Return dictionary of average R - Rsep value for each of the probes.
#
# shots     = list of shots in for. If only one shot still enter as list.
# r_probe   = radial position of probe holder tip (cm).
# location  = location along probe where you want the R - Rsep value (cm)
# startTime = start of time range to be averaged over (ms).
# endTime   = end of time range (ms).
def avg_Rsep(shots, r_probe, location, writeToFile=False, filename="IDidNotEnterAFilename.txt",
             server='atlas.gat.com', startTime=2500, endTime=5000):
    time = startTime
    count = 0

    # Dict to hold each r-rsep value for calculating statistical values later.
    rminrsep_values = {'ad': [], 'au': [], 'bd': [], 'bu': [], 'cd': [], 'cu': []}

    # sums = {'ad':0, 'au':0, 'bd':0, 'bu':0, 'cd':0, 'cu':0}

    # Radial position of each probe.
    rad_pos = {
        'ad': geo.calc_R_measAD(r_probe, location),
        'au': geo.calc_R_measAU(r_probe, location),
        'bd': geo.calc_R_measBD(r_probe, location),
        'bu': geo.calc_R_measBU(r_probe, location),
        'cd': geo.calc_R_measCD(r_probe, location),
        'cu': geo.calc_R_measCU(r_probe, location)}

    # Make an MDSplus connection only once.
    # This makes the EFIT data grabs much faster.
    MDSplusCONN = mds.Connection(server)

    for shot in shots:
        time = startTime
        while time <= endTime:

            print("Shot:     " + str(shot))
            print("Location: " + str(location))
            print("Time:     " + str(time))

            # Lines from Zeke's code.
            parmDICT = loadg.read_g_file_mds(shot, time,
                                             connection=MDSplusCONN,
                                             write2file=False)
            Rs, Zs = np.meshgrid(parmDICT['R'], parmDICT['Z'])
            Zes = np.copy(parmDICT['lcfs'][:, 1][13:-12])
            Res = np.copy(parmDICT['lcfs'][:, 0][13:-12])
            f_Rs = scinter.interp1d(Zes, Res, assume_sorted=False)

            # R_Sep for each z location of the three probes.
            rSep = {}
            rSep['a'] = f_Rs(-0.18)
            rSep['b'] = f_Rs(-0.1546)
            rSep['c'] = f_Rs(-0.2054)

            # rSep is in meters. Convert to cm.
            for value in rSep:
                rSep[value] = rSep[value] * 100.0

            # Put into dictionary of lists containing rminrsep values for
            # each probe.
            rminrsep_values['ad'].append(rad_pos['ad'] - rSep['a'])
            rminrsep_values['au'].append(rad_pos['au'] - rSep['a'])
            rminrsep_values['bd'].append(rad_pos['bd'] - rSep['b'])
            rminrsep_values['bu'].append(rad_pos['bu'] - rSep['b'])
            rminrsep_values['cd'].append(rad_pos['cd'] - rSep['c'])
            rminrsep_values['cu'].append(rad_pos['cu'] - rSep['c'])

            # Sums for calculating avg R - Rsep.
#            sums['ad'] += rad_pos['ad'] - rSep['a']
#            sums['au'] += rad_pos['au'] - rSep['a']
#            sums['bd'] += rad_pos['bd'] - rSep['b']
#            sums['bu'] += rad_pos['bu'] - rSep['b']
#            sums['cd'] += rad_pos['cd'] - rSep['c']
#            sums['cu'] += rad_pos['cu'] - rSep['c']

            # Next time step.
            time += 500
            count += 1
            print("\n")

    # Put rminrsep values into np.array to calculate mean and std. dev.
    avg_rminrsep = {'ad': 0, 'ad_err': 0, 'au': 0, 'au_err': 0, 'bd': 0, 'bd_err': 0, 'bu': 0,
                    'bu_err': 0, 'cd': 0, 'cd_err': 0, 'cu': 0, 'cu_err': 0}
    for key in rminrsep_values:
        rminrsep_values[key] = np.array(rminrsep_values[key])

        # Can use same key for both since they have the same names.
        avg_rminrsep[key] = np.average(rminrsep_values[key])

    # Calc std. dev. for each list of rminrsep values.
    avg_rminrsep['ad_err'] = np.std(rminrsep_values['ad'])
    avg_rminrsep['au_err'] = np.std(rminrsep_values['au'])
    avg_rminrsep['bd_err'] = np.std(rminrsep_values['bd'])
    avg_rminrsep['bu_err'] = np.std(rminrsep_values['bu'])
    avg_rminrsep['cd_err'] = np.std(rminrsep_values['cd'])
    avg_rminrsep['cu_err'] = np.std(rminrsep_values['cu'])

    # Calculate avg R - Rsep
#    avg_RminRseps = {
#        'location':location,
#        'ad':float(sums['ad']) / count,
#        'au':float(sums['au']) / count,
#        'bd':float(sums['bd']) / count,
#        'bu':float(sums['bu']) / count,
#        'cd':float(sums['cd']) / count,
#        'cu':float(sums['cu']) / count}

    # print avg_RminRseps['ad']

    # Save to a txt file.
    if writeToFile:
        # filename = raw_input('Enter filename: ')
        f = open(filename, 'a')
        f.write("Shot: " + str(shot) + "\n")
        for key, value in avg_rminrsep.iteritems():
            f.write(key + " ")
            f.write(str(value))
            f.write("\n")
        f.write("\n")
        f.close()

    return avg_rminrsep


def avg_Rsep_all(shots, r_probe, locations, writeToFile=False,
                 filename="IDidNotEnterAFilename.txt", Etree='EFIT01', server='atlas.gat.com',
                 startTime=2500, endTime=5000, step=500):
    time = startTime
    count = 0

    # main_dict = {} ### not used? Is there a reason to keep?
    tmp_dict = {}

    # Dict to hold each r-rsep value for calculating statistical values later.
    for shot in shots:
        tmp_dict[str(shot)] = {}
        time = startTime
        while time <= endTime:
            tmp_dict[str(shot)][str(time)] = {}
            for loc in locations:
                tmp_dict[str(shot)][str(time)][str(loc)] = {}
            time += step

    # Make an MDSplus connection only once.
    # This makes the EFIT data grabs much faster.
    MDSplusCONN = mds.Connection(server)

    for shot in shots:
        time = startTime
        while time <= endTime:

            print("Shot:     " + str(shot), end='\n')
            # print("Location: " + str(location))
            print("Time:     " + str(time), end='\n')

            # Lines from Zeke's code.
            parmDICT = loadg.read_g_file_mds(shot, time,
                                             connection=MDSplusCONN,
                                             write2file=False,
                                             tree=Etree)
            Rs, Zs = np.meshgrid(parmDICT['R'], parmDICT['Z'])
            Z_axis = parmDICT['ZmAxis']
            R_axis = parmDICT['RmAxis']
            Zes = np.copy(parmDICT['lcfs'][:, 1][13:-12])
            Res = np.copy(parmDICT['lcfs'][:, 0][13:-12])
            f_Rs = scinter.interp1d(Zes, Res, assume_sorted=False)

            # translate up to outbard midplane
            Rs_trunc = Rs > R_axis
            f_psiN = scinter.Rbf(Rs[Rs_trunc],
                                 Zs[Rs_trunc],
                                 parmDICT['psiRZn'][Rs_trunc], function='linear')
            f_Romp = scinter.interp2d(parmDICT['psiRZn'][Rs_trunc], Zs[Rs_trunc], Rs[Rs_trunc])

            # R_Sep for each z location of the three probes.
            rSep = {}
            rSep['a'] = f_Rs(-0.18)
            rSep['b'] = f_Rs(-0.1546)
            rSep['c'] = f_Rs(-0.2054)

            # The R value of the separatrix at the omp. Convert to cm.
            rSep_omp = f_Rs(Z_axis) * 100.0

            # rSep is in meters. Convert to cm.
            for value in rSep:
                rSep[value] = rSep[value] * 100.0

            for location in locations:
                # print("Location: " + str(location))
                # Radial position of each probe.
                rad_pos = {
                    'ad': geo.calc_R_measAD(r_probe, location),
                    'au': geo.calc_R_measAU(r_probe, location),
                    'bd': geo.calc_R_measBD(r_probe, location),
                    'bu': geo.calc_R_measBU(r_probe, location),
                    'cd': geo.calc_R_measCD(r_probe, location),
                    'cu': geo.calc_R_measCU(r_probe, location)
                }

                # Now to get the radial position at the omp for each probe.
                # Give this function a radial position and Z coordinate, and it will
                # return the psiN value (what flux surface it is on). Convert from cm
                # to m as well.
                psiN_AD = f_psiN((rad_pos['ad'] / 100.0), -0.18)
                psiN_AU = f_psiN(rad_pos['au'] / 100.0, -0.18)
                psiN_BD = f_psiN(rad_pos['bd'] / 100.0, -0.1546)
                psiN_BU = f_psiN(rad_pos['bu'] / 100.0, -0.1546)
                psiN_CD = f_psiN(rad_pos['cd'] / 100.0, -0.2054)
                psiN_CU = f_psiN(rad_pos['cu'] / 100.0, -0.2054)

                # Give this function a psiN value and z value (in this case the
                # z of the omp), and it will return the R value. Convert to cm.
                R_omp_AD = f_Romp(psiN_AD,
                                  np.around(Z_axis, decimals=3)) * 100.0
                R_omp_AU = f_Romp(np.around(psiN_AU, decimals=3), Z_axis) * 100.0
                R_omp_BD = f_Romp(psiN_BD, Z_axis) * 100.0
                R_omp_BU = f_Romp(psiN_BU, Z_axis) * 100.0
                R_omp_CD = f_Romp(psiN_CD, Z_axis) * 100.0
                R_omp_CU = f_Romp(psiN_CU, Z_axis) * 100.0

                # Now add the omp values into the rad_pos dictionary.
                rad_pos['ad_omp'] = R_omp_AD
                rad_pos['au_omp'] = R_omp_AU
                rad_pos['bd_omp'] = R_omp_BD
                rad_pos['bu_omp'] = R_omp_BU
                rad_pos['cd_omp'] = R_omp_CD
                rad_pos['cu_omp'] = R_omp_CU

                # Put into dictionary of lists containing rminrsep values for
                # each probe.
                rminrsep_values = {'ad': [], 'au': [], 'bd': [], 'bu': [], 'cd': [], 'cu': [],
                                   'ad_omp': [], 'au_omp': [], 'bd_omp': [], 'bu_omp': [],
                                   'cd_omp': [], 'cu_omp': []}
                rminrsep_values['ad'].append(rad_pos['ad'] - rSep['a'])
                rminrsep_values['au'].append(rad_pos['au'] - rSep['a'])
                rminrsep_values['bd'].append(rad_pos['bd'] - rSep['b'])
                rminrsep_values['bu'].append(rad_pos['bu'] - rSep['b'])
                rminrsep_values['cd'].append(rad_pos['cd'] - rSep['c'])
                rminrsep_values['cu'].append(rad_pos['cu'] - rSep['c'])

                # Add in R_omp - Rsep_omp, where it is translated to the outboard midplane.
                rminrsep_values['ad_omp'].append(rad_pos['ad_omp'] - rSep_omp)
                rminrsep_values['au_omp'].append(rad_pos['au_omp'] - rSep_omp)
                rminrsep_values['bd_omp'].append(rad_pos['bd_omp'] - rSep_omp)
                rminrsep_values['bu_omp'].append(rad_pos['bu_omp'] - rSep_omp)
                rminrsep_values['cd_omp'].append(rad_pos['cd_omp'] - rSep_omp)
                rminrsep_values['cu_omp'].append(rad_pos['cu_omp'] - rSep_omp)

                # Store in dictionary so these values can be used outside location loop.
                tmp_dict[str(shot)][str(time)][str(location)] = rminrsep_values

            # Next time step.
            time += step
            count += 1
            print("\r")

    # Create dictionaries to hold all the r-rep values and the averages.
    all_rminrsep = {}
    avg_rminrsep = {}
    probe_list = ['ad', 'au', 'bd', 'bu', 'cd', 'cu',
                  'ad_omp', 'au_omp', 'bd_omp', 'bu_omp', 'cd_omp', 'cu_omp']
    for probe in probe_list:
        all_rminrsep[probe] = {}
        avg_rminrsep[probe] = {}
        avg_rminrsep[probe + '_err'] = {}
        for loc in locations:
            all_rminrsep[probe][str(loc)] = []
            avg_rminrsep[probe][str(loc)] = 0
            avg_rminrsep[probe + '_err'][str(loc)] = 0

    # Place all the r-rsep values of a specific probe at a specific location
    # into a single list in all_rminrsep.
    for probe in probe_list:
        for loc in locations:
            for shot in shots:
                time = startTime
                while time <= endTime:
                    all_rminrsep[probe][str(loc)].append(tmp_dict[str(shot)][str(time)][str(loc)][probe])
                    time += step

    # Use all_rminrsep to get averages and std. devs.
    for probe in probe_list:
        for loc in locations:
            all_rminrsep[probe][str(loc)] = np.array(all_rminrsep[probe][str(loc)])
            tmp_avg = np.nanmean(all_rminrsep[probe][str(loc)])
            avg_rminrsep[probe][str(loc)] = tmp_avg
            tmp_std = np.nanstd(all_rminrsep[probe][str(loc)])
            avg_rminrsep[probe + '_err'][str(loc)] = tmp_std

    # Save to a txt file.
    if writeToFile:
        # filename = raw_input('Enter filename: ')
        f = open(filename, 'a')
        f.write("Shot: " + str(shot) + "\n")
        for key, value in avg_rminrsep.iteritems():
            f.write(key + " ")
            f.write(str(value))
            f.write("\n")
        f.write("\n")
        f.close()

    return avg_rminrsep
