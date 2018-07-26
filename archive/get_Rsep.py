# In separate terminal use
# ssh -Y -p 2039 -L 8000:atlas.gat.com:8000 username@cybele.gat.com
# where you substitute in your username. The MDSplus connection is
# then MDSplus.Connection('localhost').
from __future__ import print_function
import numpy as np
import scipy.interpolate as scinter

# Normal load_gfile_d3d except needed to cast 'g['lcfs']' as an int.
import MDSplus as mds
#import EFIT.load_gfile_d3d as loadg
import meas_locations as geo


# Copy and pasted from the EFOT repository for convienence.
def read_g_file_mds(shot, time, tree='EFIT01', exact=False, connection=None,
                    write2file=True, gpath='.'):
    import MDSplus

    # in case those are passed in as strings
    shot = int(shot)
    time = int(time)

    # Connect to server, open tree and go to g-file
    if connection is None:
        connection = MDSplus.Connection('atlas.gat.com')
    connection.openTree(tree, shot)
    base = 'RESULTS:GEQDSK:'

    # get time slice
    signal = 'GTIME'
    k = np.argmin(np.abs(connection.get(base + signal).data() - time))
    time0 = int(connection.get(base + signal).data()[k])

    if (time != time0):
        if exact:
            raise RuntimeError(tree + ' does not exactly contain time ' + str(time) + '  ->  Abort')
        else:
            print('Warning: ' + tree + ' does not exactly contain time ' + str(time) + ' the closest time is ' + str(time0))
            print('Fetching time slice ' + str(time0))
            time = time0

    # store data in dictionary
    g = {'shot': shot, 'time': time}

    # get header line
    header = connection.get(base + 'ECASE').data()[k]

    # get all signals, use same names as in read_g_file
    translate = {'MW': 'NR', 'MH': 'NZ', 'XDIM': 'Xdim', 'ZDIM': 'Zdim', 'RZERO': 'R0',
                 'RMAXIS': 'RmAxis', 'ZMAXIS': 'ZmAxis', 'SSIMAG': 'psiAxis', 'SSIBRY': 'psiSep',
                 'BCENTR': 'Bt0', 'CPASMA': 'Ip', 'FPOL': 'Fpol', 'PRES': 'Pres',
                 'FFPRIM': 'FFprime', 'PPRIME': 'Pprime', 'PSIRZ': 'psiRZ', 'QPSI': 'qpsi',
                 'NBBBS': 'Nlcfs', 'LIMITR': 'Nwall'}
    for signal in translate:
        g[translate[signal]] = connection.get(base + signal).data()[k]

    g['R1'] = connection.get(base + 'RGRID').data()[0]
    g['Zmid'] = 0.0

    RLIM = connection.get(base + 'LIM').data()[:, 0]
    ZLIM = connection.get(base + 'LIM').data()[:, 1]
    g['wall'] = np.vstack((RLIM, ZLIM)).T

    RBBBS = connection.get(base + 'RBBBS').data()[k][:int(g['Nlcfs'])]
    ZBBBS = connection.get(base + 'ZBBBS').data()[k][:int(g['Nlcfs'])]
    g['lcfs'] = np.vstack((RBBBS, ZBBBS)).T

    KVTOR = 0
    RVTOR = 1.7
    NMASS = 0
    RHOVN = connection.get(base + 'RHOVN').data()[k]

    # convert floats to integers
    for item in ['NR', 'NZ', 'Nlcfs', 'Nwall']:
        g[item] = int(g[item])

    # convert single (float32) to double (float64) and round
    for item in ['Xdim', 'Zdim', 'R0', 'R1', 'RmAxis', 'ZmAxis', 'psiAxis', 'psiSep', 'Bt0', 'Ip']:
        g[item] = np.round(np.float64(g[item]), 7)

    # convert single arrays (float32) to double arrays (float64)
    for item in ['Fpol', 'Pres', 'FFprime', 'Pprime', 'psiRZ', 'qpsi', 'lcfs', 'wall']:
        g[item] = np.array(g[item], dtype=np.float64)

    # write g-file to disk
    if write2file:
        if not (gpath[-1] == '/'): gpath += '/'
        with open(gpath + 'g' + format(shot,'06d') + '.' + format(time,'05d'), 'w') as f:
            if ('EFITD' in header[0]) and (len(header) == 6):
                for item in header: f.write(item)
            else:
                f.write('  EFITD    xx/xx/xxxx    #' + str(shot) + '  ' + str(time) + 'ms        ')

            f.write('   3 ' + str(g['NR']) + ' ' + str(g['NZ']) + '\n')
            f.write('% .9E% .9E% .9E% .9E% .9E\n'%(g['Xdim'], g['Zdim'], g['R0'], g['R1'], g['Zmid']))
            f.write('% .9E% .9E% .9E% .9E% .9E\n'%(g['RmAxis'], g['ZmAxis'], g['psiAxis'], g['psiSep'], g['Bt0']))
            f.write('% .9E% .9E% .9E% .9E% .9E\n'%(g['Ip'], 0, 0, 0, 0))
            f.write('% .9E% .9E% .9E% .9E% .9E\n'%(0,0,0,0,0))
            write_array(g['Fpol'], f)
            write_array(g['Pres'], f)
            write_array(g['FFprime'], f)
            write_array(g['Pprime'], f)
            write_array(g['psiRZ'].flatten(), f)
            write_array(g['qpsi'], f)
            f.write(str(g['Nlcfs']) + ' ' + str(g['Nwall']) + '\n')
            write_array(g['lcfs'].flatten(), f)
            write_array(g['wall'].flatten(), f)
            f.write(str(KVTOR) + ' ' + format(RVTOR, ' .9E') + ' ' + str(NMASS) + '\n')
            write_array(RHOVN, f)

    # Construct (R,Z) grid for psiRZ
    g['dR'] = g['Xdim']/(g['NR'] - 1)
    g['R'] = g['R1'] + np.arange(g['NR'])*g['dR']

    g['dZ'] = g['Zdim']/(g['NZ'] - 1)
    NZ2 = int(np.floor(0.5*g['NZ']))
    g['Z'] = g['Zmid'] + np.arange(-NZ2, NZ2+1)*g['dZ']

    # normalize psiRZ
    g['psiRZn'] = (g['psiRZ'] - g['psiAxis']) / (g['psiSep'] - g['psiAxis'])

    return g


# Return dictionary of average R - Rsep value for each of the probes.
#
# shots     = list of shots in for. If only one shot still enter as list.
# r_probe   = radial position of probe holder tip (cm).
# location  = location along probe where you want the R - Rsep value (cm)
# startTime = start of time range to be averaged over (ms).
# endTime   = end of time range (ms).

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

    print("Shots to be analyzed: ", end='')
    for shot in shots:
        print(str(shot), end='')
        if (len(shots) > 1):
            print(", ", end='')
    print("\n")

    total_runs = float(len(shots) * (endTime-startTime)/step)
    count = 0

    for shot in shots:
        time = startTime
        while time <= endTime:

            #print("Shot:     " + str(shot), end='\n')
            # print("Location: " + str(location))
            #print("Time:     " + str(time), end='\n')

            # The "\033[F" means go up to the previous line. This is all so there
            # can be a progress bar and not spew a bunch of stuff at u.
            if count != 0:
                print("\033[F \033[F \033[F \033[F \033[F", end='')
            print("\rCurrent Shot: ", end='')
            print(str(shot), end='')
            print("  ", end='')
            print("Time: ", end='')
            print(str(time))

            print("\r[", end='')
            print(u"\u2588"*int((count/total_runs)*15.0), end='')
            print(" "*int(15-((count/total_runs)*15.0)-1), end='')
            print("] " + str(int(count/total_runs*100.0))+"% \n\n")

            # Lines from Zeke's code.
            #parmDICT = loadg.read_g_file_mds(shot, time,
            #                                 connection=MDSplusCONN,
            #                                 write2file=False,
            #                                 tree=Etree)
            parmDICT = read_g_file_mds(shot, time,
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
            rSep['a'] = f_Rs(-0.188)
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
                                   'cd_omp': [], 'cu_omp': [],
                                   'ad_psiN': [], 'au_psiN': [], 'bd_psiN': [], 'bu_psiN': [],
                                   'cd_psiN': [], 'cu_psiN': []}
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

                # Also store in the psiN values.
                rminrsep_values["ad_psiN"].append(psiN_AD)
                rminrsep_values["au_psiN"].append(psiN_AU)
                rminrsep_values["bd_psiN"].append(psiN_BD)
                rminrsep_values["bu_psiN"].append(psiN_BU)
                rminrsep_values["cd_psiN"].append(psiN_CD)
                rminrsep_values["cu_psiN"].append(psiN_CU)

                # Store in dictionary so these values can be used outside location loop.
                tmp_dict[str(shot)][str(time)][str(location)] = rminrsep_values

            # Next time step.
            time += step
            count += 1
            print("\r")
        #print("\n\n")
    # Create dictionaries to hold all the r-rsep values and the averages.
    all_rminrsep = {}
    avg_rminrsep = {}
    probe_list = ['ad', 'au', 'bd', 'bu', 'cd', 'cu',
                  'ad_omp', 'au_omp', 'bd_omp', 'bu_omp', 'cd_omp', 'cu_omp',
                  'ad_psiN', 'au_psiN', 'bd_psiN', 'bu_psiN', 'cd_psiN', 'cu_psiN']
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
