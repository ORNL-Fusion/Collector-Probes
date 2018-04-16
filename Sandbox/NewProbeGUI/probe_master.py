import pandas as pd
import numpy  as np
import create_df
import pull_mdsplus as pull


class Probe:
    """
    Probe class to hold the rbs data and LAMS data and to save it all.
    """
    def __init__(self, letter, number, start, end, step):
        self.number     = int(number)
        self.letter     = letter
        self.time_start = start
        self.time_end   = end
        self.time_step  = step

    def get_rbs(self, remote=True):
        self.rbs_data = create_df.get_rbs(self.number, self.letter, self.time_start,
                                          self.time_end, self.time_step, remote, True)
    def get_lams(self, remote=True):
        if remote:
            input("SSH link r2d2 to localhost. Press enter to continue...")
            conn = pull.thin_connect(self.number, server='localhost')
        lams_dict_U    = pull.pull_lams(conn, self.number, self.letter + 'U', verbal=True)
        lams_dict_D    = pull.pull_lams(conn, self.number, self.letter + 'D', verbal=True)
        lams_df_U      = pd.DataFrame(lams_dict_U)
        lams_df_D      = pd.DataFrame(lams_dict_D)
        lams_df        = pd.concat((lams_df_U, lams_df_D), axis=1)
        self.lams_data = lams_df

    def save_excel(self, filename):
        writer = pd.ExcelWriter(filename)
        self.rbs_data.to_excel(writer, 'RBS Data')
        self.lams_data.to_excel(writer, 'LAMS Data')
        writer.save()

    def to_hdf5(self, filename):
        pass

def get_probe_data(probes_to_get, time_start=2500, time_end=5000, time_step=500, remote=True):
    """
    Function to pull probes and return the DataFrames in a list.
    """

    p_list = []
    for probe in probes_to_get:
        p = Probe(probe[0], probe[1:], time_start, time_end, time_step)
        p.get_rbs(remote)
        p.get_lams(remote)
        p_list.append(p)

    ans = input("Save to Excel files (y/n)? ")
    if ans == 'y':
        for p in p_list:
            p.save_excel(p.letter + str(p.number) + '.xlsx')

    return p_list
