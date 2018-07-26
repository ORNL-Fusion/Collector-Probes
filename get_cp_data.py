import pandas as pd
import numpy  as np
import create_df
import pull_mdsplus as pull


print("\nScript related to retrieving collector probe data from r2d2 and atlas.")
print("and atlas. This requires an account on Cybele/Iris as well as R2D2. Try")
print("'help(get_cp_data.get_probe_data)' for info and an example.")

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
        self.lams_avail = True
        self.rbs_avail  = True

    def get_rbs_and_lams(self, remote=True):
        """
        Pull all the RBS relevant data stored on r2d2 and return it as a DataFrame.
        """
        self.rbs_data, self.lams_data = create_df.get_rbs_and_lams(self.number,
                                                                   self.letter,
                                                                   self.time_start,
                                                                   self.time_end,
                                                                   self.time_step,
                                                                   remote,
                                                                   True)

    def save_excel(self, filename):
        """
        Save all the RBS and LAMS data into a Excel workbook with each on its
        own sheet.
        """
        writer = pd.ExcelWriter(filename)
        if self.rbs_data is not None:
            self.rbs_data.to_excel(writer, 'RBS Data')
        if self.lams_data is not None:
            self.lams_data.to_excel(writer, 'LAMS Data')
        writer.save()

    def to_hdf5(self, filename):
        """
        Not yet implemented.
        """
        pass

def get_probe_data(probes_to_get, time_start=2500, time_end=5000, time_step=500, remote=True):
    """
    Function to pull probes and return the Probe class objects in a list. The
    user can either use the defined functions in the Probe class above (such as
    saving in Excel or hdf5), or use the data in their own way.

    probes_to_get:
        List of strings indicating which probes to retrieve. An example input
        could be: probes_to_get = ['A15', 'B8', 'C8'].
    time_start:
        The beginning of the time range for which the average value of Rsep
        (the radial position of the separatrix), psin and Rsep omp will be used.
        A good time range would be during the flat top.
    time_end:
        End time of the above.
    time_step:
        Time step for the above. Having five or so times seems to work fine.
    remote:
        Set to True if accessing outside DIII-D network. This requires ssh linking
        either r2d2 or atlas to your localhost. An example way of doing this is
        in a new terminal try:
            ssh -Y -p 2039 -L 8000:r2d2.gat.com:8000 username@cybele.gat.com
        for either r2d2 or atlas, as required. Once you enter your Cybele
        password at the prompt, switch back over and continue this script.

    Example session:
        import get_cp_data as cp
        myprobes = ['A15', 'B8', 'C8']
        mylist = cp.get_probe_data(myprobes, time_start=3000, time_end=4500,
                                   time_step=500, remote=True)

        # You will be prompted to save the data to Excel. If you wish to save
        # after the fact

        for myprobe in mylist:
            myprobe.save_excel(myprobe.letter + str(myprobe.number) + '.xlsx')

        # If you want the probe data to use as your own in the terminal...

        my_rbs_data = mylist[0].rbs_data
        my_lams_data = mylist[0].lams_data


    """

    # The list to hold the returned Probe objects.
    p_list = []
    for probe in probes_to_get:
        p = Probe(probe[0], probe[1:], time_start, time_end, time_step)
        try:
            p.get_rbs_and_lams(remote)
        except:
            print("Error retrieving data.")
        p_list.append(p)

    # Ask if user wants to save data in an Excel file.
    ans = input("Save to Excel files (y/n)? ")
    if ans == 'y':
        try:
            for p in p_list:
                p.save_excel(p.letter + str(p.number) + '.xlsx')
        except Exception as e:
            print("Error saving probe data: " + p.letter + str(p.number))
            print("Error message: " + e)

    # Return the list of Probe objects for the user to use in a terminal.
    return p_list
