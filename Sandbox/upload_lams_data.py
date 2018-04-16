import MDSplus as mds
import pandas  as pd
import numpy   as np


print("Will upload the LAMS csv files to MDSplus on R2D2.")
print("  p_type: One of AD, AU, BD, BU, CD, CU")
print("  conn:   An MDSplus Connection to R2D2. Optional.")
print("  p_num:  Which probe number.")

# Currently can't get this to work with anything other than 2.
def upload_thin(p_num, p_type='AD', conn=None):
    """
    Upload the csv files of the LAMS data to the R2D2 MDSplus repository.

    conn: An MDSplus connection to R2D2.
    """

    if conn is None:
        conn = mds.Connection('localhost')

    # Open the tree.
    conn.openTree('dp_probes', p_num)

    if p_num < 10:
        p_num = '0' + str(p_num)
    else:
        p_num = str(p_num)

    # Pull in the LAMS data into a pandas DataFrame.
    df = pd.read_csv('CP_LAMS_Dictionary/' + p_type + p_num + '_dict.csv')

    # Put in ScanDist.
    path = '\\DP_PROBES::TOP.' + p_type[0] + '.' + p_type + '.LAMS:SCANDIST'
    print(path)
    data = mds.Float64Array(df['ScanDist'])
    conn.put(path, "$", data)

    for iso in ['180', '182', '183', '184', '186']:
        path = '\\DP_PROBES::TOP.' + p_type[0] + '.' + p_type + '.LAMS:EF' + iso
        data = mds.Float64Array(df['EF' + iso])
        conn.put(path, '$', data)
        data = mds.Float64Array(df['EF' + iso + 'err'])
        conn.put(path + '_ERR', '$', data)

    conn.closeAllTrees()

def upload_thicc(p_num, p_type):
    """
    Require a thick connection. If connecting via ssh tunneling, this just
    means putting the line:

    export dp_probes_path='localhost::/mnt/syn/MDSplus/dp_probes_data'

    in your bashrc or bash_profile file. And then making sure you have an
    ssh tunnel open to R2D2.
    """

    # Grab the tree.
    tree = mds.Tree('dp_probes', p_num)

    # Put number into correct string format.
    if p_num < 10:
        p_num = '0' + str(p_num)
    else:
        p_num = str(p_num)

    # Pull in the LAMS data into a pandas DataFrame.
    df = pd.read_csv('CP_LAMS_Dictionary/' + p_type + p_num + '_dict.csv')

    # Put in ScanDist.
    path = '\\DP_PROBES::TOP.' + p_type[0] + '.' + p_type + '.LAMS:SCANDIST'
    node = tree.getNode(path)
    data = mds.Float64Array(df['ScanDist'])
    node.putData(data)

    # Put in the EF's.
    for iso in ['180', '182', '183', '184', '186']:
        path = '\\DP_PROBES::TOP.' + p_type[0] + '.' + p_type + '.LAMS:EF' + iso
        data = mds.Float64Array(df['EF' + iso])
        node = tree.getNode(path)
        node.putData(data)
        data = mds.Float64Array(df['EF' + iso + 'err'])
        node = tree.getNode(path)
        node.putData(data)

    tree.write()
