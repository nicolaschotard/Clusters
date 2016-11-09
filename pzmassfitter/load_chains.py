#########
# Load Output chains from MyMC
#########

import readtxtfile
import numpy as np

def loadChains(chainfilenames, trim=False):

    chainfiles = [readtxtfile.readtxtfile(x) for x in chainfilenames]

    takelength = len(chainfiles[0])
    if trim is True:
        takelength = np.min(np.array([len(chainfile) for chainfile in chainfiles]))

    rawdata = [np.row_stack([map(float, y) for y in x[1:takelength]]) for x in chainfiles]

    columns = chainfiles[0][0]

    chain = {}

    for i, col in enumerate(columns):
        chain[col] = np.row_stack([x[:,i] for x in rawdata])

    return chain



