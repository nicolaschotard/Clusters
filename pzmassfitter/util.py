"""Some useful classes or functions used in pzmassfitter."""


from __future__ import print_function
import re
import numpy as np


class VarContainer(dict):
    """Basic Python Object that is a dual between an object and a dictionary."""
    def __getattr__(self, name):
        return self[name]

    def __setattr__(self, name, value):
        self[name] = value


def readtxtfile(filename):
    """Read a text file."""
    infile = open(filename)
    if infile is None:
        return None

    try:
        rows = []
        for line in infile:
            if not re.match('^#', line):
                tokens = line.split()
                if len(tokens) != 0:
                    rows.append([float(l) for l in line.split()])
        return np.row_stack(rows)
    except ValueError as error:
        print('Cannot convert to floats, returning list')
        print(error)
        rows = []
        infile = open(filename)
        for line in infile:
            if not re.match('^#', line):
                rows.append(line.split())
        return rows


def loadchains(chainfilenames, trim=False):
    """Load Output chains from MyMC."""
    chainfiles = [readtxtfile(x) for x in chainfilenames]

    takelength = len(chainfiles[0])
    if trim is True:
        takelength = np.min(np.array([len(chainfile) for chainfile in chainfiles]))

    rawdata = [np.row_stack([[float(yy) for yy in y] for y in x[1:takelength]]) for x in chainfiles]

    columns = chainfiles[0][0]

    chain = {}

    for i, col in enumerate(columns):
        chain[col] = np.row_stack([x[:, i] for x in rawdata])

    return chain


def matchById(firstcat, othercat, otherid='SeqNr', selfid='SeqNr'):
    """Returns a subset of the first catalog, that matches the order of the second catalog."""
    order = {}
    for i, fid in enumerate(firstcat[selfid]):
        order[fid] = i

    keeporder = []
    for oid in othercat[otherid]:
        if oid in order:
            keeporder.append(order[oid])

    keep = np.array(keeporder)
    matched = firstcat[keep]
    print("INFO: %i matched galaxies kept" % len(matched))
    return matched
