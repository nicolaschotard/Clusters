import re
from numpy import row_stack


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
        return row_stack(rows)
    except ValueError, e:
        print 'Cannot convert to floats, returning list'
        print e
        rows = []
        infile = open(filename)
        for line in infile:
            if not re.match('^#', line):
                rows.append(line.split())
        return rows
