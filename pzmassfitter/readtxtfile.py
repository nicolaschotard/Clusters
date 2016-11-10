import re
from numpy import *

########################################

def readtxtfile(filename):

    input = open(filename);
    if input is None:
        return None

    try:
        rows = []
        for line in input:
            if not re.match('^#', line):
                tokens = line.split()
                if len(tokens) != 0:
                    rows.append(map(float, line.split()))
        return row_stack(rows)
    except ValueError,e:
        print 'Cannot convert to floats, returning list'
        print e
        rows = []
        input = open(filename)
        for line in input:
            if not re.match('^#', line):
                rows.append(line.split())            
        return rows

#####################################
