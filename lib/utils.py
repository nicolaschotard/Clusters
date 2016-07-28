import numpy as N

def get_from_butler(butler, key, filt, patch, tract=0):
    """Return selected data from a butler"""
    dataId = {'tract': tract, 'filter': filt, 'patch': patch}
    return butler.get(key, dataId=dataId)

def from_list_to_array(d):
    """Transform lists (of dict of list) into numpy arrays"""
    if type(d) in [list, N.ndarray]:
        return N.array(d)
    for k in d:
        if type(d[k]) == list:
            d[k] = N.array(d[k])
        elif type(d[k]) == dict:
            from_list_to_array(d[k])
    return d
