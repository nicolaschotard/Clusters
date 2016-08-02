import yaml
import numpy as N

def load_config(config):
    return yaml.load(open(config))
    
def get_from_butler(butler, key, filt, patch, tract=0, dic=False):
    """Return selected data from a butler"""
    dataId = {'tract': tract, 'filter': filt, 'patch': patch}
    print dataId
    b = butler.get(key, dataId=dataId)
    print b
    return b if not dic else {n: b.get(n) for n in b.getSchema().getNames()}

def add_magnitudes(d, getMagnitude):
    Kfluxes = [k for k in d if k.endswith('_flux')]
    Ksigmas = [k+'Sigma' for k in fluxes]
    for kf, ks in zip(Kfluxes, Ksigmas):
        m, dm = N.array([getMagnitude(f, s) for f, s in zip(d[kf]. d[ks])]).T
        d[kf.replace('_flux', '_mag')] = m
        d[kf.replace('_fluxSigma', '_magSigma')] = dm

def add_position(d, wcs):
    d['x_src'], d['y_Src'] = N.array([wcs.skyToPixel(ra, dec)
                                      for ra, dec in zip(d["coord_ra"],
                                                         d["coord_dec"])]).T

def add_extra_info(d):

    # take the first filter, and the first patch
    f = d.keys()[0]
    p = d[f].keys()[0]

    # get the calib objects
    mag, wcs = d[f][p]['calexp'].getCalib().getmagnitude, d[f][p]['calexp'].getWcs()

    # compute all magnitudes
    for f in d:
        for p in d[f]:
            print "INFO:     adding magnitude for", f, p 
            add_magnitudes(d[f][p], mag)

    # compute all positions
    for f in d:
        for p in d[f]:
            print "INFO:     adding position for", f, p 
            add_position(d[f][p], wcs)
    
def get_all_data(path, patches, filters, add_extra=False):
    """
    Get butler data for a list of patches, for a list of filters
    Return a dictionnary with patches as keys
    """
    print "INFO: Loading data from", path, ", pathes:", patches, ", filters:", filters
    import lsst.daf.persistence as dafPersist
    butler = dafPersist.Butler(path)
    data = {f: get_filter_data(butler, path, patches, f) for f in filters}
    return data if not add_extra else add_extra_info(data)

def get_filter_data(butler, path, patches, f):
    """
    Get butler data for a list of patches, for a given filter
    Return a dictionnary with patches as keys
    """
    print "INFO: loading filter", f
    return {p: get_patch_data(butler, p, f) for p in patches}

def get_patch_data(butler, p, f):
    print "INFO:   loading patch", p
    meas = get_from_butler(butler, 'deepCoadd_meas', f, p, dic=True)
    forced = get_from_butler(butler, 'deepCoadd_forced_src', f, p, dic=True)
    calexp = get_from_butler(butler, 'deepCoadd_calexp', f, p, dic=False)
    return {'meas': meas, 'forced': forced, 'calexp':calexp}
    

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
