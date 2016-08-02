import yaml
import numpy as N
import lsst.afw.geom as afwGeom
from astropy.table import Table, Column, vstack

def load_config(config):
    """
    Load the configuration file, and return the corresponding dictionnary
    """
    return yaml.load(open(config))

def get_astropy_table(cat):
    """
    Convert an afw data table into a simple astropy table
    """
    schema = cat.getSchema() 
    dic = {n: cat.get(n) for n in schema.getNames()}
    tab = Table(dic)
    for k in schema.getNames():
        tab[k].description=schema[k].asField().getDoc()
        tab[k].unit=schema[k].asField().getUnits()
    return tab

def get_from_butler(butler, key, filt, patch, tract=0, table=False):
    """
    Return selected data from a butler for a given key, tract, patch and filter
    Either retrun the object or the astropy table version of it
    """
    dataId = {'tract': tract, 'filter': filt, 'patch': patch}
    b = butler.get(key, dataId=dataId)
    return b if not table else get_astropy_table(b)

def add_magnitudes(t, getMagnitude):
    """
    Compute magnitude for all fluxes of a given table and add the corresponding
    new columns
    """
    Kfluxes = [k for k in t.columns if k.endswith('_flux')]
    Ksigmas = [k+'Sigma' for k in Kfluxes]
    for kf, ks in zip(Kfluxes, Ksigmas):
        m, dm = N.array([getMagnitude(f, s) for f, s in zip(t[kf], t[ks])]).T
        t.add_columns([Column(name=kf.replace('_flux', '_mag'), data=m,
                              description='Magnitude', unit='mag'),
                       Column(name=ks.replace('_fluxSigma', '_magSigma'), data=dm,
                              description='Magnitude error', unit='mag')])

def add_position(t, wcs):
    """
    Compute the x/y position in pixel for all sources and add new columns to 
    the astropy table
    """
    x, y = N.array([wcs.skyToPixel(afwGeom.geomLib.Angle(ra), 
                                   afwGeom.geomLib.Angle(dec))
                    for ra, dec in zip(t["coord_ra"], t["coord_dec"])]).T
    t.add_columns([Column(name='x_Src', data=x,
                          description='x coordinate', unit='pixel'),
                   Column(name='y_Src', data=y,
                          description='y coordinate', unit='pixel')])

def add_filter_column(t, f):
    t.add_column(Column(name='filter', data=[f]*len(t), description='Filter name'))
    
def add_extra_info(d):
    """
    Add magnitude and position to all tables
    """
    # take the first filter, and the first patch
    f = d.keys()[0]
    p = d[f].keys()[0]

    # get the calib objects
    wcs = d[f][p]['calexp'].getWcs()
    getmag = d[f][p]['calexp'].getCalib().getMagnitude

    # redefine the magnitude function to make it 'work' for negative flux or sigma
    def mag(flux, sigma):
        if flux <= 0 or sigma <= 0:
            return N.nan, N.nan
        else:
            return getmag(flux, sigma)

    # compute all magnitudes and positions
    for f in d:
        for p in d[f]:
            for e in ['meas', 'forced']:
                print "INFO:     adding magnitude for", f, p, e
                add_magnitudes(d[f][p][e], mag)
                add_filter_column(d[f][p][e], f):
            print "INFO:     adding position for", f, p 
            add_position(d[f][p]['forced'], wcs)

    return d
    
def get_all_data(path, patches, filters, add_extra=False):
    """
    Get butler data for a list of patches and filters
    Return a dictionnary with filters as keys
    """
    print "INFO: Loading data from", path, " pathes:", patches, " filters:", filters
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
    """
    Get bulter data for a given set of patch and filter
    """
    print "INFO:   loading patch", p
    meas = get_from_butler(butler, 'deepCoadd_meas', f, p, table=True)
    forced = get_from_butler(butler, 'deepCoadd_forced_src', f, p, table=True)
    calexp = get_from_butler(butler, 'deepCoadd_calexp', f, p, table=False)
    return {'meas': meas, 'forced': forced, 'calexp':calexp}
    
def from_list_to_array(d):
    """
    Transform lists (of dict of list) into numpy arrays
    """
    if type(d) in [list, N.ndarray]:
        return N.array(d)
    for k in d:
        if type(d[k]) == list:
            d[k] = N.array(d[k])
        elif type(d[k]) == dict:
            from_list_to_array(d[k])
    return d

def stack_tables(d):
    """
    Stack the astropy tables across all patches
    Return a new dictionnary of the form:
    d = {u: 
          'forced': table,
          'meas': table,
         g: 
          'forced': table,
          'meas': table
         ...
        }
    """
    {'forced': vstack([vstack([d[f][p]['forced'] for p in d[f]]) for f in d])
    patches_stack = {f: {'forced': vstack([d[f][p]['forced'] for p in d[f]]),
                'meas': vstack([d[f][p]['meas'] for p in d[f]])}
            for f in d}
    filter_stack = {'forced': vstack([patches_stack[f]['forced'] for f in d])}
    mall=vstack([mr, mg])
    return {f: {'forced': vstack([d[f][p]['forced'] for p in d[f]]),
                'meas': vstack([d[f][p]['meas'] for p in d[f]])}
            for f in d}


def filter_table(t):

    # Select galaxies (and reject stars)
    filt = t['base_ClassificationExtendedness_flag'] == 0 # keep galaxy
    filt &= t['base_ClassificationExtendedness_value'] >= 0.5 # keep galaxy
    
            # Select sources which have a proper flux value in r, g and i bands
            # Notice that it would not be strictly necessary with forced photometry
            if N.any([forced[f][i].get(fluxFlagKey) for f in filters]):
                rejected['flag_flux'] += 1
                continue
    
            # Check the flux value, which must be > 0
            fluxes = {f: forced[f][i].get(fluxKey) for f in filters}
            fluxes_sigma = {f: forced[f][i].get(fluxSigmaKey) for f in filters}
            if any([fluxes[f] <= 0. for f in fluxes]):
                rejected['pos_flux'] += 1
                continue
    
            # Check the signal to noise (stn) value, which must be > 10
            stns = [forced[f][i].get(fluxKey)/forced[f][i].get(fluxSigmaKey) for f in filters
                    if forced[f][i].get(fluxSigmaKey) != 0]
            if any([stn < 10. for stn in stns]):
                rejected['stn'] += 1
                continue
    
            # Gauss regulerarization flag?
            if meas['r'][i].get(regaussFlagKey) or meas['r'][i].get(regaussFlagKey):
                rejected['gauss'] += 1
                continue
  
    return t[filt]

def keep_galaxies(table, key_colnames):
    if table['base_ClassificationExtendedness_flag'] == 0 \
       or ['base_ClassificationExtendedness_value'] < 0.5:
        return False
    else:
        return True

    
def select_data(data):

    # Initialize some lists
    print "INFO: Initializing variables"
    mags = {f: [] for f in filters}
    mags_sigma = {f: [] for f in filters}
    ell = {f: {'e1': [], 'e2': []} for f in filters}
    coords = {'ra': [], 'dec': [], 'id': []}
    resolution = {f: [] for f in filters}
    xSrc, ySrc = [], []
    # Loop over deblended sources in the filters forcedPhotCoadd catalogs
    ejected = {'star': 0, 'flag_flux':0, 'pos_flux': 0, 'stn': 0, 'gauss': 0}
    for i in range(len(forced[filters[0]])):
        
        # Select galaxies (and reject stars)
        if meas['r'][i].get(extFlagKey) or meas['r'][i].get(extKey) < 0.5:
            rejected['star'] += 1
            continue
        
        # Select sources which have a proper flux value in r, g and i bands
        # Notice that it would not be strictly necessary with forced photometry
        if N.any([forced[f][i].get(fluxFlagKey) for f in filters]):
            rejected['flag_flux'] += 1
            continue
        
        # Check the flux value, which must be > 0
        fluxes = {f: forced[f][i].get(fluxKey) for f in filters}
        fluxes_sigma = {f: forced[f][i].get(fluxSigmaKey) for f in filters}
        if any([fluxes[f] <= 0. for f in fluxes]):
            rejected['pos_flux'] += 1
            continue
        
        # Check the signal to noise (stn) value, which must be > 10
        stns = [forced[f][i].get(fluxKey)/forced[f][i].get(fluxSigmaKey) for f in filters
                if forced[f][i].get(fluxSigmaKey) != 0]
        if any([stn < 10. for stn in stns]):
            rejected['stn'] += 1
            continue
        
        # Gauss regulerarization flag?
        if meas['r'][i].get(regaussFlagKey) or meas['r'][i].get(regaussFlagKey):
            rejected['gauss'] += 1
            continue
        
        # Get filter dependent values
        for f in filters:
        
            # Need to use a calibobject in order to convert flux to magnitude
            m, sm = calib[f].getMagnitude(fluxes[f], fluxes_sigma[f])
            mags[f].append(m)
            mags_sigma[f].append(sm)
            
            # Get ellipticities
            ell[f]['e1'].append(meas[f][i].get(e1Key))
            ell[f]['e2'].append(meas[f][i].get(e2Key))
        
            # Get resolution
            resolution[f].append(meas[f][i].get(resKey))
        
        ra = forced[filters[0]][i].get(raKey)
        dec = forced[filters[0]][i].get(decKey)
        x, y = wcs.skyToPixel(ra, dec)
        coords['ra'].append(float(ra))
        coords['dec'].append(float(dec))
        coords['id'].append(int(meas[filters[0]][i].get("id")))
        xSrc.append(x)
        ySrc.append(y)
