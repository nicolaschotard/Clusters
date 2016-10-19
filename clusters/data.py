"""Data builder and parser for the Clusters package."""

import os
import sys
import yaml
import numpy
from astropy.wcs import WCS, utils
from astropy.coordinates import SkyCoord, Angle
from astropy.table import Table, Column, vstack
from progressbar import Bar, ProgressBar, Percentage, ETA


def load_config(config):
    """Load the configuration file, and return the corresponding dictionnary.

    :param config: Name of the configuration file.
    :type config: str.
    :returns: the configuration elements in a python dictionnary
    """
    return yaml.load(open(config))


def shorten(doc):
    """Hack to go around an astropy/hdf5 bug. Cut in half words longer than 18 chars."""
    return " ".join([w if len(w) < 18 else (w[:len(w) / 2] + ' - ' + w[len(w) / 2:])
                     for w in doc.split()])


def get_astropy_table(cat, **kwargs):
    """Convert an afw data table into a simple astropy table.

    :param cat: an afw data table
    :return: the corresponding astropy.table.Table
    """
    schema = cat.getSchema()
    dic = cat.getColumnView().extract(*kwargs['keys'] if 'keys' in kwargs else "*")
    tab = Table(dic)
    for k in tab.keys():
        tab[k].description = shorten(schema[k].asField().getDoc())
        tab[k].unit = schema[k].asField().getUnits()
    return tab


def get_from_butler(butler, catalog, filt, patch, **kwargs):
    """
    Return selected data from a butler for a given catalog, tract, patch and filter.

    :param string butler: a butler instance
    :param string catalog: name of a catalog
    :param string filter: name a a filter
    :param string patch: name of a patch

    Possible kwargs are

    :param int tract: tract in which to look for data (default will be 0)
    :param bool table: If True, an astropy tab le will be returned (default is False)
    :param list keys: a list of keys to get from the catalog

    Either retrun the object or the astropy table version of it
    """
    tract = 0 if 'tract' not in kwargs else kwargs['tract']
    table = False if 'table' not in kwargs else kwargs['table']
    dataid = {'tract': tract, 'filter': filt, 'patch': patch}
    afwt = butler.get(catalog, dataId=dataid)
    return afwt if not table else get_astropy_table(afwt, **kwargs)


def add_magnitudes(table, getmagnitude):
    """Compute magns for all fluxes of a given table. Add the corresponding new columns."""
    kfluxes = [k for k in table.columns if k.endswith('_flux')]
    ksigmas = [k + 'Sigma' for k in kfluxes]
    for kflux, ksigma in zip(kfluxes, ksigmas):
        mag, dmag = numpy.array([getmagnitude(f, s)
                                 for f, s in zip(table[kflux], table[ksigma])]).T
        table.add_columns([Column(name=kflux.replace('_flux', '_mag'), data=mag,
                                  description='Magnitude', unit='mag'),
                           Column(name=ksigma.replace('_fluxSigma', '_magSigma'), data=dmag,
                                  description='Magnitude error', unit='mag')])


def add_position_and_deg(table, wcs):
    """Compute the x/y position in pixel for all sources. Add new columns to the table."""
    # Add the x / y position in pixel
    xsrc, ysrc = skycoord_to_pixel(SkyCoord(table["coord_ra"].tolist(),
                                            table["coord_dec"].tolist(), unit='rad'), wcs)
    table.add_columns([Column(name='x_Src', data=xsrc,
                              description='x coordinate', unit='pixel'),
                       Column(name='y_Src', data=ysrc,
                              description='y coordinate', unit='pixel')])

    # Add a new column to have to coordinates in degree
    table.add_columns([Column(name='coord_ra_deg',
                              data=Angle(table['coord_ra'].tolist(), unit='rad').degree,
                              description='RA coordinate', unit='degree'),
                       Column(name='coord_dec_deg',
                              data=Angle(table['coord_dec'].tolist(), unit='rad').degree,
                              description='DEC coordinate', unit='degree')])


def add_filter_column(table, filt):
    """Add a new column containing the filter name."""
    table.add_column(Column(name='filter', data=[filt] * len(table), description='Filter name'))


def add_patch_column(table, patch):
    """Add a new column containing the patch name."""
    table.add_column(Column(name='patch', data=[patch] * len(table), description='Patch name'))


def add_extra_info(data):
    """Add magnitude and position to all tables."""
    # get the wcs
    wcs = WCS(get_wcs(data))

    # Shorcut to magnitude function
    filt = data.keys()[0]
    patch = data[filt].keys()[0]
    getmag = data[filt][patch]['calexp'].getCalib().getMagnitude

    def mag(flux, sigma):
        """Redefine the magnitude function. Negative flux or sigma possible."""
        if flux <= 0 or sigma <= 0:
            return numpy.nan, numpy.nan
        else:
            return getmag(flux, sigma)

    # compute all magnitudes and positions
    for filt in data:  # loop on filters
        for patch in data[filt]:  # loop on patches
            for cat in ['deepCoadd_meas', 'deepCoadd_forced_src']:  # loop on catalogs
                print "INFO: adding extra info for", filt, patch, cat
                add_magnitudes(data[filt][patch][cat], mag)
                add_filter_column(data[filt][patch][cat], filt)
                add_patch_column(data[filt][patch][cat], patch)
                add_position_and_deg(data[filt][patch][cat], wcs)

    return data


def get_wcs(d):
    """Get the wcs dictionnary from the butler."""
    # take the first filter, and the first patch
    filt = d.keys()[0]
    patch = d[filt].keys()[0]
    return d[filt][patch]['calexp'].getWcs().getFitsMetadata().toDict()


def save_wcs(wcs, output):
    """Save the wcs dictionnary into a valid astropy Table format."""
    table = Table({k: [wcs[k]] for k in wcs})
    table.write(output, path='wcs', compression=True,
                append=True, serialize_meta=True)


def load_wcs(wcs):
    """Get back the right wcs format from the hdf5 table."""
    return WCS({k: wcs[k].item() for k in wcs.keys()})


def skycoord_to_pixel(coords, wcs, unit='deg'):
    """Transform sky coordinates (ra, dec) to pixel coordinates (x, y) given a wcs.

    :param coords: Coordinates. Multiple formats accepted:

     - [ra, dec]
     - [[ra1, ra2], [dec1, dec2]]
     - or a SkyCoord object

    :param wcs: an astropy.wcs.WCS object
    :return: A list of (x, y) coordinates in pixel units
    """
    if not isinstance(coords, SkyCoord):
        coords = SkyCoord(coords[0], coords[1], unit=unit)
    return utils.skycoord_to_pixel(coords, wcs)


def pixel_to_skycoord(xsrc, ysrc, wcs):
    """Transform pixel coordinates (x, y) to sky coordinates (ra, dec in deg) given a wcs.

    :param float x: x coordinate
    :param float y: y coordinate
    :param wcs: an astropy.wcs.WCS object
    :return: an astropy.coordinates.SkyCoord object.
    """
    return utils.pixel_to_skycoord(xsrc, ysrc, wcs)


def get_all_data(path, patches, filters, add_extra=False, keys=None, show=False):
    """
    Get butler data for a list of patches and filters.

    Return a dictionnary with filters as keys
    """
    import lsst.daf.persistence as dafPersist
    if show:
        patches = [patches[0]]
        filters = filters[0]
    else:
        print "INFO: Loading data from", path, " pathes:", patches, " filters:", filters
    butler = dafPersist.Butler(path)
    d = {f: get_filter_data(butler, patches, f, keys=keys) for f in filters}
    out = stack_tables(d) if not add_extra else stack_tables(add_extra_info(d))
    if show:
        print "INFO: Available list of keys for the deepCoadd_forced_src catalog"
        table = Table(numpy.transpose([[k, out['deepCoadd_forced_src'][k].description,
                                        out['deepCoadd_forced_src'][k].unit] for k in
                                       sorted(out['deepCoadd_forced_src'].keys())]).tolist(),
                      names=["Keys", "Description", "Units"])
        print " -> All saved in deepCoadd_forced_src_keys.txt"
        table.write("deepCoadd_forced_src_keys.txt", format='ascii', comment="#")

        print "INFO: Available list of keys for the deepCoadd_meas catalog"
        table = Table(numpy.transpose([[k, out['deepCoadd_meas'][k].description,
                                        out['deepCoadd_meas'][k].unit]
                                       for k in sorted(out['deepCoadd_meas'].keys())]).tolist(),
                      names=["Keys", "Description", "Units"])
        print " -> All saved in deepCoadd_meas_keys.txt"
        table.write("deepCoadd_meas_keys.txt", format='ascii', comment="#")

        sys.exit()

    out['wcs'] = get_wcs(d)
    return out


def get_filter_data(butler, patches, filt, keys=None):
    """
    Get butler data for a list of patches, for a given filter.

    Return a dictionnary with patches as keys
    """
    print "INFO: loading filter", filt
    return {patch: get_patch_data(butler, patch, filt, keys=keys) for patch in patches}


def get_patch_data(butler, p, f, keys=None):
    """Get bulter data for a given set of patch and filter."""
    keys = {} if keys is None else keys
    print "INFO:   loading patch", p
    mkeys = "*" if 'deepCoadd_meas' not in keys else keys['deepCoadd_meas']
    fkeys = "*" if 'deepCoadd_forced_src' not in keys else keys['deepCoadd_forced_src']
    meas = get_from_butler(butler, 'deepCoadd_meas', f, p, table=True, keys=mkeys)
    forced = get_from_butler(butler, 'deepCoadd_forced_src', f, p, table=True, keys=fkeys)
    calexp = get_from_butler(butler, 'deepCoadd_calexp', f, p, table=False)
    return {'deepCoadd_meas': meas, 'deepCoadd_forced_src': forced, 'calexp': calexp}


def merge_dicts(*dict_args):
    """Merge two dictionnary.

    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def get_ccd_data(path, save=False, keys="*", show=False):
    """Get the 'forced_src' catalogs.

    .. todo::

     - add the list of needed keys in the config file for the following catalog:

       - forced_src
       - deepCoadd_meas
       - deepCoadd_forced_src

     - clean this function and add it to clusters_data
    """
    import lsst.daf.persistence as dafPersist
    butler = dafPersist.Butler(path)
    idkeys = ['ccd', 'date', 'filter', 'object', 'runId', 'visit']
    print "INFO: Getting list of available data"
    dids = [merge_dicts(dict(zip(idkeys, v)), {'tract': 0})
            for v in butler.queryMetadata("forced_src", format=idkeys)]
    dids = [d for d in dids if os.path.exists(path + "/forced/%s/%s/%s/%s/%s/FORCEDSRC-%i-%i.fits" %
                                              (d['runId'], d['object'], d['date'], d['filter'],
                                               d['tract'], d['visit'], d['ccd']))]
    if show:
        print "INFO: Available list of keys for the forced_src catalog"
        table = get_astropy_table(butler.get("forced_src", dataId=dids[0]), keys=keys)
        table = Table(numpy.transpose([[k, table[k].description, table[k].unit]
                                       for k in sorted(table.keys())]).tolist(),
                      names=["Keys", "Description", "Units"])
        print " -> All saved in forced_src_keys.txt"
        table.write("forced_src_keys.txt", format='ascii', comment="#")
        sys.exit()

    print "INFO: Getting the data from the butler for %i fits files" % len(dids)
    pbar = ProgressBar(widgets=[Percentage(), Bar(), ETA()], maxval=len(dids)).start()

    def get_tables(did, i):
        """Get a table and add a few keys."""
        table = get_astropy_table(butler.get("forced_src", dataId=did), keys=keys)
        for key in did:
            table.add_column(Column(name=key, data=[did[key]] * len(table)))
        pbar.update(i + 1)
        return table

    tables = [get_tables(did, i) for i, did in enumerate(dids)]
    pbar.finish()
    if save:
        save_tables(tables)
    return tables


def save_tables(tables):
    """Save a list of astropy tables in an hdf5 file."""
    pbar = ProgressBar(widgets=[Percentage(), Bar(), ETA()], maxval=len(tables)).start()
    for i, table in enumerate(tables):
        if i == 0:
            table.write('forced_src.hdf5', path='%i' % i, compression=True,
                        serialize_meta=True, overwrite=True)
        else:
            table.write('forced_src.hdf5', path='%i' % i, compression=True,
                        serialize_meta=True, append=True)
        pbar.update(i + 1)
    pbar.finish()


def filter_and_stack_tables(tables, **kwargs):
    """Apply filter on a list of astropy table and vertically stack the resulting table list."""
    if len(kwargs) == 0:
        raise IOError("You should at least give one filter, e.g. **{'objectId': 2199694371871}")
    for k in kwargs:
        tables = [tables[i][tables[i][k] == kwargs[k]] for i in range(len(tables))]
    return vstack(tables)


def vstack2(tables):
    """Verticaly stack large amount of astropy tables."""
    pbar = ProgressBar(widgets=[Percentage(), Bar(), ETA()], maxval=len(tables)).start()
    table = tables.pop(0)
    for i in range(len(tables)):
        table = vstack([table] + [tables.pop(0)])
        pbar.update(i + 1)
    pbar.finish()
    return table


def from_list_to_array(d):
    """Transform lists (of dict of list) into numpy arrays."""
    if isinstance(d, (list, numpy.ndarray)):
        return numpy.array(d)
    for k in d:
        if isinstance(d[k], list):
            d[k] = numpy.array(d[k])
        elif isinstance(d[k], dict):
            from_list_to_array(d[k])
    return d


def stack_tables(data):
    """
    Stack the astropy tables across all patches.

    Return a new dictionnary of the form:
    d = {
    u: {'deepCoadd_forced_src': table, 'deepCoadd_meas': table},
    g: {'deepCoadd_forced_src': table, 'deepCoadd_meas': table}, ...
    }
    """
    print "Info: Stacking the data (patches, filters) into a single astropy table"
    return {'deepCoadd_meas': vstack([vstack([data[filt][patch]['deepCoadd_meas']
                                              for patch in data[filt]]) for filt in data]),
            'deepCoadd_forced_src': vstack([vstack([data[filt][patch]['deepCoadd_forced_src']
                                                    for patch in data[filt]]) for filt in data])}


def write_data(d, output, overwrite=False):
    """Write astropy 'deepCoadd_forced_src' and 'deepCoadd_meas' tables in an hdf5 file."""
    d['deepCoadd_forced_src'].write(output, path='deepCoadd_forced_src', compression=True,
                                    serialize_meta=True, overwrite=overwrite)
    d['deepCoadd_meas'].write(output, path='deepCoadd_meas', compression=True,
                              append=True, serialize_meta=True)
    save_wcs(d['wcs'], output)


def read_data(data_file, path=None):
    """Read astropy tables from an hdf5 file.

    :param string data_file: Name of the hdf5 file to load
    :param string path: Path (key) of the table to load
    :return: A dictionnary containing the following keys and values:

     - meas: the 'deepCoadd_meas' catalog (an astropy table)
     - forced: the 'deepCoad_forced_src' catalog (an astropy table)
     - wcs: the 'wcs' of these catalogs (an ``astropy.wcs.WCS`` object)
    """
    if path is None:
        try:
            return {'deepCoadd_meas': Table.read(data_file, path='deepCoadd_meas'),
                    'deepCoadd_forced_src': Table.read(data_file, path='deepCoadd_forced_src'),
                    'wcs': load_wcs(Table.read(data_file, path='wcs'))}
        except IOError:
            return Table.read(data_file)
    else:
        return Table.read(data_file, path=path)


def filter_table(t):
    """Apply a few quality filters on the data tables."""
    # Get the initial number of filter
    nfilt = len(t['deepCoadd_meas'].group_by('id').groups[0])

    # Select galaxies (and reject stars)
    filt = t['deepCoadd_meas']['base_ClassificationExtendedness_flag'] == 0  # keep galaxy
    filt &= t['deepCoadd_meas']['base_ClassificationExtendedness_value'] >= 0.5  # keep galaxy

    # Gauss regulerarization flag
    filt &= t['deepCoadd_meas']['ext_shapeHSM_HsmShapeRegauss_flag'] == 0

    # Make sure to keep primary sources
    filt &= t['deepCoadd_meas']['detect_isPrimary'] == 1

    # Check the flux value, which must be > 0
    filt &= t['deepCoadd_forced_src']['modelfit_CModel_flux'] > 0

    # Select sources which have a proper flux value
    filt &= t['deepCoadd_forced_src']['modelfit_CModel_flag'] == 0

    # Check the signal to noise (stn) value, which must be > 10
    filt &= (t['deepCoadd_forced_src']['modelfit_CModel_flux'] /
             t['deepCoadd_forced_src']['modelfit_CModel_fluxSigma']) > 10

    # Only keeps sources with the 5 filters
    dmg = t['deepCoadd_meas'][filt].group_by('id')
    dfg = t['deepCoadd_forced_src'][filt].group_by('objectId')

    # Indices difference is a quick way to get the lenght of each group
    filt = (dmg.groups.indices[1:] - dmg.groups.indices[:-1]) == nfilt

    return {'deepCoadd_meas': dmg.groups[filt],
            'deepCoadd_forced_src': dfg.groups[filt], 'wcs': t['wcs']}


def getdata(config, output='all_data.hdf5', output_filtered='filtered_data.hdf5', overwrite=False):
    """Shortcut function to get all the data from a bulter, fitler them, and save same."""
    if not overwrite:
        if os.path.exists(output) or os.path.exists(output_filtered):
            raise IOError("Output(s) already exist(s). Remove them or use overwrite=True.")
    if isinstance(config, str):
        config = load_config(config)
    d = get_all_data(config['butler'], config['patch'],
                     config['filter'], add_extra=True,
                     keys=config['keys'] if 'keys' in config else {})
    write_data(d, output, overwrite=overwrite)
    df = filter_table(d)
    write_data(df, output_filtered, overwrite=overwrite)
    return d, df


def correct_for_extinction(ti, te, mag='modelfit_CModel_mag', ext='sfd', ifilt="i_new"):
    """
    Compute extinction-corrected magnitude.

    :param table ti: input data table to fill with extinction-corrected magnitudes
    :param table te: input extinction table
    :param str mag: magnitude key from the catalog
    :param str ext: type of extinction map
    :param str ifilt: the 'i' filter you want to use (i_old or i_new)
    :return: a new column in the input table 'mag'+_extcorr.
    """
    # get the list of filter
    filters = list(set(te['filter']))

    # replace the 'i' filter by the one asked from the user
    for i, filt in enumerate(filters):
        if filt == 'i':
            filters[i] = ifilt

    # name of the new key
    magext = mag + '_extcorr'

    # Compute the corrected magnitude for each filter
    mcorr = numpy.zeros(len(ti[mag]))
    for f in filters:
        filt = ti['filter'] == (f if 'i' not in f else 'i')
        mcorr[filt] = ti[mag][filt] - te['albd_%s_%s' % (f, ext)][filt]

    # Add the new corrected-magnitudes column to the input data table
    ti.add_columns([Column(name=magext, data=mcorr, unit='mag',
                           description='Extinction corrected magnitude (i=%s, ext=%s)' %
                           (ifilt, ext))])


def filter_around(data, config, **kwargs):
    """Apply a circulat filter on the catalog around the center of the cluster.

    :param table data: The astropy table containing your data
    :param dict config: The analysis configuration file, which must contains the cluster coordinates

    List of available kwargs:

    :param float exclude_inner: Cut galaxies inside this radius [0]
    :param float exclude_outer: Cut galaxies outside this radius [inf]
    :param str unit: Unit of the input cuts [degree]
    :param bool plot: Do not produce a figure if set to False
    :return: A filter data table containing galaxie inside [exclude_inner, exclude_outer]
    """
    coord = SkyCoord(ra=[config['ra']], dec=[config['dec']], unit='deg')
    coords = SkyCoord(data['coord_ra_deg'], data['coord_dec_deg'], unit='deg')
    separation = coord.separation(coords)
    unit = kwargs.get('unit', 'degree')
    if hasattr(separation, unit):
        separation = getattr(separation, unit)
    else:
        arglist = "\n" + ", ".join(sorted([a for a in dir(separation) if not a.startswith('_')]))
        raise AttributeError("Angle instance has no attribute %s. Available attributes are: %s" %
                             (unit, arglist))
    data_around = data[(separation >= kwargs.get('exclude_inner', 0)) &
                       (separation < kwargs.get('exclude_outer', numpy.inf))]
    if kwargs.get('plot', True):
        title = "%s, %.2f < d < %.2f %s cut" % \
                (config['cluster'], kwargs.get('exclude_inner', 0),
                 kwargs.get('exclude_outer', numpy.inf), unit)
        plot_coordinates(data, data_around,
                         cluster_coord=(config['ra'], config['dec']), title=title)
    return data_around


def plot_coordinates(all_data, filtered_data, cluster_coord=None, title=None):
    """Plot a map of coordinates before and after a circular cut around the cluster center."""
    import pylab
    fig = pylab.figure()
    ax = fig.add_subplot(111, xlabel='ra', ylabel='dec')
    ax.scatter(all_data['coord_ra_deg'], all_data['coord_dec_deg'],
               color='k', label='All data', s=10)
    ax.scatter(filtered_data['coord_ra_deg'], filtered_data['coord_dec_deg'], color='b',
               label="Filtered data", s=10)
    if cluster_coord is not None:
        ax.scatter([cluster_coord[0]], [cluster_coord[1]], color='r', label='Cluster center',
                   marker='x', s=60)
    if title is not None:
        ax.set_title(title)
    ax.legend(loc='lower left', scatterpoints=1, frameon=False, fontsize='small')
    pylab.show()
