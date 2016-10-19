"""Data builder and parser for the Clusters package."""

import os
import yaml
import numpy
import h5py
from astropy.wcs import WCS, utils
from astropy.coordinates import SkyCoord, Angle
from astropy.table import Table, Column, vstack
from progressbar import Bar, ProgressBar, Percentage, ETA
from termcolor import colored


class Catalogs(object):

    """Load data from a LSST stack butler path"""

    def __init__(self, path):
        """."""
        # Load the bulter
        import lsst.daf.persistence as dafPersist
        print "INFO: Loading data from", path
        self.butler = dafPersist.Butler(path)

        # Initialize data dictionnaries
        self.dataids = {}
        self.catalogs = {}
        self.keys = {}
        self._getmag = self.wcs = None

    def _load_dataids(self, catalog, **kwargs):
        """Get the 'forced_src' catalogs."""
        print "INFO: Getting list of available data for", catalog
        if 'deepCoadd' in catalog:  # The deepCoadd* catalogs
            deepcoadd = [cat for cat in self.dataids if 'deepCoadd' in cat]
            if len(deepcoadd):
                dataids = self.dataids[deepcoadd[0]]
            else:
                filters = kwargs.get('filter', ['u', 'g', 'r', 'i', 'i2', 'z'])
                skymap = self.butler.get("deepCoadd_skyMap")
                dataids = [dict(tract=tract.getId(), patch="%d,%d" % patch.getIndex(), filter=f)
                           for tract in skymap for patch in tract for f in filters]
        else:  # The other catalogs
            keys = self.butler.getKeys(catalog)
            if 'tract' in keys:
                keys.pop('tract')
            dataids = [merge_dicts(dict(zip(keys, v)), {'tract': 0})
                       for v in self.butler.queryMetadata("forced_src", format=keys)]

        if len(dataids) == 0:
            raise IOError("No dataIds. Check the catalog, the config file, and path to the bulter.")

        # Specific selection make by the user?
        for kwarg in kwargs:
            if not kwarg in dataids[0]:
                continue
            print "INFO: Selecting data ids according to the '%s' selection" % kwarg
            print "  - input: %i data ids" % len(dataids)
            if not isinstance(kwargs[kwarg], list):
                kwargs[kwarg] = [kwargs[kwarg]]
            dataids = [dataid for dataid in dataids if dataid[kwarg] in kwargs[kwarg]]
            print "  - selected: %i data ids" % len(dataids)

        # Select the ccd/visit according to the input list of patch if given
        if 'deepCoadd' not in catalog:
            if 'patch' in kwargs and 'filter' in kwargs:
                print "INFO: Selecting visit/ccd according to the input list of patches"
                print "  - input: %i data ids" % len(dataids)
                dids = [{'filter': filt, 'patch': patch, 'tract': 0}
                        for filt in kwargs['filter'] for patch in kwargs['patch']
                        if self.butler.datasetExists('deepCoadd',
                                                     dataId={'filter': filt, 'patch':
                                                             patch, 'tract': 0})]
                coadds = [self.butler.get('deepCoadd', did) for did in dids]
                ccds = [coadd.getInfo().getCoaddInputs().ccds for coadd in coadds]
                ccds_visits = [numpy.transpose([ccd.get('visit'), ccd.get('ccd')]) for ccd in ccds]
                ccds_visits = [list(c) for c in numpy.concatenate(ccds_visits)]
                dataids = [dataid for dataid in dataids
                           if [dataid['visit'], dataid['ccd']] in ccds_visits]
                print "  - selected: %i data ids" % len(dataids)

        # Only keep dataids with data
        self.dataids[catalog] = [dataid for dataid in dataids if
                                 self.butler.datasetExists(catalog, dataId=dataid)]
        print "INFO: %i data ids finally kept" % len(self.dataids[catalog])


    def _load_catalog_dataid(self, catalog, dataid, astropy_table=True, **kwargs):
        """Load a catalog from a 'dataId' set of parameter."""
        cat = self.butler.get(catalog, dataId=dataid)
        return cat if not astropy_table else get_astropy_table(cat, **kwargs)

    def _load_catalog(self, catalog, **kwargs):
        """Load a given catalog."""
        self._load_dataids(catalog, **kwargs)
        print "INFO: Getting the data from the butler for %i fits files" % \
            len(self.dataids[catalog])
        pbar = progressbar(len(self.dataids[catalog]))
        tables = [self._get_tables(catalog, did, i, pbar)
                  for i, did in enumerate(self.dataids[catalog])]
        pbar.finish()
        print "INFO: Stacking the %i tables together" % len(tables)
        self.catalogs[catalog] = vstack2(tables)

    def _get_tables(self, catalog, did, i, pbar):
        """Get a table and add a few keys."""
        table = self._load_catalog_dataid(catalog, did, **{'keys': self.keys[catalog]})
        for key in did:
            table.add_column(Column(name=key, data=[did[key]] * len(table)))
        pbar.update(i + 1)
        return table

    def _match_ids(self):
        """Select in the 'forced_src' catalog the source that are in the deepCoad catalogs"""
        deepcoadd = [cat for cat in self.catalogs if 'deepCoadd' in cat]
        if len(deepcoadd):
            if 'forced_src' in self.catalogs:
                print colored("INFO: Matching 'forced_src' and 'deepCoadd' catalogs", "green")
                print "  - Initially %i sources in the forced-src catalog" % \
                    len(self.catalogs['forced_src'])
                filt = numpy.where(numpy.in1d(self.catalogs['forced_src']['objectId'],
                                              self.catalogs[deepcoadd[0]]['id']))[0]
                self.catalogs['forced_src'] = self.catalogs['forced_src'][filt]
                print "  - %i sources in the forced-src catalog after selection\n" % \
                    len(self.catalogs['forced_src'])
            else:
                print "WARNING: forced_src catalogs not loaded. No match possible."
        else:
            print "WARNING: No deepCoadd* catalog loaded. No match possible."

    def _add_new_columns(self):
        """Compute magns for all fluxes of a given table. Add the corresponding new columns.
        Compute the x/y position in pixel for all sources. Add new columns to the table."""
        print "INFO: Adding magnitude and coordinates columns to all catalog tables"
        for catalog in self.catalogs:
            print "  - for", catalog
            # skip wcs key
            if catalog == 'wcs':
                continue
            # Shortuc to the table
            table = self.catalogs[catalog]

            # Add magnitudes
            kfluxes = [k for k in table.columns if k.endswith('_flux')]
            ksigmas = [k + 'Sigma' for k in kfluxes]
            print "    -> getting magnitudes"
            for kflux, ksigma in zip(kfluxes, ksigmas):
                if kflux.replace('_flux', '_mag') in table.keys():
                    continue
                mag, dmag = numpy.array([self._mag(f, s)
                                         for f, s in zip(table[kflux], table[ksigma])]).T
                table.add_columns([Column(name=kflux.replace('_flux', '_mag'), data=mag,
                                          description='Magnitude', unit='mag'),
                                   Column(name=ksigma.replace('_fluxSigma', '_magSigma'), data=dmag,
                                          description='Magnitude error', unit='mag')])

            if 'x_Src' in table.keys():
                return
            # Add the x / y position in pixel
            print "    -> getting pixel coordinates"
            xsrc, ysrc = skycoord_to_pixel(SkyCoord(table["coord_ra"].tolist(),
                                                    table["coord_dec"].tolist(), unit='rad'),
                                           self.wcs)
            table.add_columns([Column(name='x_Src', data=xsrc,
                                      description='x coordinate', unit='pixel'),
                               Column(name='y_Src', data=ysrc,
                                      description='y coordinate', unit='pixel')])

            # Add a new column to have to coordinates in degree
            print "    -> getting degree coordinates"
            table.add_columns([Column(name='coord_ra_deg',
                                      data=Angle(table['coord_ra'].tolist(), unit='rad').degree,
                                      description='RA coordinate', unit='degree'),
                               Column(name='coord_dec_deg',
                                      data=Angle(table['coord_dec'].tolist(), unit='rad').degree,
                                      description='DEC coordinate', unit='degree')])

    def _load_calexp(self, **kwargs):
        """Load the 'calexp' info in order to get the WCS and the magnitudes."""
        print colored("\nINFO: Loading the 'calexp' info", 'green')
        calcat = 'deepCoadd_calexp'
        self._load_dataids(calcat, **kwargs)
        print "INFO: Getting the %s catalog for one dataId" % calcat
        calexp = self._load_catalog_dataid(calcat, self.dataids[calcat][0], astropy_table=False)
        print "INFO: Getting the magnitude function"
        self._getmag = calexp.getCalib().getMagnitude
        print "INFO: Getting the wcs function"
        wcs = calexp.getWcs().getFitsMetadata().toDict()
        self.wcs = WCS(wcs)
        self.catalogs['wcs'] = Table({k: [wcs[k]] for k in wcs})

    def _mag(self, flux, sigma):
        """Redefine the magnitude function. Negative flux or sigma accepted."""
        return (numpy.nan, numpy.nan) if (flux <= 0 or sigma <= 0) else self._getmag(flux, sigma)

    def load_catalogs(self, catalogs, keys=None, **kwargs):
        """Load a list of catalogs, e.g.,

        :param str/list catalogs: A catalog name, or a list of catalogs (see below)
        :param dict keys: A dictionnary of keys to load for each catalog

        Available kwargs are:

        :param bool update: Set to True if you want to update an already loaded catalog
        :param bool show: Set to True to get all available keys of a (list of) catalog(s)
        :param bool matchid: Will only keep objects which are in the deepCoad catalogs (to be used
                             when loading the forced_src and deepCoadd catalogs)

        Examples of catalogs that you can load:

         - 'deepCoadd_meas',
         - 'deepCoadd_forced_src',
         - 'deepCoadd_calexp',
         - 'forced_src'
        """
        if 'show' in kwargs:
            self.show_keys(catalogs)
            return
        keys = {} if keys is None else keys
        catalogs = [catalogs] if isinstance(catalogs, str) else catalogs
        for catalog in catalogs:
            if catalog in self.catalogs and 'update' not in kwargs:
                print "WARNING: %s is already loaded. Use 'update' to reload it." % catalog
                continue
            if 'calexp' in catalog:
                print "WARNING: Skipping %s. This is not a regular catalog (no schema).\n" % catalog
                continue
            print colored("INFO: Loading the %s catalog" % catalog, 'green')
            self.keys[catalog] = keys.get(catalog, "*")
            self._load_catalog(catalog, **kwargs)
            print "INFO: %s loaded." % catalog
        self._filter_around(**kwargs)
        if 'matchid' in kwargs:
            self._match_ids()
        self._load_calexp(**kwargs)
        self._add_new_columns()
        print "\nINFO: Done loading the data."

    def _filter_around(self, **kwargs):
        """Filter out galaxies outside a given radius."""
        if not 'radius' in kwargs:
            return
        print colored("\nINFO: Select galaxies inside a radius of %s around cluster center" % \
                      kwargs['radius'], "green")
        radius = float(kwargs['radius'].split()[0])
        unit = kwargs['radius'].split()[1]
        for catalog in self.catalogs:
            print "  -", catalog
            self.catalogs[catalog] = filter_around(self.catalogs[catalog], kwargs,
                                                   exclude_outer=radius, unit=unit)

    def show_keys(self, catalogs=None):
        """Show all the available keys."""
        if catalogs is None:
            catalogs = self.catalogs.keys()
        catalogs = [catalogs] if isinstance(catalogs, str) else catalogs
        if len(catalogs) == 0:
            print "WARNING: No catalog loaded nor given."
            return
        for cat in catalogs:
            print "INFO: Available list of keys for the %s catalog" % cat
            if cat not in self.dataids:
                self._load_dataids(cat)
            table = get_astropy_table(self.butler.get(cat, dataId=self.dataids[cat][0]), keys="*")
            ktable = Table(numpy.transpose([[k, table[k].description, table[k].unit]
                                            for k in sorted(table.keys())]).tolist(),
                           names=["Keys", "Description", "Units"])
            print "  -> All saved in %s_keys.txt" % cat
            ktable.write("%s_keys.txt" % cat, format='ascii', comment="#")

    def save_catalogs(self, output_name, overwrite=False):
        """Save the catalogs into an hdf5 file."""
        if not output_name.endswith('.hdf5'):
            output_name += '.hdf5'
        print "INFO: Saving the catalogs in", output_name
        for i, cat in enumerate(self.catalogs):
            print "  - saving", cat
            if i == 0:
                self.catalogs[cat].write(output_name, path=cat, compression=True,
                                         serialize_meta=True, overwrite=overwrite)
            else:
                self.catalogs[cat].write(output_name, path=cat, compression=True,
                                         serialize_meta=True, append=True)
        print "INFO: Saving done."


def progressbar(maxnumber, prefix='loading'):
    """Create and return a standard progress bar."""
    return ProgressBar(widgets=['  - %s ' % prefix, Percentage(), Bar(marker='>'), ETA()],
                       term_width=60, maxval=maxnumber).start()


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


def add_magnitudes(t, getmagnitude):
    """Compute magns for all fluxes of a given table. Add the corresponding new columns."""
    kfluxes = [k for k in t.columns if k.endswith('_flux')]
    ksigmas = [k + 'Sigma' for k in kfluxes]
    for kf, ks in zip(kfluxes, ksigmas):
        m, dm = numpy.array([getmagnitude(f, s) for f, s in zip(t[kf], t[ks])]).T
        t.add_columns([Column(name=kf.replace('_flux', '_mag'), data=m,
                              description='Magnitude', unit='mag'),
                       Column(name=ks.replace('_fluxSigma', '_magSigma'), data=dm,
                              description='Magnitude error', unit='mag')])


def add_position_and_deg(t, wcs):
    """Compute the x/y position in pixel for all sources. Add new columns to the table."""
    # Add the x / y position in pixel
    x, y = skycoord_to_pixel(SkyCoord(t["coord_ra"].tolist(),
                                      t["coord_dec"].tolist(), unit='rad'), wcs)
    t.add_columns([Column(name='x_Src', data=x,
                          description='x coordinate', unit='pixel'),
                   Column(name='y_Src', data=y,
                          description='y coordinate', unit='pixel')])

    # Add a new column to have to coordinates in degree
    t.add_columns([Column(name='coord_ra_deg',
                          data=Angle(t['coord_ra'].tolist(), unit='rad').degree,
                          description='RA coordinate', unit='degree'),
                   Column(name='coord_dec_deg',
                          data=Angle(t['coord_dec'].tolist(), unit='rad').degree,
                          description='DEC coordinate', unit='degree')])


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


def pixel_to_skycoord(x, y, wcs):
    """Transform pixel coordinates (x, y) to sky coordinates (ra, dec in deg) given a wcs.

    :param float x: x coordinate
    :param float y: y coordinate
    :param wcs: an astropy.wcs.WCS object
    :return: an astropy.coordinates.SkyCoord object.
    """
    return utils.pixel_to_skycoord(x, y, wcs)


def merge_dicts(*dict_args):
    """Merge two dictionnary.

    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def vstack2(tables):
    """Verticaly stack large amount of astropy tables."""
    pbar = progressbar(len(tables), prefix='stacking')
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


def write_data(d, output, overwrite=False):
    """Write astropy 'deepCoadd_forced_src' and 'deepCoadd_meas' tables in an hdf5 file."""
    d['deepCoadd_forced_src'].write(output, path='deepCoadd_forced_src', compression=True,
                                    serialize_meta=True, overwrite=overwrite)
    d['deepCoadd_meas'].write(output, path='deepCoadd_meas', compression=True,
                              append=True, serialize_meta=True)
    save_wcs(d['wcs'], output)


def read_hdf5(hdf5_file, path=None):
    """Read astropy tables from an hdf5 file.

    :param string data_file: Name of the hdf5 file to load.
    :param string path: Path (key) of the table to load.
    :return: A dictionnary containing the table(s). Keys are the path names.
    """
    if path is None:
        paths = hdf5_paths(hdf5_file)
        return {path: Table.read(hdf5_file, path=path) for path in paths}
    else:
        return {path: Table.read(hdf5_file, path=path)}


def hdf5_paths(hdf5_file):
    """Get all the available paths of an hdf5 file.

    :param str hdf5_file: Name of the hdf5 file to load
    :return: The available list of paths in the input hdf5 file
    """
    hdf5_content = h5py.File(hdf5_file, 'r')
    paths = hdf5_content.keys()
    hdf5_content.close()
    return paths


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


#def getdata(config, output='all_data.hdf5', output_filtered='filtered_data.hdf5', overwrite=False):
#    """Shortcut function to get all the data from a bulter, fitler them, and save same."""
#    if not overwrite:
#        if os.path.exists(output) or os.path.exists(output_filtered):
#            raise IOError("Output(s) already exist(s). Remove them or use overwrite=True.")
"    if isinstance(config, str):
#        config = load_config(config)
#    
#    d = get_all_data(config['butler'], config['patch'],
#                     config['filter'], add_extra=True,
#                     keys=config['keys'] if 'keys' in config else {})
#    write_data(d, output, overwrite=overwrite)
#    df = filter_table(d)
#    write_data(df, output_filtered, overwrite=overwrite)
#    return d, df


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
    :param bool plot: Produce a figure if set to False
    :return: A filter data table containing galaxie inside [exclude_inner, exclude_outer]
    """
    coord = SkyCoord(ra=[config['ra']], dec=[config['dec']], unit='deg')
    datag = data.group_by('filter')
    same_length = len(set([len(g) for g in datag.groups])) == 1
    if same_length:
        coords = SkyCoord(ra=list(datag.groups[0]['coord_ra']),
                          dec=list(datag.groups[0]['coord_dec']), unit='rad')
    else:
        coords = SkyCoord(ra=list(numpy.round(data['coord_ra'], 3)),
                          dec=list(numpy.round(data['coord_dec'], 3)), unit='rad')
    separation = coord.separation(coords)
    unit = kwargs.get('unit', 'degree')
    if hasattr(separation, unit):
        separation = getattr(separation, unit)
    else:
        arglist = "\n" + ", ".join(sorted([a for a in dir(separation) if not a.startswith('_')]))
        raise AttributeError("Angle instance has no attribute %s. Available attributes are: %s" %
                             (unit, arglist))
    filt = (separation >= kwargs.get('exclude_inner', 0)) & \
           (separation < kwargs.get('exclude_outer', numpy.inf))
    data_around = vstack([group[filt] for group in datag.groups]) if same_length else data[filt]
    if kwargs.get('plot', False):
        title = "%s, %.2f < d < %.2f %s cut" % \
                (config['cluster'], kwargs.get('exclude_inner', 0),
                 kwargs.get('exclude_outer', numpy.inf), unit)
        plot_coordinates(data, data_around,
                         cluster_coord=(config['ra'], config['dec']), title=title)
    if 'pbar' in kwargs:
        kwargs['pbar'].update(kwargs['i'] + 1)
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


def plot_patches(catalog, clust_coords=None):
    import pylab
    colors = pylab.cm.jet(pylab.linspace(0, 1, len(set(catalog['patch']))))
    fig = pylab.figure()
    ax = fig.add_subplot(111, xlabel="RA", ylabel="DEC")
    for i, path in enumerate(sorted(set(catalog['patch']))):
        filt = catalog['patch'] == path
        ra = catalog['coord_ra_deg'][filt]
        dec = catalog['coord_dec_deg'][filt]
        ax.scatter(ra, dec, color=colors[i], label=path, s=10)
        ax.vlines(min(ra), min(dec), max(dec), color='k')
        ax.vlines(max(ra), min(dec), max(dec), color='k')
        ax.hlines(min(dec), min(ra), max(ra), color='k')
        ax.hlines(max(dec), min(ra), max(ra), color='k')
        ax.legend(loc='best', numpoints=1, frameon=False)
    if clust_coords is not None:
        ax.scatter(clust_coords['ra'], clust_coords['dec'], s=100, marker='s', color='k')
    pylab.show()
