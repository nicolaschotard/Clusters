"""Some utilities, mostly related to the data."""


from __future__ import print_function
import os
import numpy as np
import h5py
from astropy.wcs import WCS, utils
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column, vstack
from astropy.units import Quantity
from progressbar import Bar, ProgressBar, Percentage, ETA
import yaml


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
    c = yaml.load(open(config))

    # if the user did not provide a 'zphot' key or 'mass' key option,
    # makes sure a default name is setup in the config dictionnary
#    if 'zphot' not in c.keys() : c['zphot'] = {'zphot_ref':{}}
#    if 'mass' not in c.keys() : c['mass'] = {'zconfig':'zphot_ref'}

    return c


def load_wcs(wcs):
    """Get back the right wcs format from the hdf5 table."""
    return WCS({k: wcs[k].item() if not isinstance(wcs[k].item(), bytes)
                else wcs[k].item().decode()
                for k in wcs.keys()})


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

    :param float xsrc: x coordinate of the source
    :param float ysrc: y coordinate of the source
    :param wcs: an astropy.wcs.WCS object
    :return: an astropy.coordinates.SkyCoord object.
    """
    return utils.pixel_to_skycoord(xsrc, ysrc, wcs)


def read_hdf5(hdf5_file, path=None, dic=True):
    """Read astropy tables from an hdf5 file.

    :param string data_file: Name of the hdf5 file to load.
    :param string/list path: Path(s) (key) of the table to load.
    :return: A dictionnary containing the table(s). Keys are the path names.
    """
    if path is None:
        paths = hdf5_paths(hdf5_file)
        return {path: Table.read(hdf5_file, path=path) for path in paths}
    else:
        if isinstance(path, list):
            return {p: Table.read(hdf5_file, path=p) for p in path} if dic \
                else (Table.read(hdf5_file, path=p) for p in path)
        else:
            return {path: Table.read(hdf5_file, path=path)} if dic \
                else Table.read(hdf5_file, path=path)


def merge_dicts(*dict_args):
    """Merge two dictionnary.

    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def concatenate_dicts(*dicts):
    """Concatenate dictionnaries containing numpy arrays."""
    return {k: np.concatenate([d.pop(k) for d in dicts]) for k in dicts[0].keys()}


def read_hdf5(hdf5_file, path=None, dic=True):
    """Read astropy tables from an hdf5 file.
    :param string data_file: Name of the hdf5 file to load.
    :param string/list path: Path(s) (key) of the table to load.
    :return: A dictionnary containing the table(s). Keys are the path names.
    """
    if path is None:
        paths = hdf5_paths(hdf5_file)
        return {path: Table.read(hdf5_file, path=path) for path in paths}
    else:
        if isinstance(path, list):
            return {p: Table.read(hdf5_file, path=p) for p in path} if dic \
                else (Table.read(hdf5_file, path=p) for p in path)
        else:
            return {path: Table.read(hdf5_file, path=path)} if dic \
                else Table.read(hdf5_file, path=path)

            
def hdf5_paths(hdf5_file):
    """Get all the available paths of an hdf5 file.

    :param str hdf5_file: Name of the hdf5 file to load
    :return: The available list of paths in the input hdf5 file
    """
    hdf5_content = h5py.File(hdf5_file, 'r')
    paths = list(hdf5_content.keys())
    hdf5_content.close()
    return paths


def filter_table(cats):
    """Apply a few quality filters on the data tables."""
    # == Get the initial number of filter
    nfilt = len(set(cats['deepCoadd_meas']['filter']))

    # == Filter the deepCoadd catalogs

    # Select galaxies (and reject stars)
    # keep galaxy
    filt = cats['deepCoadd_meas']['base_ClassificationExtendedness_flag'] == 0

    # keep galaxy
    filt &= cats['deepCoadd_meas']['base_ClassificationExtendedness_value'] >= 0.5

    # Gauss regulerarization flag
    filt &= cats['deepCoadd_meas']['ext_shapeHSM_HsmShapeRegauss_flag'] == 0

    # Make sure to keep primary sources
    filt &= cats['deepCoadd_meas']['detect_isPrimary'] == 1

    # Check the flux value, which must be > 0
    filt &= cats['deepCoadd_forced_src']['modelfit_CModel_flux'] > 0

    # Select sources which have a proper flux value
    filt &= cats['deepCoadd_forced_src']['modelfit_CModel_flag'] == 0

    # Check the signal to noise (stn) value, which must be > 10
    filt &= (cats['deepCoadd_forced_src']['modelfit_CModel_flux'] /
             cats['deepCoadd_forced_src']['modelfit_CModel_fluxSigma']) > 10

    # == Only keeps sources with the 'nfilt' filters
    dmg = cats['deepCoadd_meas'][filt].group_by('id')
    dfg = cats['deepCoadd_forced_src'][filt].group_by(
        'id' if 'id' in cats['deepCoadd_forced_src'].keys() else 'objectId')

    # Indices difference is a quick way to get the lenght of each group
    filt = (dmg.groups.indices[1:] - dmg.groups.indices[:-1]) == nfilt

    output = {'deepCoadd_meas': dmg.groups[filt],
              'deepCoadd_forced_src': dfg.groups[filt], 'wcs': cats['wcs']}

    # == Filter the forced_src catalog: only keep objects present in the other catalogs
    if "forced_src" not in cats.keys():
        return output

    filt = np.where(np.in1d(cats['forced_src']['objectId'],
                            output['deepCoadd_meas']['id']))[0]
    output['forced_src'] = cats['forced_src'][filt]

    return output


def correct_for_extinction(data, extinction, mag='modelfit_CModel_mag', ext='sfd', ifilt="i_new"):
    """
    Compute extinction-corrected magnitude.

    :param table data: input data table to fill with extinction-corrected magnitudes, which also
                       contains the extinction values
    :param table extinction: input extinction table
    :param str mag: magnitude key from the catalog
    :param str ext: type of extinction map
    :param str ifilt: the 'i' filter you want to use (i_old or i_new)
    :return: a new column in the input table 'mag'+_extcorr.
    """
    # get the list of filter
    filters = list(set(extinction['filter']))

    # replace the 'i' filter by the one asked from the user
    for i, filt in enumerate(filters):
        if filt == 'i':
            filters[i] = ifilt

    # name of the new key
    magext = mag + '_extcorr'

    # available dust maps
    dustmaps = set([k.split('_')[-1]
                    for k in extinction.keys() if k.startswith('albd_')])
    if ext not in dustmaps:
        raise IOError("ERROR: The selected dustmap (%s) must be in" %
                      ext, dustmaps)

    # Compute the corrected magnitude for each filter
    mcorr = np.zeros(len(data[mag]))
    for f in filters:
        filt = data['filter'] == (f if 'i' not in f else 'i')
        mcorr[filt] = data[mag][filt] - \
            extinction['albd_%s_%s' % (f, ext)][filt]

    # Add the new corrected-magnitudes column to the input data table
    data.add_columns([Column(name=magext, data=mcorr, unit='mag',
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
    :param plot: Produce a figure if specified
    :return: A filter data table containing galaxie inside [exclude_inner, exclude_outer]
    """
    plot = kwargs.get('plot', True)

    datag = data.group_by('filter')
    same_length = len(set([len(g) for g in datag.groups])) == 1
    if same_length:
        ra = Quantity(datag.groups[0]['coord_ra'].tolist(), 'rad')
        dec = Quantity(datag.groups[0]['coord_dec'].tolist(), 'rad')
    else:
        ra = Quantity(np.round(data['coord_ra'].tolist(), 3), 'rad')
        dec = Quantity(np.round(data['coord_dec'].tolist(), 3), 'rad')
    sep = SkyCoord(Quantity([config['ra']], 'deg'),
                   Quantity([config['dec']], 'deg')).separation(SkyCoord(ra, dec))
    if hasattr(sep, kwargs.get('unit', 'degree')):
        sep = getattr(sep, kwargs.get('unit', 'degree'))
    else:
        raise AttributeError("Angle instance has no attribute %s. Available attributes are: %s" %
                             (kwargs.get('unit', 'degree'),
                              "\n" + ", ".join(sorted([a for a in dir(sep)
                                                       if not a.startswith('_')]))))
    filt = (sep >= kwargs.get('exclude_inner', 0)) & \
           (sep < kwargs.get('exclude_outer', np.inf))
    data_around = vstack(
        [group[filt] for group in datag.groups]) if same_length else data[filt]
    if plot:
        title = "%s, %.2f < d < %.2f %s cut" % \
                (config['cluster'], kwargs.get('exclude_inner', 0),
                 kwargs.get('exclude_outer', np.inf), kwargs.get('unit', 'degree'))
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
    ax.legend(loc='lower left', scatterpoints=1,
              frameon=False, fontsize='small')
    pylab.show()


def plot_patches(catalog, clust_coords=None):
    """Plot patches in the RA/DEC parameter space."""
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
        ax.scatter(clust_coords['ra'], clust_coords['dec'],
                   s=100, marker='s', color='k')
    pylab.show()


def overwrite_or_append(filename, path, table, overwrite=False):
    """
    Overwrites or append new path/table to existing file or creates new file

    The overwrite keyword of data.write(file,path) does not overwrites
    only the data in path, but the whole file, i.e. (we lose all other
    paths in the process) --> need to do it by hand
    """

    if not os.path.isfile(filename):
        print("Creating", filename)
        table.write(filename, path=path, compression=True, serialize_meta=True)
    else:
        data = read_hdf5(filename)
        if path in data.keys():
            if overwrite:
                print("Overwriting path =", path, " in", filename)
                data[path] = table  # update the table with new values
                os.remove(filename)  # delete file
                for p in data.keys():  # rewrite all paths/tables to file
                    data[p].write(filename, path=p, compression=True, serialize_meta=True,
                                  append=True)
            else:
                raise IOError(
                    "Path already exists in hdf5 file. Use --overwrite to overwrite.")
        else:
            print("Adding", path, " to", filename)
            table.write(filename, path=path, compression=True, serialize_meta=True,
                        append=True)
