"""Shear analysis."""

import numpy
import pylab
import seaborn
from astropy.table import Table, Column
from . import data as cdata

def compute_shear(e1, e2, distx, disty):
    """Compute the shear."""
    phi = numpy.arctan2(disty, distx)
    gamt = - (e1 * numpy.cos(2.0 * phi) + e2 * numpy.cos(2.0 * phi))
    gamc = - e1 * numpy.sin(2.0 * phi) + e2 * numpy.cos(2.0 * phi)
    dist = numpy.sqrt(distx**2 + disty**2)
    return gamt, gamc, dist


def analysis(table, xclust, yclust):
    """Computethe shear.

    :param string data_file: Name of the hdf5 file to load
    :param string path: Path (key) of the table to load
    :return: A dictionnary containing the following keys and values:

     - meas: the 'deepCoadd_meas' catalog (an astropy table)
     - forced: the 'deepCoad_forced_src' catalog (an astropy table)
     - wcs: the 'wcs' of these catalogs (an ``astropy.wcs.WCS`` object)
    """
    e1r = table["ext_shapeHSM_HsmShapeRegauss_e1"][table['filter'] == 'r']
    e2r = table["ext_shapeHSM_HsmShapeRegauss_e2"][table['filter'] == 'r']
    e1i = table["ext_shapeHSM_HsmShapeRegauss_e1"][table['filter'] == 'i']
    e2i = table["ext_shapeHSM_HsmShapeRegauss_e2"][table['filter'] == 'i']
    distx = table["x_Src"][table['filter'] == 'r'] - xclust
    disty = table["y_Src"][table['filter'] == 'r'] - yclust

    # Quality cuts
    # magnitude cut
    filt = table['modelfit_CModel_mag'][table['filter'] == 'r'] < 23.5
    # resolution cut
    filt &= table['ext_shapeHSM_HsmShapeRegauss_resolution'][table['filter'] == 'i'] > 0.4
    # ellipticity cut
    filt &= (abs(e1i) < 1) & (abs(e2i) < 1)
    # er ~= ei
    filt &= (abs(e1r - e1i) < 0.5) & (abs(e2r - e2i) < 0.5)

    # Apply cuts
    e1i, e2i, distx, disty = [x[filt] for x in [e1i, e2i, distx, disty]]

    # Comput the shear
    gamt, gamc, dist = compute_shear(e1i, e2i, distx, disty)

    # Make some plots
    plot_shear(gamt, gamc, dist)


def xy_clust(config, wcs):
    return cdata.skycoord_to_pixel([config['ra'], config['dec']], wcs)
        
def compare_shear(catalogs, xclust, yclust, qcut=None):
    """Compare shear mesured on the coadd and shear measured on indivial ccd.

    For now, do:
    from Clusters import data2
    from Clusters import shear
    config = data2.load_config('MACSJ2243.3-0935.yaml')
    catalogs = data2.read_hdf5('test_data2.hdf5')
    xc, yc = shear.xy_clust(config, data2.load_wcs(catalogs['wcs']))
    tables = shear.compare_shear([catalogs['deepCoadd_meas'], catalogs['forced_src']], xc, yc)
    """

    # Compute shear and distance for all srouces in both catalogs
    # And add that info into the tables
    tables = []
    print "INFO: %i catalogs to load" % len(catalogs)
    for i, cat in enumerate(catalogs):
        filtr = cat['filter'] == 'r'
        filti = (cat['filter'] == 'i') | (cat['filter'] == 'i2')
        if 'objectId' in cat.keys():
            objectids = cat["objectId"][filti]
        else:
            objectids = cat["id"][filti]
        e1i = cat["ext_shapeHSM_HsmShapeRegauss_e1"][filti]
        e2i = cat["ext_shapeHSM_HsmShapeRegauss_e2"][filti]
        e1r = cat["ext_shapeHSM_HsmShapeRegauss_e1"][filtr]
        e2r = cat["ext_shapeHSM_HsmShapeRegauss_e2"][filtr]
        distx = cat["x_Src"][filti] - xclust
        disty = cat["y_Src"][filti] - yclust

        # Quality cuts
        # resolution cut
        filt = cat['ext_shapeHSM_HsmShapeRegauss_resolution'][filti] > 0.3

        # ellipticity cut
        filt &= (abs(e1i) < 1) & (abs(e2i) < 1)

        if qcut is not None and qcut[i] is True:
            # Select galaxies (and reject stars)
            filt &= cat['base_ClassificationExtendedness_flag'][filti] == 0  # keep galaxy
            filt &= cat['base_ClassificationExtendedness_value'][filti] >= 0.5  # keep galaxy

            # Gauss regulerarization flag
            filt &= cat['ext_shapeHSM_HsmShapeRegauss_flag'][filti] == 0

            # Make sure to keep primary sources
            filt &= cat['detect_isPrimary'][filti] == 1

            # magnitude cut
            filt &= cat['modelfit_CModel_mag'][filtr] < 23.5

            # er ~= ei
            filt &= (abs(e1r - e1i) < 0.5) & (abs(e2r - e2i) < 0.5)

        # Apply cuts
        e1i, e2i, distx, disty, objectids = [x[filt] for x in [e1i, e2i, distx,
                                                               disty, objectids]]

        tshear, cshear, dist = compute_shear(e1i, e2i, distx, disty)

        tables.append(Table([Column(name='Tshear', data=tshear, description='Tangential shear'),
                             Column(name='Cshear', data=cshear, description='Cross shear'),
                             Column(name='Distance', data=dist, description='Distance to center'),
                             Column(name='objectId', data=objectids, description='Object ID')]))
    print "INFO: Done loading shear data"

    ids = tables[0]['objectId'][numpy.argsort(tables[0]['objectId'])].tolist()
    tshear_coadd = numpy.array(tables[0]['Tshear'][numpy.argsort(tables[0]['objectId'])].tolist())
    tshear_ccds = [tables[1]['Tshear'][tables[1]['objectId'] == oid].tolist() for oid in ids]
    
    fig = pylab.figure()
    ax = fig.add_subplot(111, xlabel='Distance', ylabel='T-shear')
    ax.scatter(tables[1]['Distance'], tables[1]['Tshear'], color='k')
    ax.scatter(tables[0]['Distance'], tables[0]['Tshear'], color='r')
    ax.set_title("%i sources" % len(tables[0]['Distance']))

    fig = pylab.figure()
    ax = fig.add_subplot(111, xlabel='T-shear (coadd)', ylabel='T-shear (ccd)')
    for i, tshear_ccd in enumerate(tshear_ccds):
        ax.scatter([tshear_coadd[i]] * len(tshear_ccd), tshear_ccd, color='k')
    means = numpy.array([numpy.mean(tshear_ccd) for tshear_ccd in tshear_ccds])
    stds = numpy.array([numpy.std(tshear_ccd) for tshear_ccd in tshear_ccds])
    ax.scatter(tshear_coadd, means, color='r')
    ax.errorbar(tshear_coadd, means, yerr=stds, color='r', ls='None', capsize=20)
    ax.plot([numpy.min(tshear_coadd), numpy.max(tshear_coadd)],
            [numpy.min(tshear_coadd), numpy.max(tshear_coadd)])
    filt = numpy.isfinite(means) & numpy.isfinite(stds) & numpy.isfinite(tshear_coadd)
    stds = numpy.array([(s if s != 0. else numpy.median(stds[filt])) for s in stds])
    chi2 = sum((means[filt] - tshear_coadd[filt])**2 / stds[filt]**2)
    print chi2 / sum(filt)
    ax.set_title("%i sources" % len(tables[0]['Distance']))
    pylab.show()
    return tables


def plot_shear(gamt, gamc, dist, drange=(0, 8500), nbins=8):
    """Plot shear."""
    dval, step = numpy.linspace(drange[0], drange[1], nbins, retstep=True)

    plot_hist([gamt, gamc], ['Gamt', 'Gamc'])
    plot_hist([gamt, gamc], ['Gamt', 'Gamc'], nbins=80, xarange=(-0.2, 0.2))

    plot_scatter([dist, dist], [gamt, gamc],
                 ['Dist', 'Dist'], ['Gamt', 'Gamc'], yarange=(-1, 1))

    masks = [(dist > d - step / 2) & (dist <= d + step / 2) for d in dval]
    tshear = [numpy.mean(gamt[mask]) for mask in masks]
    cshear = [numpy.mean(gamc[mask]) for mask in masks]
    tsheare = [numpy.std(gamt[mask]) / numpy.sqrt(len(gamt[mask])) for mask in masks]
    csheare = [numpy.std(gamc[mask]) / numpy.sqrt(len(gamc[mask])) for mask in masks]

    plot_scatter([dval, dval], [tshear, cshear],
                 ['Distance to cluster center (px)', 'Distance to cluster center (px)'],
                 ['Tangential shear', 'Cross shear'], yerrs=[tsheare, csheare],
                 xarange=(-500, 9000), yarange=(-0.06, 0.08))

    pylab.show()


def plot_hist(xs, labels, nbins=200, xarange=(-2, 2)):
    """Plot multiple histograms in subplots."""
    fig = pylab.figure(figsize=(15, 8))
    for i, x in enumerate(xs):
        ax = fig.add_subplot(1, len(xs), i + 1, xlabel=labels[i])
        ax.hist(x, bins=nbins, range=xarange)


def plot_scatter(xs, ys, xlabels, ylabels, **kwargs):
    """Plot multiple scatter plots in subplots.

    :param list xs: List of arrays for x axis
    :param list ys: List of arrays for y axis
    :param str xlabels: List of x labels
    :param str ylabels: List of y labels

    List of available kwargs:
    :param list yerrs: List of arrays, error on the y axis
    :param list xarange: Range for x axis (min,max)
    :param list yarange: Range for y axis (min,max)
    """
    fig = pylab.figure(figsize=(15, 8))
    for i, x in enumerate(xs):
        ax = fig.add_subplot(1, len(xs), i + 1, xlabel=xlabels[i], ylabel=ylabels[i])
        ax.axhline(0, color='k', ls=':')
        ax.scatter(x, ys[i], s=1, color='b')
        if 'yerrs' in kwargs:
            ax.errorbar(x, ys[i], yerr=kwargs['yerrs'][i])
        if 'xarange' in kwargs:
            ax.set_xlim(kwargs['xarange'])
        if 'xarange' in kwargs:
            ax.set_ylim(kwargs['yarange'])

