"""Shear analysis."""

import numpy
import pylab
import seaborn


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
        ax = fig.add_subplot(1, len(xs), i+1, xlabel=labels[i])
        ax.hist(x, bins=nbins, range=xarange)

def plot_scatter(xs, ys, xlabels, ylabels, yerrs=None, xarange=None, yarange=None):
    """Plot multiple histogramsscatter plots in subplots."""
    fig = pylab.figure(figsize=(15, 8))
    for i, x in enumerate(xs):
        ax = fig.add_subplot(1, len(xs), i+1, xlabel=xlabels[i], ylabel=ylabels[i])
        ax.axhline(0, color='k', ls=':')
        ax.scatter(x, ys[i], s=1, color='b')
        if yerrs is not None:
            ax.errorbar(x, ys[i], yerr=yerrs[i])
        if xarange is not None:
            ax.set_xlim(xarange)
        if yarange is not None:
            ax.set_ylim(yarange)

