"""Shear analysis."""

import numpy
import pylab
import seaborn


def compute_shear(e1r, e2r, e1i, e2i, distx, disty):
    """Compute the shear."""
    phi = numpy.arctan2(disty, distx)
    gamt = - (e1i * numpy.cos(2.0 * phi) + e2i * numpy.cos(2.0 * phi))
    gamc = - e1i * numpy.sin(2.0 * phi) + e2i * numpy.cos(2.0 * phi)
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
    filt = table['modelfit_CModel_mag'][table['filter'] == 'r'] < 23.5  # magnitude cut
    filt &= table['ext_shapeHSM_HsmShapeRegauss_resolution'][table['filter'] == 'i'] \
            > 0.4  # resolution cut
    filt &= (abs(e1i) < 1) & (abs(e2i) < 1)  # ellipticity cut
    filt &= (abs(e1r - e1i) < 0.5) & (abs(e2r - e2i) < 0.5)  # er ~= ei

    # Apply cuts
    e1r, e2r, e1i, e2i, distx, disty = [x[filt] for x in [e1r, e2r, e1i, e2i, distx, disty]]

    # Comput the shear
    gamt, gamc, dist = compute_shear(e1r, e2r, e1i, e2i, distx, disty)

    # Make some plots
    plot_shear(gamt, gamc, dist)


def plot_shear(gamt, gamc, dist, drange=(0, 8500), nbins=8):
    """Plot shear."""

    dval, step = numpy.linspace(drange[0], drange[1], nbins, retstep=True)

    fig = pylab.figure(figsize=(15, 8))
    ax0 = fig.add_subplot(121, xlabel='Gamt')
    ax1 = fig.add_subplot(122, xlabel='Gamc')
    ax0.hist(gamt, bins=200, range=[-2, 2])
    ax1.hist(gamc, bins=200, range=[-2, 2])

    fig = pylab.figure(figsize=(15, 8))
    ax0 = fig.add_subplot(121, xlabel='Gamt')
    ax1 = fig.add_subplot(122, xlabel='Gamc')
    ax0.hist(gamt, bins=80, range=[-0.2, 0.2])
    ax1.hist(gamc, bins=80, range=[-0.2, 0.2])

    fig = pylab.figure(figsize=(15, 8))
    ax0 = fig.add_subplot(121, xlabel='Dist', ylabel='Gamt')
    ax1 = fig.add_subplot(122, xlabel='Dist', ylabel='Gamc')
    ax0.scatter(dist, gamt, s=1, color='b')
    ax0.set_ylim([-1., 1.])
    ax1.scatter(dist, gamc, s=1, color='b')
    ax1.set_ylim([-1., 1.])

    masks = [(dist > d - step / 2) & (dist <= d + step / 2) for d in dval]
    tshear = [numpy.mean(gamt[mask]) for mask in masks]
    cshear = [numpy.mean(gamc[mask]) for mask in masks]
    tsheare = [numpy.std(gamt[mask])/numpy.sqrt(len(gamt[mask])) for mask in masks]
    csheare = [numpy.std(gamc[mask])/numpy.sqrt(len(gamc[mask])) for mask in masks]

    fig = pylab.figure(figsize=(15, 8))
    ax0 = fig.add_subplot(121, ylabel='Tangential shear',
                          xlabel='Distance to cluster center (px)')
    ax1 = fig.add_subplot(122, ylabel='Cross shear',
                          xlabel='Distance to cluster center (px)')
    ax0.errorbar(dval, tshear, yerr=tsheare)
    ax1.errorbar(dval, cshear, yerr=csheare)
    ax0.set_xlim(-500, 9000)
    ax0.set_ylim(-0.06, 0.08)
    ax1.set_xlim(-500, 9000)
    ax1.set_ylim(-0.06, 0.08)

    pylab.show()
