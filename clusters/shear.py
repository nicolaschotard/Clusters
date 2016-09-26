"""Shear analysis."""

import numpy
import pylab
import seaborn


def compute_shear():
    """Compute the shear."""
    pass


def quality_cuts(table, mag='modelfit_CModel_mag'):
    """Apply quality cuts."""

    # keep the 'r' filter
    filt = table[mag] == 'r'

    # resolution cut
    filt = ['ext_shapeHSM_HsmShapeRegauss_resolution']
    indx = np.where((np.asarray(magR) < 23.5) & (gMinusR > yMin1) & (gMinusR < yMax1))
    indx2 = np.where((np.asarray(magR) >= 23.5) | (gMinusR < yMin1) | (gMinusR > yMax1))
    
    rMinusI = np.asarray(magR)-np.asarray(magI)
    yMax2 = -0.033*np.asarray(magR) + 1.54
    yMin2 = -0.033*np.asarray(magR) + 1.35
    indx1p = np.where((np.asarray(magR) < 23.5) & (rMinusI > yMin2) & (rMinusI < yMax2))
    indx3 = np.where((np.asarray(magR) >= 23.5) | (rMinusI < yMin2) | (rMinusI > yMax2))
    indx4 = np.intersect1d(indx2, indx3)
    indx5 = np.where((np.asarray(magR)) < 25)
    indx6 = np.intersect1d(indx4, indx5)
    
    # Define masks in order to reject cluster galaxies from the selection
    mask1 = np.zeros(gMinusR.shape, dtype='bool')
    mask1[indx] = True
    mask2 = np.zeros(gMinusR.shape, dtype='bool')
    mask2[indx1p] = True
    mask3 = np.zeros(gMinusR.shape, dtype='bool')
    mask3[indx5] = True
    mask4 = ~(mask1*mask2) * mask3
    
    indx6 = np.where(mask4)
    indx7 = np.where((mask4) & (np.asarray(res_i)>0.4) & (np.fabs(np.asarray(e1I))<1) & (np.fabs(np.asarray(e2I))<1))
    
    de1 = np.fabs(np.asarray(e1R)-np.asarray(e1I))
    de2 = np.fabs(np.asarray(e2R)-np.asarray(e2I))
    
    mask5 = np.zeros(gMinusR.shape, dtype='bool')
    mask5[indx7] = True
    indx8 = np.where((mask5) & (de1<0.5) & (de2<0.5))


def analysis(table, xclust, yclust):
    """Computethe shear.

    :param string data_file: Name of the hdf5 file to load
    :param string path: Path (key) of the table to load
    :return: A dictionnary containing the following keys and values:

     - meas: the 'deepCoadd_meas' catalog (an astropy table)
     - forced: the 'deepCoad_forced_src' catalog (an astropy table)
     - wcs: the 'wcs' of these catalogs (an ``astropy.wcs.WCS`` object)
    """
    ell1 = table["ext_shapeHSM_HsmShapeRegauss_e1"]
    ell2 = table["ext_shapeHSM_HsmShapeRegauss_e2"]
    xsrc = table["x_Src"][table['filter'] == 'r']
    ysrc = table["y_Src"][table['filter'] == 'r']

    magr = table['modelfit_CModel_mag'][table['filter'] == 'r']
    magi = table['modelfit_CModel_mag'][table['filter'] == 'i']

    resi = data['ext_shapeHSM_HsmShapeRegauss_resolution'][table['filter'] == 'i']

    e1r = ell1[table['filter'] == 'r']
    e2r = ell2[table['filter'] == 'r']
    e1i = ell1[table['filter'] == 'i']
    e2i = ell2[table['filter'] == 'i']

    # magnitude cuts
    filt = magr < 23.5
    
    filt &= resi > 0.4
    filt &= numpy.fabs(e1i) < 1
    filt &= numpy.fabs(e2i) < 1
    
    de1 = numpy.fabs(e1R - e1I)
    de2 = numpy.fabs(e2R - e2I)

    filt &= de1 < 0.5) & (de2 < 0.5)

    mod = numpy.hypot(e1r, e2r)
    tg = e2r / e1r

    print numpy.mean(mod), numpy.median(mod)

    phi = numpy.arctan2(ysrc - yclust, xsrc - xclust)

    gamT = - (e1i * numpy.cos(2.0 * phi) + e2i * numpy.cos(2.0 * phi))
    gamC = - e1i * numpy.sin(2.0 * phi) + e2i * numpy.cos(2.0 * phi)
    dist = numpy.sqrt((xsrc - xclust) * (xsrc - xclust) + (ysrc - yclust) * (ysrc - yclust))

    indx = (gamT < 0.2) & (gamT > -0.2)
    print numpy.median(gamT[indx]), numpy.median(gamC)

    fig, (ax0, ax1) = pylab.subplots(ncols=2, figsize=(15, 8))
    ax0.hist(gamT, bins=200, range=[-2, 2])
    ax1.hist(gamC, bins=200, range=[-2, 2])

    fig, (ax0, ax1) = pylab.subplots(ncols=2, figsize=(15, 8))
    ax0.hist(gamT, bins=80, range=[-0.2, 0.2])
    ax1.hist(gamC, bins=80, range=[-0.2, 0.2])

    fig, (ax0, ax1) = pylab.subplots(ncols=2, figsize=(15, 8))
    ax0.scatter(dist, gamT, s=1, color='b')
    ax0.set_ylim([-1., 1.])
    ax1.scatter(dist, gamC, s=1, color='b')
    ax1.set_ylim([-1., 1.])

    dval, step = numpy.linspace(0, 8500, num=8, retstep=True)

    tShear = []
    cShear = []
    tShearE = []
    cShearE = []
    for d in dval :
        mask = (dist > d - step / 2) & (dist <= d + step / 2)
        tShear.append(numpy.mean(gamT[mask]))
        cShear.append(numpy.mean(gamC[mask]))
        tShearE.append(numpy.std(gamT[mask])/numpy.sqrt(len(gamT[mask])))
        cShearE.append(numpy.std(gamC[mask])/numpy.sqrt(len(gamC[mask])))

    fig, (ax0, ax1) = pylab.subplots(ncols=2, figsize=(15, 8))
    ax0.errorbar(dval,tShear,yerr=tShearE)
    ax1.errorbar(dval,cShear,yerr=cShearE)
    ax0.tick_params(labelsize=15)
    ax0.set_xlabel("Distance to cluster center (px)",fontsize=15)
    ax0.set_ylabel("Tangential shear",fontsize=15)
    ax0.set_xlim(-500, 9000)
    ax0.set_ylim(-0.06, 0.08)
    ax1.tick_params(labelsize=15)
    ax1.set_xlabel("Distance to cluster center (px)",fontsize=15)
    ax1.set_ylabel("Cross shear",fontsize=15)
    ax1.set_xlim(-500, 9000)
    ax1.set_ylim(-0.06, 0.08)
    pylab.show()
