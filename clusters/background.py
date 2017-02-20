"""Tools to fit the red sequence and extract background galaxies around a cluster."""

import math
from scipy import optimize, special
from astropy.cosmology import Planck15 as cosmo
from astropy import units as u
import numpy as N
import pylab as P
from . import data as cdata


def color_histo(mags):
    """Plot color histograms."""
    filt = (mags['g'] - mags['r']) > 1.2
    for i, filt1 in enumerate('gri'):
        for j, filt2 in enumerate('gri'):
            if i >= j:
                continue
            fig = P.figure()
            ax = fig.add_subplot(111)
            ax.hist((mags[filt1] - mags[filt2])[filt],
                    bins=100, label='%s - %s' % (filt1, filt2))
            ax.legend(loc='best')
    P.show()


def color_mag_plot(mags):
    """Plot color / mag diagrams."""
    filt = (mags['g'] - mags['r']) > 1.2
    for fref in enumerate('gri'):
        for i, filt1 in enumerate('gri'):
            for j, filt2 in enumerate('gri'):
                if i >= j:
                    continue
                fig = P.figure()
                ax = fig.add_subplot(111)
                ax.scatter(mags[fref][filt], (mags[filt1] - mags[filt2])[filt],
                           s=1, label='%s - %s' % (filt1, filt2))
                ax.set_xlabel(fref)
                ax.set_ylabel('%s - %s' % (filt1, filt2))
                ax.legend(loc='best')
    P.show()


def fit_red_sequence(color, mag, **kwargs):
    r"""
    Fit red sequence (RS) band in galaxy color plots, i.e m(i)-m(j) vs. m(k).

    :param list color: A list of color mag_i - mag_j (ordinate)
    :param list mag: List of magnitude (abciss)
    :param \*\*kwargs:

     - minc (float): lower cut on the color axis: color > minc (1.0)
     - maxc (float): upper cut of the color axis: color < maxc  (2.0)
     - minm (float): lower cut on the mag axis: mag > minm (20.0)
     - maxm (float): upper cut of the mag axis: mag < maxm  (23.5)
     - islope (float): first guess for the red sequence band slope (-0.04)
     - nbins (int): Number of bins used in the fits (40)
     - plot (bool): if True plot stuff
     - verbose (bool): if True print information to screen
    :return:

     - slope of the red sequence band,
     - ordinate at the origin of the red sequence band +/- 1.5 sigma

    fitRedSequence is also producing some control plots
    """
    # Set the default
    minc = kwargs.get('minc', 1.0)
    maxc = kwargs.get('maxc', 2.0)
    minm = kwargs.get('minm', 20.0)
    maxm = kwargs.get('maxm', 23.5)
    islope = kwargs.get('islope', -0.04)
    nbins = kwargs.get('nbins', 40)
    plot = kwargs.get('plot', False)
    verb = kwargs.get('verbose', False)

    magref = minm  # Arbitrary reference magnitude for projection
    diffref = 0.5 * (minc + maxc)  # Arbitrary reference ordinate for projection

    # Project color on an axis perpendicular to the RS band and at an
    # arbitrary magnitude (magref)
    # The idea is that the projection of the RS band on this axis is a gaussian
    # over some background represented by an Exponentially Modified Gaussian

    alpha = math.atan(islope)
    dy = N.cos(alpha) * ((N.asarray(magref) - mag) * islope + color - diffref)

    idx = N.where((color > minc) & (color < maxc) & (mag < maxm) & (mag > minm))
    fig, (ax0) = P.subplots(ncols=1)
    n, bins, patches = ax0.hist(dy[idx], bins=nbins, color='b')

    x = N.asarray([0.5 * (bins[i + 1] - bins[i]) + bins[i] for i in range(len(n))])

    # Fit a gaussian for the RS projection plus an Exponentially Modified
    # Gaussian distribution for the background
    def func(p, z):
        """Function to fit."""
        return p[0] * N.exp(-(p[1] - z)**2 / (2 * p[2]**2)) + \
            p[6] * p[5] * N.exp(0.5 * p[5] * (2 * p[3] + p[5] * p[4]**2 - 2 * z)) * \
            special.erfc((p[3] + p[5] * p[4] ** 2 - z) / (math.sqrt(2) * p[4]))

    def dist(p, z, y):
        return (func(p, z) - y) / (N.sqrt(y) + 1.)

    p0 = [n.max(), 0.2, 0.1, -3.0, 1., 1., 40.]  # Initial parameter values
    p2, cov, infodict, mesg, ier = optimize.leastsq(dist, p0[:], args=(x, n), full_output=True)
    ss_err = (infodict['fvec']**2).sum()

    if verb:
        print "mean %f - sigma %f" % (p2[1], p2[2])
        print "Reduced chi2 = ", ss_err / (nbins + 6 - 1)
        print p2
    # Superimpose fitted curve over color projection
    ax0.plot(bins, func(p2, bins), color='g')
    ax0.tick_params(labelsize=20)
    ax0.set_xlabel("Projected red sequence")

    # Minimize the width of the gaussian by varying the RS slope around the
    # initial guess

    nsteps = 80
    step = 0.001
    slope = islope - 0.5 * nsteps * step

    val = []
    sigma = []
    sigmamin = 999.0

    def lfunc(p, z):
        return p[0] * N.exp(-(p[1] - z) ** 2 / (2 * p[2] ** 2)) + \
            p[6] * p[5] * N.exp(0.5 * p[5] * (2 * p[3] + p[5] * p[4] ** 2 - 2 * z)) * \
            special.erfc((p[3] + p[5] * p[4] ** 2 - z) / (math.sqrt(2) * p[4]))

    def ldist(p, z, y):
        return (lfunc(p, z) - y) / (N.sqrt(y) + 1.)

    for i in range(nsteps):
        slope += step
        # RS slope is always negative
        if slope > 0.0:
            break
        alpha = math.atan(slope)
        dy = N.cos(alpha) * ((magref - mag) * slope + color - diffref)
        idx = N.where((color > minc) & (color < maxc) & (mag < maxm) & (mag > minm))
        n, bins, = N.histogram(dy[idx], bins=nbins)

        x = N.asarray([0.5 * (bins[i + 1] - bins[i]) + bins[i] for i in range(len(n))])
        p0 = p2  # Start fit with parameters fitted at the previous step
        p1, cov, infodict, mesg, ier = optimize.leastsq(ldist, p0[:], args=(x, n), full_output=True)
        val.append(slope)
        sigma.append(p1[2])
        if p1[2] < sigmamin:
            sigmamin = p1[2]
            param = p1

    fig, (ax0) = P.subplots(ncols=1)
    ax0.scatter(val, sigma, s=5, color='b')
    bestslope = val[N.argmin(N.asarray(sigma))]

    # Fit a parabola on the (slope, sigma) distribution to find the RS slope
    # corresponding to the minimum sigma
    def parabola(p, z):
        return p[0] * z * z + p[1] * z + p[2]

    def pdist(p, z, y):
        return parabola(p, z) - y

    p0 = [1., 1., 1.]
    p1, cov, infodict1, mesg, ier = optimize.leastsq(pdist, p0[:],
                                                     args=(N.asarray(val),
                                                           N.asarray(sigma)),
                                                     full_output=True)
    fitslope = -0.5 * p1[1] / p1[0]

    ax0.plot(N.asarray(val), parabola(p1, N.asarray(val)), color='r')
    ax0.tick_params(labelsize=20)
    ax0.set_xlabel("Red sequence slope")
    ax0.set_ylabel("Sigma")

    ss_err = (infodict1['fvec']**2).sum()
    ss_tot = ((N.asarray(sigma) - N.asarray(sigma).mean())**2).sum()
    rsquared = 1 - (ss_err / ss_tot)
    if verb:
        print "R^2 = ", rsquared
    if rsquared < 0.9:
        if verb:
            print "Bad fit - take absolute minimun instead of fitted value"
        fitslope = bestslope
    if verb:
        print "Fitted minimum: %f" % fitslope

    # Plot RS projection corresponding to the optimal slope
    alpha = math.atan(slope)
    dy = N.cos(alpha) * ((magref - mag) * fitslope + color - diffref)
    idx = N.where((color > minc) & (color < maxc) & (mag < maxm) & (mag > minm))
    fig, (ax2) = P.subplots(ncols=1)
    n, bins, patches = ax2.hist(dy[idx], bins=nbins, color='b')
    x = N.asarray([0.5 * (bins[i + 1] - bins[i]) + bins[i] for i in range(len(n))])

    def hfunc(p, z):
        """Function to fit."""
        return p[0] * N.exp(-(p[1] - z) ** 2 / (2 * p[2] ** 2)) + \
            p[6] * p[5] * N.exp(0.5 * p[5] * (2 * p[3] + p[5] * p[4] ** 2 - 2 * z)) * \
            special.erfc((p[3] + p[5] * p[4] ** 2 - z) / (math.sqrt(2) * p[4]))

    def hdist(p, z, y):
        return (hfunc(p, z) - y) / (N.sqrt(y) + 1.)

    p0 = param
    p1, cov, infodict, mesg, ier = optimize.leastsq(hdist, p0[:], args=(x, n), full_output=True)
    ss_err = (infodict['fvec']**2).sum()
    ss_tot = ((n - n.mean()) ** 2).sum()
    rsquared = 1 - (ss_err / ss_tot)
    if verb:
        print "mean %f - sigma %f" % (p1[1], p1[2])
        print "Reduced chi2 = %f - R^2 = %f" % (ss_err / (nbins + 6 - 1), rsquared)
    ax2.plot(bins, hfunc(p1, bins), color='r')
    ax2.tick_params(labelsize=20)

    # Compute the ordinates at origin corresponding to a +/- 1.5 sigma
    # interval around the best gaussian mean value
    alpha = math.atan(fitslope)
    b0 = (p1[1] - fitslope * magref) / math.cos(alpha) + diffref
    b1 = (p1[1] - 1.5 * p1[2] - fitslope * magref) / math.cos(alpha) + diffref
    b2 = (p1[1] + 1.5 * p1[2] - fitslope * magref) / math.cos(alpha) + diffref

    if verb:
        print"Ordinate at origin of the RS band - middle : %f, lower : %f, upper : %f" %(b0, b1, b2)

    # plot fitted RS band over color plot
    fig, (ax3) = P.subplots(ncols=1)
    ax3.scatter(mag, color, s=1, color='b')
    ax3.set_xlim([minm - 3.0, maxm + 3.0])
    ax3.set_ylim([minc - 1.5, maxc + 0.5])
    x = N.linspace(minm - 3.0, maxm + 3.0)
    ymin = fitslope * x + b1
    ymax = fitslope * x + b2
    ymid = fitslope * x + b0
    ax3.plot(x, ymin, color='r')
    ax3.plot(x, ymax, color='r')
    ax3.plot(x, ymid, color='g')

    if plot is True:
        P.show()

    params = [[fitslope, b1], [fitslope, b2]]
    return params


def zphot_cut(zclust, zdata, **kwargs):
    r"""
    Redshif selection of the galaxies used for analysis, using both:
    - hard cut, z_cl+0.1 < z_best < 1.25 (cf WtGIII)
    - cut from pdz. \int_0^z_cl p(z) dz < x%

    :param float plot: if keywords exists, plot stuff for visual inspection
    :param float thresh: tolerance x% for the pdz cut method.
    Returns bool arrays, where False means the object does not pass the cut
    """
    plot = kwargs.get('plot', False)
    thresh = kwargs.get('thresh')
    zmin = kwargs.get('zmin')
    zmax = kwargs.get('zmax')
    zbest = zdata['Z_BEST']
    pdz = zdata['pdz']
    zbins = zdata['zbins'][0]  # all objects have same zbins, take 0th.

    # WtGIII hard cuts
    filt1 = (zbest > zmin) & (zbest < zmax)

    # pdz_based cut
    cut = (zbins < zclust + 0.1)
    # probability for the cluster to be located below zclust + 0.1
    filt2 = N.array([N.trapz(pdzi[cut], zbins[cut]) * 100. < thresh for pdzi in pdz])

    if plot:
        fig = P.figure()
        ax = fig.add_subplot(121, xlabel='ZBEST')
        ax.hist(zbest, bins=100)
        ax.axvline(zclust + 0.1, color='r')
        ax.axvline(1.25, color='r')
        ax.set_title("%i galaxies" % len(zbest))
        P.show()

    return (filt1, filt2)


def red_sequence_cut(config, data, **kwargs):
    """
    Identify RS galaxies using color-magnitude diagram.

    First do a radial cut on catalogue and identify the RS from the inner galaxies
    --> increase the contrast of the RS
    Then go back and apply the cut on the entire catalogue as some RS galaxies are
    located far away from the centre

    Returns bool array, where False means the object does not pass the cut

    List of available kwargs:

    :param float mag_cut: rband magnitude cut - default is 25
    :param float plot: if keywords exists, plot stuff for visual inspection
    """

    mcut = kwargs.get('mag_cut', 25.)
    plot = kwargs.get('plot', False)

    da = cosmo.angular_diameter_distance(config['redshift'])  # Mpc - using Planck15 cosmo
    rcut_rs = 1 * u.Mpc
    sub_sample = cdata.filter_around(data, config, exclude_outer=N.arctan(rcut_rs / da).value,
                                     unit='rad', plot=plot)

    color_gr = sub_sample['modelfit_CModel_mag'][sub_sample['filter'] == 'g'] \
               - sub_sample['modelfit_CModel_mag'][sub_sample['filter'] == 'r']
    mag = sub_sample['modelfit_CModel_mag'][sub_sample['filter'] == 'r']
    params = fit_red_sequence(color_gr, mag, plot=plot)  # slopes and intercepts of the RS band

    # apply cut to entire dataset
    color_gr = data['modelfit_CModel_mag'][data['filter'] == 'g'] \
               - data['modelfit_CModel_mag'][data['filter'] == 'r']
    mag = data['modelfit_CModel_mag'][data['filter'] == 'r']
    lower_bound = params[0][0] * mag + params[0][1]
    upper_bound = params[1][0] * mag + params[1][1]
    filt = ((color_gr < lower_bound) & (mag < mcut)) | ((color_gr > upper_bound) & (mag < mcut))

    return filt


def get_zphot_background(config, zdata, zspec=None, z_config=None, thresh=None, zmin=None, zmax=None, plot=None):
    """Return flag based on zphot criterion for galaxy selection."""

    # Spectroscopic against photometric redshifts
    if zspec is not None:
        print "INFO: Checking photo/spectro redshifts consistency"

    # Photometric redshift cut
    print "INFO: Flagging foreground/uncertain objects using redshift information from ", z_config
    z_flag1, z_flag2 = zphot_cut(config['redshift'], zdata, thresh=thresh,
                                 zmin=zmin, zmax=zmax, plot=plot)
    print "INFO: %i galaxies have been kept after the hard redshift cut" %(sum(z_flag1))
    print "INFO: %i galaxies have been kept after the pdz redshift cut" %(sum(z_flag2))

    return (z_flag1, z_flag2)

def get_rs_background(config, data):
    """Return flag based on RS criterion for galaxy selection."""

    print "INFO: Flagging red sequence galaxies"
    rs_flag = red_sequence_cut(config, data)
    print "INFO: %i galaxies have been flagged as RS" %(sum(~rs_flag))

    return rs_flag
