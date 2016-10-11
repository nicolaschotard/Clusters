"""Tools to fit the red sequence and extract background galaxies around a cluster."""

from scipy import optimize
from scipy import special
import numpy as N
import pylab as P
import seaborn
import math

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
    print "R^2 = ", rsquared
    if rsquared < 0.9:
        print "Bad fit - take absolute minimun instead of fitted value"
        fitslope = bestslope
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

    print"Ordinate at origin of the RS band - middle : %f, lower : %f, upper : %f" %\
        (b0, b1, b2)

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

    P.show()


def zphot_cut(zclust, zdata):
    """."""
    zbest = zdata['Z_BEST']
    cbest = zdata['CHI_BEST']
    error = (zbest - zdata['Z_BEST68_LOW'] + zdata['Z_BEST68_HIGH']) / 2.

    # basic filters
    filt = (zbest > 0) & (zbest < 3) & (error < 1) & (cbest < 10)
    print "INFO: Removing %i objects (over %i) from the list" % \
        (len(zbest[~filt]), len(zbest))
    print "       - (zbest > 0) & (zbest < 5) & (error < 1) & (cbest < 10)"
    zbest, cbest, error = [a[filt] for a in [zbest, cbest, error]]
    fig = P.figure()
    ax = fig.add_subplot(131, xlabel='ZBEST')
    ax.hist(zbest, bins=100)
    ax.axvline(zclust, color='r')
    ax.set_title("%i galaxies" % len(zbest))
    ax = fig.add_subplot(132, xlabel='Error')
    ax.hist(error, bins=100)
    ax.set_title("%i galaxies" % len(zbest))
    ax = fig.add_subplot(133, xlabel='CHI_BEST')
    ax.hist(cbest, bins=100)
    ax.set_title("%i galaxies" % len(zbest))

    P.show()

    return zdata['objectId'][filt]


def get_background(config, data, zdata=None, zspec=None):
    """Apply different cuts to the data in order to get the background galaxies."""
    print config['cluster'], len(data)

    # Cut data futher than a given radius around the center of the cluster
    cdata.filter_around(data, config, exclude_outer=20, exclude_inner=3, unit='arcmin')

    # Red sequence
    print "INFO; Getting red sequence"

    # Spectroscopic against photometric redshifts
    if zspec is not None:
        print "INFO: Checking photo/spectro redshifts consitancy"

    # Photometric redshift cut
    if zdata is not None:
        zdata = cdata.read_data(zdata)
        print "INFO: A redshift cut will be applied."
        zphot_cut(config['redshift'], zdata)
