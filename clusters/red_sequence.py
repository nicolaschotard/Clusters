from scipy import optimize
from scipy import special
import numpy as N
import pylab as P
import seaborn
import math

def color_histo(mags):
    """Plot color histograms."""
    filt = (mags['g'] - mags['r']) > 1.2
    for i, filt1 in enumerate('gri'):
        for j, filt2 in enumerate('gri'):
            if i >= j:
                continue
            fig, ax = P.subplots(ncols=1)
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
                fig, ax = P.subplots(ncols=1)
                ax.scatter(mags[fref][filt], (mags[filt1]-mags[filt2])[filt],
                           s=1, label='%s - %s' % (filt1, filt2))
                ax.set_xlabel(fref)
                ax.set_ylabel('%s - %s' % (filt1, filt2))
                ax.legend(loc='best')
    P.show()

def fit_red_sequence(diffmag, mag, mindiff=1.0, maxdiff=2.0,
                     minmag=20.0, maxmag=23.5, inislope=-0.04):
    """
    Fit red sequence (RS) band in galaxy color plots, i.e m(i)-m(j) vs. m(k).

    Input :
    diffmag : mag_i-mag_j axis of the color plot
    mag : abcissa of the color plot
    mindiff, maxdiff : select mindiff < diffmag < maxdiff
    minmag, maxmag : select minmag < mag < maxmag
    inislope : first guess for the RS band slope

    Output :
    slope of the red sequence band
    ordinate at the origin of the red sequence band +/- 1.5 sigma

    fitRedSequence is also producing some control plots
    """
    magref = minmag  # Arbitrary reference magnitude for projection
    diffref = 0.5 * (mindiff + maxdiff)  # Arbitrary reference ordinate for projection
    appslope = inislope

    # Project diffmag on an axis perpendicular to the RS band and at an
    # arbitrary magnitude (magref)
    # The idea is that the projection of the RS band on this axis is a gaussian
    # over some background represented by an Exponentially Modified Gaussian

    alpha = math.atan(appslope)
    dy = N.cos(alpha) * ((N.asarray(magref) - mag) * appslope + diffmag - diffref)

    nbins = 40
    idx = N.where((diffmag > mindiff) & (diffmag < maxdiff) & (mag < maxmag) & (mag > minmag))
    fig, (ax0) = P.subplots(ncols=1)
    n, bins, patches = ax0.hist(dy[idx], bins=nbins, color='b')

    x = N.asarray([0.5 * (bins[i + 1] - bins[i]) + bins[i] for i in range(len(n))])

    # Fit a gaussian for the RS projection plus an Exponentially Modified
    # Gaussian distribution for the background
    def func(p, z):
        """Function to fit"""
        return p[0] * N.exp(-(p[1] - z)**2 / (2 * p[2]**2)) + \
            p[6] * p[5] * N.exp(0.5 * p[5] * (2 * p[3] + p[5] * p[4]**2 - 2 * z)) * \
            special.erfc((p[3] + p[5] * p[4]**2 - z)/(math.sqrt(2) * p[4]))

    dist = lambda p, z, y: (func(p, z) - y)/(N.sqrt(y) + 1.)
    p0 = [n.max(), 0.2, 0.1, -3.0, 1., 1., 40.]  # Initial parameter values
    p2, cov, infodict, mesg, ier = optimize.leastsq(dist, p0[:], args=(x, n), full_output=True)
    ss_err = (infodict['fvec']**2).sum()
    print "mean %f - sigma %f" % (p2[1], p2[2])
    print "Reduced chi2 = ", ss_err / (nbins + 6 - 1)
    print p2
    # Superimpose fitted curve over diffmag projection
    ax0.plot(bins, func(p2, bins), color='g')
    ax0.tick_params(labelsize=20)
    ax0.set_xlabel("Projected red sequence")

    # Minimize the width of the gaussian by varying the RS slope around the
    # initial guess

    nsteps = 80
    step = 0.001
    slope = appslope - 0.5 * nsteps * step

    val = []
    sigma = []
    sigmamin = 999.0

    for i in range(nsteps):
        slope = slope + step
        # RS slope is always negative
        if slope > 0.0:
            break
        alpha = math.atan(slope)
        dy = N.cos(alpha) * ((magref - mag) * slope + diffmag - diffref)
        nbins = 40
        idx = N.where((diffmag > mindiff) & (diffmag < maxdiff) & (mag < maxmag) & (mag > minmag))
        n, bins, = N.histogram(dy[idx], bins=nbins)

        x = N.asarray([0.5 * (bins[i + 1] - bins[i]) + bins[i] for i in range(len(n))])

        func = (lambda p, z: p[0] * N.exp(-(p[1]-z)**2/(2*p[2]**2))+
                p[6]*p[5]*N.exp(0.5*p[5]*(2*p[3]+p[5]*p[4]**2-2*z))*
                special.erfc((p[3]+p[5]*p[4]**2-z)/(math.sqrt(2)*p[4])))
        dist = lambda p, z, y: (func(p, z) - y)/(N.sqrt(y)+1.)
        p0 = p2 # Start fit with parameters fitted at the previous step
        p1, cov, infodict, mesg, ier = optimize.leastsq(dist, p0[:], args=(x, n), full_output=True)
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
    func = lambda p, z: p[0]*z*z + p[1]*z + p[2]
    dist = lambda p, z, y: (func(p, z) - y)
    p0 = [1., 1., 1.]
    p1,cov,infodict1,mesg,ier = optimize.leastsq(dist, p0[:], args=(N.asarray(val), N.asarray(sigma)), full_output=True)
    fitslope = -0.5 * p1[1] / p1[0]

    ax0.plot(N.asarray(val), func(p1,N.asarray(val)), color='r')
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
    nbins = 40
    alpha = math.atan(slope)
    dy = N.cos(alpha) * ((magref - mag) * fitslope + diffmag - diffref)
    idx = N.where((diffmag > mindiff) & (diffmag < maxdiff) & (mag < maxmag) & (mag > minmag))
    fig, (ax2) = P.subplots(ncols=1)
    n, bins, patches = ax2.hist(dy[idx], bins=nbins, color='b')
    x = N.asarray([0.5*(bins[i+1]-bins[i])+bins[i] for i in range(len(n))])

    def func(p, z):
        """Function to fit."""
        return p[0]*N.exp(-(p[1]-z)**2/(2*p[2]**2)) + \
            p[6]*p[5]*N.exp(0.5*p[5]*(2*p[3]+p[5]*p[4]**2-2*z)) * \
            special.erfc((p[3]+p[5]*p[4]**2-z)/(math.sqrt(2)*p[4]))

    dist = lambda p, z, y: (func(p, z) - y) / (N.sqrt(y) + 1.)
    p0 = param
    p1, cov, infodict, mesg, ier = optimize.leastsq(dist, p0[:], args=(x, n), full_output=True)
    ss_err = (infodict['fvec']**2).sum()
    ss_tot = ((n-n.mean())**2).sum()
    rsquared = 1 - (ss_err / ss_tot)
    print "mean %f - sigma %f" % (p1[1], p1[2])
    print "Reduced chi2 = %f - R^2 = %f" % (ss_err / (nbins + 6-1), rsquared)
    ax2.plot(bins, func(p1, bins), color='r')
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
    ax3.scatter(mag, diffmag, s=1, color='b')
    ax3.set_xlim([minmag - 3.0, maxmag + 3.0])
    ax3.set_ylim([mindiff - 1.5, maxdiff + 0.5])
    x = N.linspace(minmag - 3.0, maxmag + 3.0)
    ymin = fitslope * x + b1
    ymax = fitslope * x + b2
    ymid = fitslope * x + b0
    ax3.plot(x, ymin, color='r')
    ax3.plot(x, ymax, color='r')
    ax3.plot(x, ymid, color='g')

    P.show()
