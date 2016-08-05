from scipy import optimize
from scipy import special
import numpy as N
import pylab as P
import seaborn
import math

def red_sequence(mags, fref, f1, f2, clip=2, colorcut=1):
    m, c = mags[fref][(mags[f1]-mags[f2])>colorcut], (mags[f1]-mags[f2])[(mags[f1]-mags[f2])>colorcut]
    op=Optimizer.LinearRegression(m, c, clip=clip)
    fig, ax = P.subplots(ncols=1)
    ax.plot(N.linspace(min(mags[fref]),max(mags[fref]), 1000),
           N.linspace(min(mags[fref]),max(mags[fref]), 1000)*op.a+op.b, 'g', ls='--')
    ax.plot(N.linspace(min(mags[fref]),max(mags[fref]), 1000),
           N.linspace(min(mags[fref]),max(mags[fref]), 1000)*op.a+op.b+0.2, 'g')
    ax.plot(N.linspace(min(mags[fref]),max(mags[fref]), 1000),
           N.linspace(min(mags[fref]),max(mags[fref]), 1000)*op.a+op.b-0.2, 'g')
    ax.scatter(mags[fref], mags[f1]-mags[f2], s=1, color='k')
    ax.set_xlabel(fref)
    ax.set_ylabel('%s - %s' % (f1, f2))
    print op.a, op.b
    return op

def color_histo(mags):
    filt = (mags['g'] - mags['r']) > 1.2
    for i, f1 in enumerate('gri'):
        for j, f2 in enumerate('gri'):
            if i >= j :
                continue
            fig, ax = P.subplots(ncols=1)
            ax.hist((mags[f1]-mags[f2])[filt],
                    bins=100, label='%s - %s' % (f1, f2))
            ax.legend(loc='best')
    P.show()

def color_mag_plot(mags):
    filt = (mags['g'] - mags['r']) > 1.2
    for k, fref in enumerate('gri'):
        for i, f1 in enumerate('gri'):
            for j, f2 in enumerate('gri'):
                if i >= j:
                    continue
                fig, ax = P.subplots(ncols=1)
                ax.scatter(mags[fref][filt], (mags[f1]-mags[f2])[filt],
                           s=1, label='%s - %s' % (f1, f2))
                ax.set_xlabel(fref)
                ax.set_ylabel('%s - %s' % (f1, f2))
                ax.legend(loc='best')
    P.show()

def fitRedSequence(diffMag, mag, minDiff=1.0, maxDiff=2.0, minMag=20.0, maxMag=23.5, iniSlope=-0.04):
    """
    fit red sequence (RS) band in galaxy color plots
    color plots are typically mag_i-mag_j versus mag_k
    Input :
    diffMag : mag_i-mag_j axis of the color plot
    mag : abcissa of the color plot
    minDiff, maxDiff : select minDiff < diffMag < maxDiff
    minMag, maxMag : select minMag < mag < maxMag
    iniSlope : first guess for the RS band slope

    Output :
    slope of the red sequence band
    ordinate at the origin of the red sequence band +/- 1.5 sigma

    fitRedSequence is also producing some control plots

    """

    magRef = minMag # Arbitrary reference magnitude for projection
    diffRef = 0.5*(minDiff+maxDiff) # Arbitrary reference ordinate for projection
    appSlope = iniSlope

    # Project diffMag on an axis perpendicular to the RS band and at an
    # arbitrary magnitude (magRef)
    # The idea is that the projection of the RS band on this axis is a gaussian
    # over some background represented by an Exponentially Modified Gaussian

    alpha = math.atan(appSlope)
    dy = N.cos(alpha)*((N.asarray(magRef)-mag)*appSlope + diffMag - diffRef)

    nbins = 40
    idx = N.where( (diffMag > minDiff) & (diffMag < maxDiff) & (mag < maxMag) & (mag > minMag) )
    fig, (ax0) = P.subplots(ncols=1)
    n, bins, patches = ax0.hist(dy[idx], bins=nbins, color='b')

    max = n.max()
    x = N.asarray([0.5*(bins[i+1]-bins[i])+bins[i] for i in range(len(n))])

    # Fit a gaussian for the RS projection plus an Exponentially Modified
    # Gaussian distribution for the background
    func = (lambda p, z: p[0]*N.exp(-(p[1]-z)**2/(2*p[2]**2))+
            p[6]*p[5]*N.exp(0.5*p[5]*(2*p[3]+p[5]*p[4]**2-2*z))*
            special.erfc((p[3]+p[5]*p[4]**2-z)/(math.sqrt(2)*p[4])))
    dist = lambda p, z, y: (func(p, z) - y)/(N.sqrt(y)+1.)
    p0 = [max, 0.2, 0.1, -3.0, 1., 1., 40.]  # Initial parameter values
    p2,cov,infodict,mesg,ier = optimize.leastsq(dist, p0[:], args=(x, n), full_output=True)
    ss_err=(infodict['fvec']**2).sum()
    print("mean %f - sigma %f"%(p2[1], p2[2]))
    print("Reduced chi2 = ", ss_err/(nbins+6-1))
    print p2
    # Superimpose fitted curve over diffMag projection
    ax0.plot(bins, func(p2,bins), color='g')
    ax0.tick_params(labelsize=20)
    ax0.set_xlabel("Projected red sequence")

    # Minimize the width of the gaussian by varying the RS slope around the
    # initial guess

    nSteps = 80
    step = 0.001
    slope = appSlope-0.5*(nSteps)*step

    val = []
    sigma = []
    sigmaMin = 999.0

    for i in range(nSteps) :
        slope = slope + step
        # RS slope is always negative
        if slope > 0.0 :
            break
        alpha = math.atan(slope)
        dy = N.cos(alpha)*((magRef-mag)*slope + diffMag - diffRef)
        nbins = 40
        idx = N.where( (diffMag > minDiff) & (diffMag < maxDiff) & (mag < maxMag) & (mag > minMag) )
        n, bins, = N.histogram(dy[idx], bins=nbins)

        x = N.asarray([0.5*(bins[i+1]-bins[i])+bins[i] for i in range(len(n))])

        func = (lambda p, z: p[0]*N.exp(-(p[1]-z)**2/(2*p[2]**2))+
                p[6]*p[5]*N.exp(0.5*p[5]*(2*p[3]+p[5]*p[4]**2-2*z))*
                special.erfc((p[3]+p[5]*p[4]**2-z)/(math.sqrt(2)*p[4])))
        dist = lambda p, z, y: (func(p, z) - y)/(N.sqrt(y)+1.)
        p0 = p2 # Start fit with parameters fitted at the previous step
        p1,cov,infodict,mesg,ier = optimize.leastsq(dist, p0[:], args=(x, n), full_output=True)
        val.append(slope)
        sigma.append(p1[2])
        if p1[2] < sigmaMin :
            sigmaMin = p1[2]
            param = p1

    fig, (ax0) = P.subplots(ncols=1)
    ax0.scatter(val, sigma, s=5, color='b')
    bestSlope = val[N.argmin(N.asarray(sigma))]

    # Fit a parabola on the (slope, sigma) distribution to find the RS slope
    # corresponding to the minimum sigma
    func = lambda p, z: p[0]*z*z + p[1]*z + p[2]
    dist = lambda p, z, y: (func(p, z) - y)
    p0 = [1., 1., 1.]
    p1,cov,infodict1,mesg,ier = optimize.leastsq(dist, p0[:], args=(N.asarray(val), N.asarray(sigma)), full_output=True)
    fitSlope = -0.5*p1[1]/p1[0]

    ax0.plot(N.asarray(val), func(p1,N.asarray(val)), color='r')
    ax0.tick_params(labelsize=20)
    ax0.set_xlabel("Red sequence slope")
    ax0.set_ylabel("Sigma")

    ss_err=(infodict1['fvec']**2).sum()
    ss_tot=((N.asarray(sigma)-N.asarray(sigma).mean())**2).sum()
    rsquared=1-(ss_err/ss_tot)
    print "R^2 = ", rsquared
    if rsquared < 0.9 :
        print "Bad fit - take absolute minimun instead of fitted value"
        fitSlope = bestSlope
    print("Fitted minimum: %f"%fitSlope)

    # Plot RS projection corresponding to the optimal slope
    nbins = 40
    alpha = math.atan(slope)
    dy = N.cos(alpha)*((magRef-mag)*fitSlope + diffMag - diffRef)
    idx = N.where( (diffMag > minDiff) & (diffMag < maxDiff) & (mag < maxMag) & (mag > minMag) )
    fig, (ax2) = P.subplots(ncols=1)
    n, bins, patches = ax2.hist(dy[idx], bins=nbins, color='b')
    x = N.asarray([0.5*(bins[i+1]-bins[i])+bins[i] for i in range(len(n))])

    func = (lambda p, z: p[0]*N.exp(-(p[1]-z)**2/(2*p[2]**2))+
                p[6]*p[5]*N.exp(0.5*p[5]*(2*p[3]+p[5]*p[4]**2-2*z))*
                special.erfc((p[3]+p[5]*p[4]**2-z)/(math.sqrt(2)*p[4])))
    dist = lambda p, z, y: (func(p, z) - y)/(N.sqrt(y)+1.)
    p0 = param
    p1,cov,infodict,mesg,ier = optimize.leastsq(dist, p0[:], args=(x, n), full_output=True)
    ss_err=(infodict['fvec']**2).sum()
    ss_tot=((n-n.mean())**2).sum()
    rsquared=1-(ss_err/ss_tot)
    print("mean %f - sigma %f"%(p1[1], p1[2]))
    print("Reduced chi2 = %f - R^2 = %f"%(ss_err/(nbins+6-1), rsquared))
    ax2.plot(bins, func(p1,bins), color='r')
    ax2.tick_params(labelsize=20)

    # Compute the ordinates at origin corresponding to a +/- 1.5 sigma
    # interval around the best gaussian mean value
    alpha = math.atan(fitSlope)
    b0 = (p1[1] - fitSlope*magRef)/math.cos(alpha) + diffRef
    b1 = (p1[1]-1.5*p1[2] - fitSlope*magRef)/math.cos(alpha) + diffRef
    b2 = (p1[1]+1.5*p1[2] - fitSlope*magRef)/math.cos(alpha) + diffRef

    print("Ordinate at origin of the RS band - middle : %f, lower : %f, upper : %f"%(b0, b1, b2))

    # plot fitted RS band over color plot
    fig, (ax3) = P.subplots(ncols=1)
    ax3.scatter(mag, diffMag, s=1, color='b')
    ax3.set_xlim([minMag-3.0,maxMag+3.0])
    ax3.set_ylim([minDiff-1.5, maxDiff+0.5])
    x = N.linspace(minMag-3.0, maxMag+3.0)
    yMin = fitSlope*x + b1
    yMax = fitSlope*x + b2
    yMid = fitSlope*x + b0
    ax3.plot(x, yMin, color='r')
    ax3.plot(x, yMax, color='r')
    ax3.plot(x, yMid, color='g')

    P.show()
