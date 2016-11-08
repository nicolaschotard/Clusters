#!/usr/bin/env python
##########################
# From Joerg Dietrich 2013
# Modified by Douglas Applegate 2013
############################

import numpy as np

def positionAngle(ra1, dec1, ra2, dec2):
    """Compute the astronomical position angle of position 2 wrt to
    position 1. Inputs are expected in decimal degrees. Returns
    position angle in radians."""
    ra1rad = ra1 * np.pi / 180.
    ra2rad = ra2 * np.pi / 180.
    dec1rad = dec1 * np.pi / 180.
    dec2rad = dec2 * np.pi / 180.
    raDelta = ra2rad - ra1rad
    theta = np.arctan2(np.sin(raDelta), np.cos(dec1rad) * np.tan(dec2rad) \
                           - np.sin(dec1rad) * np.cos(raDelta))
    return theta

def positionAngleCartesian(ra1, dec1, ra2, dec2):
    ra1rad = ra1 * np.pi / 180.
    ra2rad = ra2 * np.pi / 180.
    dec1rad = dec1 * np.pi / 180.
    dec2rad = dec2 * np.pi / 180.
    raDelta = ra2rad - ra1rad
    decDelta = dec2rad - dec1rad
    return np.arctan2(decDelta, raDelta)

def greatCircleDistance(ra1, dec1, ra2, dec2):
    """Compute the great circle distance from poistion 1 to position
    2. Inputs are expected in decimal degrees. Returns great circle
    distance in decimal degrees. Accurate for small separations. This
    routine uses the Vincenty equation, which also handles antipodal
    point correctly."""
    phi1 = dec1 * np.pi / 180.
    lam1 = ra1 * np.pi / 180.
    phi2 = dec2 * np.pi / 180.
    lam2 = ra2 * np.pi / 180.
    deltaLam = np.abs(lam1 - lam2)
    n1 = (np.cos(phi2) * np.sin(deltaLam))**2
    n2 = (np.cos(phi1) * np.sin(phi2) - 
          np.sin(phi1) * np.cos(phi2) * np.cos(deltaLam))**2
    num = np.sqrt(n1 + n2)
    d1 = np.sin(phi1) * np.sin(phi2)
    d2 = np.cos(phi1) * np.cos(phi2) * np.cos(deltaLam)
    denom = d1 + d2
    return np.arctan2(num, denom) * 180. / np.pi

def duffyConcentration(m, z, overdensity):
    """Compute the Duffy et al. mass concentration relation for a halo
    with
    m           - in M_sun/h
    z           - redshift
    overdensity - wrt to the critical density
    """
    n = 1e5
    A = 5.71
    B = -0.084
    C = -0.47
    # Get a first estimate of the concentration to convert to M200crit
    c0 = A * (m / 2e12)**B * (1. + z)**C
    delta200c = overdensity / 200.
    minval = 0.2
    maxval = 5.
    ratio = minval + np.arange(n) / n * (maxval - minval)
    res = ratioResid(ratio, c0, delta200c)
    rRatio = ratio[(res**2).argmin()]
    mRatio = delta200c * rRatio**3
    m200c = m / mRatio
    # Convert input to M200c using the updated concentration
    nIter = 2
    for i in range(nIter):
        c = A * (m200c / 2e12)**B * (1. + z)**C
        res = ratioResid(ratio, c, delta200c)
        rRatio = ratio[(res**2).argmin()]
        mRatio = delta200c * rRatio**3
        m200c = m / mRatio
    return A * (m200c / 2e12)**B * (1. + z)**C

