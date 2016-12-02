#!/usr/bin/env python
##########################
# From Joerg Dietrich 2013
# Modified by Douglas Applegate 2013
############################

import numpy as np


def positionAngle(ra1, dec1, ra2, dec2):
    """Compute the astronomical position angle of position 2 wrt to position 1.

    Inputs are expected in decimal degrees. Returns position angle in radians."""
    ra1rad = ra1 * np.pi / 180.
    ra2rad = ra2 * np.pi / 180.
    dec1rad = dec1 * np.pi / 180.
    dec2rad = dec2 * np.pi / 180.
    radelta = ra2rad - ra1rad
    theta = np.arctan2(np.sin(radelta), np.cos(dec1rad) * np.tan(dec2rad) \
                           - np.sin(dec1rad) * np.cos(radelta))
    return theta


def positionAngleCartesian(ra1, dec1, ra2, dec2):
    ra1rad = ra1 * np.pi / 180.
    ra2rad = ra2 * np.pi / 180.
    dec1rad = dec1 * np.pi / 180.
    dec2rad = dec2 * np.pi / 180.
    radelta = ra2rad - ra1rad
    decdelta = dec2rad - dec1rad
    return np.arctan2(decdelta, radelta)


def greatCircleDistance(ra1, dec1, ra2, dec2):
    """Compute the great circle distance from poistion 1 to position 2.

    Inputs are expected in decimal degrees. Returns great circle distance in decimal degrees.
    Accurate for small separations. This routine uses the Vincenty equation, which also handles
    antipodal point correctly."""
    phi1 = dec1 * np.pi / 180.
    lam1 = ra1 * np.pi / 180.
    phi2 = dec2 * np.pi / 180.
    lam2 = ra2 * np.pi / 180.
    deltalam = np.abs(lam1 - lam2)
    n1 = (np.cos(phi2) * np.sin(deltalam))**2
    n2 = (np.cos(phi1) * np.sin(phi2) -
          np.sin(phi1) * np.cos(phi2) * np.cos(deltalam))**2
    num = np.sqrt(n1 + n2)
    d1 = np.sin(phi1) * np.sin(phi2)
    d2 = np.cos(phi1) * np.cos(phi2) * np.cos(deltalam)
    denom = d1 + d2
    return np.arctan2(num, denom) * 180. / np.pi
