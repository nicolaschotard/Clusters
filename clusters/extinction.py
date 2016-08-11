"""
From http://argonaut.skymaps.info/usage#function-call
"""

import json
import requests
import pylab as P
import numpy as N
import seaborn

def query(lon, lat, coordsys='gal', mode='full', limit=500000):
    '''
    Send a line-of-sight reddening query to the Argonaut web server.

    lon, lat: longitude and latitude, in degrees.
    coordsys: 'gal' for Galactic, 'equ' for Equatorial (J2000).
    mode: 'full', 'lite' or 'sfd'

    In 'full' mode, outputs a dictionary containing, among other things:

    'distmod':    The distance moduli that define the distance bins.
    'best':       The best-fit (maximum proability density)
                  line-of-sight reddening, in units of SFD-equivalent
                  E(B-V), to each distance modulus in 'distmod.' See
                  Schlafly & Finkbeiner (2011) for a definition of the
                  reddening vector (use R_V = 3.1).
    'samples':    Samples of the line-of-sight reddening, drawn from
                  the probability density on reddening profiles.
    'success':    1 if the query succeeded, and 0 otherwise.
    'converged':  1 if the line-of-sight reddening fit converged, and
                  0 otherwise.
    'n_stars':    # of stars used to fit the line-of-sight reddening.
    'DM_reliable_min':  Minimum reliable distance modulus in pixel.
    'DM_reliable_max':  Maximum reliable distance modulus in pixel.

    Less information is returned in 'lite' mode, while in 'sfd' mode,
    the Schlegel, Finkbeiner & Davis (1998) E(B-V) is returned.
    '''

    # Make sure to have less than 500000 objects (the limit).
    # Cut the list in smaller pieces if that is the case.
    def chunk(ilist, length):
        """Divide a list 'l' into smaller lists of maximal length 'num'."""
        return [ilist[i:i + length] for i in range(0, len(ilist), length)]
    
    if len(lon) >= limit:
        lons = chunk(lon, limit - 1)
        lats = chunk(lat, limit - 1)
        dicts = [query(loni, lati, coordsys=coordsys, mode=mode) for loni, lati in zip(lons, lats)]
        for dic in dicts[1:]:
            for k in dic:
                dicts[0][k].extend(dic[k])
        return dicts[0]

    url = 'http://argonaut.skymaps.info/gal-lb-query-light'

    payload = {'mode': mode}

    if coordsys.lower() in ['gal', 'g']:
        payload['l'] = lon
        payload['b'] = lat
    elif coordsys.lower() in ['equ', 'e']:
        payload['ra'] = lon
        payload['dec'] = lat
    else:
        raise ValueError("coordsys '{0}' not understood.".format(coordsys))

    req = requests.post(url, data=json.dumps(payload), headers={'content-type': 'application/json'})

    try:
        req.raise_for_status()
    except requests.exceptions.HTTPError as excep:
        print 'Response received from Argonaut:'
        print req.text
        raise excep

    return json.loads(req.text)

def from_ebv_sfd_TO_sdss_albd(ebv):
    """Return A(lbd) for the 5 SDSS filters: u, g, r, i, z"""
    coeff = {'u': 5.155, 'g': 3.793, 'r': 2.751, 'i': 2.086, 'z': 1.479}
    return {f: coeff[f] * N.array(ebv) for f in coeff}

def from_sdss_albd_TO_megacam_albd(sdss):
    """Return A(lbd) for the 6 Megecam filters: u, g, r, i_old, i_new, z"""
    megacam = {}
    megacam['u'] = sdss['u'] - 0.241 * (sdss['u'] - sdss['g'])
    megacam['g'] = sdss['g'] - 0.153 * (sdss['g'] - sdss['r'])
    megacam['r'] = sdss['r'] - 0.024 * (sdss['g'] - sdss['r'])
    megacam['z'] = sdss['z'] - 0.074 * (sdss['i'] - sdss['z'])
    megacam['i_old'] = sdss['i'] - 0.085 * (sdss['r'] - sdss['i'])
    megacam['i_new'] = sdss['i'] - 0.003 * (sdss['r'] - sdss['i'])
    return megacam

def from_ebv_sfd_TO_megacam_albd(ebv):
    """Return A(lbd) for the 6 Megacam filters: u, g, r, i, z"""
    return from_sdss_albd_TO_megacam_albd(from_ebv_sfd_TO_sdss_albd(ebv))

def plots(ra, dec, ebv, albd, title=None, figname="", filters=None):

    fig = P.figure()
    ax = fig.add_subplot(111, xlabel='RA (deg)', ylabel='DEC (deg)')
    scat = ax.scatter(ra, dec, c=ebv, cmap=(P.cm.jet))
    cbar = fig.colorbar(scat)
    cbar.set_label('E(B-V)')
    if title is not None:
        ax.set_title(title)
    fig.savefig(figname+"_ebmv_map.png")

    if filters is None:
        filters = albd.keys()
    fig = P.figure()
    ax = fig.add_subplot(111, xlabel='A(lbd)', ylabel='#')
    for filt in filters:
        ax.hist(albd[filt], histtype='step', lw=2, label='<%s>=%.2f' % (filt, N.mean(albd[filt])))
    if title is not None:
        ax.set_title(title)
    ax.legend(loc='best')
    fig.savefig(figname+"_albd.png")

    P.show()
