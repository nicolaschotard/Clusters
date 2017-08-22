"""Converting tools for extinction."""

import pylab as P
import numpy as N


def from_ebv_sfd_to_sdss_albd(ebv):
    """Return A(lbd) for the 5 SDSS filters: u, g, r, i, z."""
    coeff = {'u': 5.155, 'g': 3.793, 'r': 2.751, 'i': 2.086, 'z': 1.479}
    return {f: coeff[f] * N.array(ebv) for f in coeff}


def from_sdss_albd_to_megacam_albd(sdss):
    """Return A(lbd) for the 6 Megecam filters: u, g, r, i_old, i_new, z."""
    megacam = {}
    megacam['u'] = sdss['u'] - 0.241 * (sdss['u'] - sdss['g'])
    megacam['g'] = sdss['g'] - 0.153 * (sdss['g'] - sdss['r'])
    megacam['r'] = sdss['r'] - 0.024 * (sdss['g'] - sdss['r'])
    megacam['z'] = sdss['z'] - 0.074 * (sdss['i'] - sdss['z'])
    megacam['i_old'] = sdss['i'] - 0.085 * (sdss['r'] - sdss['i'])
    megacam['i_new'] = sdss['i'] - 0.003 * (sdss['r'] - sdss['i'])
    return megacam


def from_ebv_sfd_to_megacam_albd(ebv):
    """Return A(lbd) for the 6 Megacam filters: u, g, r, i, z."""
    return from_sdss_albd_to_megacam_albd(from_ebv_sfd_to_sdss_albd(ebv))


def plots(ra, dec, ebv, albd, title=None, figname=""):
    """Plot the extinction sky-map."""
    fig = P.figure()
    ax = fig.add_subplot(111, xlabel='RA (deg)', ylabel='DEC (deg)')
    scat = ax.scatter(ra, dec, c=ebv, cmap=(P.cm.jet), edgecolor='none')
    cbar = fig.colorbar(scat)
    cbar.set_label('E(B-V)')
    if title is not None:
        ax.set_title(title)
    fig.savefig(figname + "_ebmv_map.png")

    fig = P.figure()
    ax = fig.add_subplot(111, xlabel='A(lbd)', ylabel='#')
    ax.hist(albd, histtype='step', lw=2, label='<>=%.2f' % N.mean(albd))
    if title is not None:
        ax.set_title(title)
    ax.legend(loc='best')
    fig.savefig(figname + "_albd.png")

    P.show()
