"""Data validation utilisites and plots."""

import numpy as N
import pylab as P
import seaborn
from . import data


def load_cluster(cluster="MACSJ2243.3-0935", ifilt="i_new"):
    """Load the data for a given cluster."""
    # load astropy tables from hdf5 file
    d = data.read_hdf5(cluster + "_all.hdf5")
    # read extinction law parameters
    d2 = data.read_hdf5(cluster + "_all_extinction.hdf5", path="extinction", dic=False)

    # correct maggnitude for extinction
    data.correct_for_extinction(d['deepCoadd_forced_src'], d2, ifilt=ifilt)
    data.correct_for_extinction(d['deepCoadd_meas'], d2, ifilt=ifilt)

    return d


def get_filter_list(table):
    """Get the filter list and number of filter in a table."""
    filters = set(table['filter'])
    return filters, len(filters)


def define_selection_filter(d, cat):
    """Define and return a standard quality selection filter."""
    filt = d['deepCoadd_meas']['base_ClassificationExtendedness_flag'] == 0
    filt &= d['deepCoadd_meas']['detect_isPrimary'] == 1
    filt &= d[cat]['modelfit_CModel_flag'] == 0
    filt &= d[cat]['modelfit_CModel_flux'] > 0

    filt &= (d[cat]['modelfit_CModel_flux'] /
             d[cat]['modelfit_CModel_fluxSigma']) > 10

    return filt


def separate_star_gal(d, cat, oid, nfilters, filt=None):
    """Return two clean tables: one for the stars, the other for the galaxies."""
    filt_star = d['deepCoadd_meas']['base_ClassificationExtendedness_value'] < 0.5
    filt_gal = d['deepCoadd_meas']['base_ClassificationExtendedness_value'] > 0.5

    if filt is not None:
        filt_star &= filt
        filt_gal &= filt

    # Select stars
    star = d[cat][filt_star].group_by(oid)
    stars = star.groups[(star.groups.indices[1:] -
                         star.groups.indices[:-1]) == nfilters]

    # Select galaxies
    gal = d[cat][filt_gal].group_by(oid)
    galaxies = gal.groups[(gal.groups.indices[1:] -
                           gal.groups.indices[:-1]) == nfilters]

    return stars, galaxies


def stellarLocus(d, mag_type="modelfit_CModel_mag_extcorr", ifilt="i_new", cat='deepCoadd_forced_src'):
    """Plot stellar locus."""
    # Get number of filters in table.
    # For stellar locus we need at least the g, r and i filters
    filters, nfilters = get_filter_list(d['deepCoadd_meas'])
    assert 'i' in filters and 'g' in filters and 'r' in filters, \
        "filter list must contains gri"

    # Get the id
    oid = 'objectId' if cat == 'deepCoadd_forced_src' else 'id'

    # Define selection filter
    filt = define_selection_filter(d, cat)

    # Seprat ethe stars from the galaxies
    star, gal = separate_star_gal(d, cat, oid, nfilters, filt=filt)

    # get color magnitudes for stars and galaxies
    mrs = star[star['filter'] == 'r'][mag_type]
    mis = star[star['filter'] == 'i'][mag_type]
    mgs = star[star['filter'] == 'g'][mag_type]
    mrg = gal[gal['filter'] == 'r'][mag_type]
    mig = gal[gal['filter'] == 'i'][mag_type]
    mgg = gal[gal['filter'] == 'g'][mag_type]

    # Plot r-i Vs g-r for stars and galaxies

    # cut on magnitude
    idx_star = mrs < 23
    idx_gal = mrg < 22

    fig, (ax1) = P.subplots(ncols=1)
    ax1.scatter(mgs[idx_star] - mrs[idx_star], mrs[idx_star] - mis[idx_star], s=1,
                color='b', label="stars %d" % len(mgs[idx_star]))
    ax1.scatter(mgg[idx_gal] - mrg[idx_gal], mrg[idx_gal] - mig[idx_gal], s=1,
                color='r', label="galaxies %d" % len(mgg[idx_gal]))
    ax1.set_xlim([-0.5, 2.0])
    ax1.set_ylim([-0.5, 2.5])
    ax1.tick_params(labelsize=20)
    ax1.set_xlabel("g - r", fontsize=20)
    ax1.set_ylabel("r - i", fontsize=20)
    ax1.legend(loc="upper left")

    # Coefficients of 5-degree polynomials in (g-i) + residuals for stellar locus parametrization
    # according to Covey et al. 2007 - arXiv:0707.4473
    # The color reference system is SDSS

    #uMinusg = N.asarray([1.0636113, -1.6267818, 4.9389572, -3.2809081,
    #                     0.8725109, -0.0828035, 0.6994, 0.0149])
    gminusr = N.asarray([0.4290263, -0.8852323, 2.0740616, -1.1091553,
                         0.2397461, -0.0183195, 0.2702, 0.0178])
    rminusi = N.asarray([-0.4113500, 1.8229991, -1.9989772, 1.0662075,
                         -0.2284455, 0.0172212, 0.2680, 0.0164])
    #iMinusz = N.asarray([-0.2270331, 0.7794558, -0.7350749, 0.3727802,
    #                     -0.0735412, 0.0049808, 0.1261, 0.0071])

    def poly_5(p, x):
        return p[0] + p[1] * x + p[2] * x ** 2 + p[3] * x ** 3 + p[4] * x ** 4 + p[5] * x ** 5

    # Color corrections CFHT --> SDSS
    def g_sdss_gr(g_mega, r_mega):
        return g_mega + 0.195 * (g_mega - r_mega)

    def r_sdss_gr(r_mega, g_mega):
        return r_mega + 0.011 * (g_mega - r_mega)

    def i_sdss_ri(i_mega, r_mega):
        return i_mega + 0.079 * (r_mega - i_mega)

    def i2_sdss_ri(i2_mega, r_mega):
        return i2_mega + 0.001 * (r_mega - i2_mega)

    # u_sdss_ug = lambda u_mega, g_mega: u_mega + 0.181 * (u_mega - g_mega)
    # g_sdss_gi = lambda g_mega, i_mega: g_mega + 0.103 * (g_mega - i_mega)
    # i_sdss_gi = lambda i_mega, g_mega: i_mega + 0.044 * (g_mega - i_mega)
    # i2_sdss_gi = lambda i2_mega, g_mega: i2_mega - 0.003 * (g_mega - i2_mega)
    # z_sdss_iz = lambda z_mega, i_mega: z_mega - 0.099 * (i_mega - z_mega)

    gsdss = g_sdss_gr(mgs, mrs)
    rsdss = r_sdss_gr(mrs, mgs)
    if ifilt == "i_new":
        print "Will use new i2 (MP9702) filter"
        isdss = i2_sdss_ri(mis, mrs)
        # isdss2 = i2_sdss_gi(mis, mgs)
    else:
        print "Will use old i (MP9701) filter"
        isdss = i_sdss_ri(mis, mrs)
        # isdss2 = i_sdss_gi(mis, mgs)

    # plot stellar locus for stars and compare to Covey et al. 2007 model

    idx_star = (mrs < 23) & (mgs < 23)
    fig, (ax1, ax2) = P.subplots(ncols=2)

    ax1.scatter(gsdss[idx_star] - isdss[idx_star], rsdss[idx_star] - isdss[idx_star],
                s=1, color='b', label='%s ' % mag_type)
    ax1.set_xlim([-0.5, 4.])
    ax1.set_ylim([-0.3, 2.5])

    model = N.linspace(0.4, 3.8, 80)
    ax1.plot(model, poly_5(rminusi[:6], model), color='r', label='Covey 2007 model', lw=2)
    ax1.set_xlabel("g - i", fontsize=10)
    ax1.set_ylabel("r - i", fontsize=10)
    ax1.tick_params(labelsize=10)
    ax1.legend(loc="upper left", fontsize=10)

    ax2.scatter(gsdss[idx_star] - isdss[idx_star], gsdss[idx_star] - rsdss[idx_star],
                s=1, color='b', label='%s' % mag_type)
    ax2.set_xlim([-0.5, 4.])
    ax2.set_ylim([-0.3, 2.0])
    ax2.plot(model, poly_5(gminusr[:6], model), color='r', label='Covey 2007 model', lw=2)
    ax2.set_xlabel("g - i", fontsize=10)
    ax2.set_ylabel("g - r", fontsize=10)
    ax2.tick_params(labelsize=10)
    ax2.legend(loc="upper left", fontsize=10)

    # Plot residuals
    fig, (ax1, ax2) = P.subplots(ncols=2)

    xmin = 0.5
    xmax = 3.0
    numstep = 40
    res = []
    err = []
    colval, step = N.linspace(xmin, xmax, num=numstep, retstep=True)
    for c in colval:
        idc = ((gsdss[idx_star] - isdss[idx_star]) > (c - step / 2)) & \
              ((gsdss[idx_star] - isdss[idx_star]) < (c + step / 2))
        res.append(poly_5(rminusi[:6], c) - N.mean((rsdss[idx_star] -
                                                    isdss[idx_star])[idc]))
        err.append(N.std((rsdss[idx_star] - isdss[idx_star])[idc]) /
                   N.sqrt(len((rsdss[idx_star] - isdss[idx_star])[idc])))
    ax1.errorbar(colval, res, yerr=err, ls='None', fmt='s',
                 capsize=25, color='b', lw=1, label='%s - Residuals' % mag_type)
    ax1.set_xlim([xmin - 0.2, xmax + 0.2])
    ax1.set_ylim([-0.1, 0.1])
    ax1.set_xlabel("g - i", fontsize=10)
    ax1.set_ylabel("Model - data (r-i)", fontsize=10)
    ax1.tick_params(labelsize=10)
    ax1.legend(loc="upper left", fontsize=10)

    res = []
    err = []
    for c in colval:
        idc = ((gsdss[idx_star] - isdss[idx_star]) > (c - step / 2)) & \
              ((gsdss[idx_star] - isdss[idx_star]) < (c + step / 2))
        res.append(poly_5(gminusr[:6], c) - N.mean((gsdss[idx_star] -
                                                    rsdss[idx_star])[idc]))
        err.append(N.std((gsdss[idx_star] - rsdss[idx_star])[idc]) /
                   N.sqrt(len((gsdss[idx_star] - rsdss[idx_star])[idc])))
    ax2.errorbar(colval, res, yerr=err, ls='None', fmt='s', capsize=25,
                 color='b', lw=1, label='%s - Residuals' % mag_type)
    ax2.set_xlim([xmin - 0.2, xmax + 0.2])
    ax2.set_ylim([-0.1, 0.1])
    ax2.set_xlabel("g - i", fontsize=10)
    ax2.set_ylabel("Model - data (g-r)", fontsize=10)
    ax2.tick_params(labelsize=10)
    ax2.legend(loc="upper left", fontsize=10)
    P.tight_layout()
    P.show()


def compute_elipticities(xx, yy, xy):
    """Compute star elipticities from second momments."""
    denom = xx + 2. * N.sqrt(xx * yy - N.square(xy))
    e1 = (xx - yy) / denom
    e2 = 2.0 * xy / denom
    return e1, e2


def check_star_elipticities(d, cat='deepCoadd_meas', oid='id'):
    """
    Compute star elipticities from second momments and check if psf correction is valid.

    Also check magnitude vss radius
    """
    filters, nfilters = get_filter_list(d[cat])
    assert 'i' in filters, "'i' filter must be in the list of filters"

    # Define selection filter
    filt = define_selection_filter(d, cat)

    #& (d[cat]['filter'] == 'i')

    # Separate the stars from the galaxies
    star, gal = separate_star_gal(d, cat, oid, nfilters, filt=filt)
    #filti = d[cat]['filter'] == 'i'
    star = star[star['filter'] == 'i']
    gal = gal[gal['filter'] == 'i']

    # Get the needed variables
    moments = {'star': {'xx': star['ext_shapeHSM_HsmSourceMoments_xx'],
                        'yy': star['ext_shapeHSM_HsmSourceMoments_yy'],
                        'xy': star['ext_shapeHSM_HsmSourceMoments_xy']},
               'psfs': {'xx': star['ext_shapeHSM_HsmPsfMoments_xx'],
                        'yy': star['ext_shapeHSM_HsmPsfMoments_yy'],
                        'xy': star['ext_shapeHSM_HsmPsfMoments_xy']},
               'gal': {'xx': gal['ext_shapeHSM_HsmSourceMoments_xx'],
                       'yy': gal['ext_shapeHSM_HsmSourceMoments_yy'],
                       'xy': gal['ext_shapeHSM_HsmSourceMoments_xy']}}
    magi = {'star': star['modelfit_CModel_mag'],
            'gal': gal['modelfit_CModel_mag']}
    radius = {'star': N.sqrt(moments['star']['xx'] + moments['star']['yy']),
              'gal': N.sqrt(moments['gal']['xx'] + moments['gal']['yy'])}

    # Plot magnitude as a function of the source radius computed from second momments
    fig, (ax1, ax2) = P.subplots(ncols=2)
    ax1.scatter(radius['star'], magi['star'], s=1, color='b',
                label='Stars %d' % len(magi['star']))
    ax1.scatter(radius['gal'], magi['gal'], s=1, color='r',
                label='Galaxies %d' % len(magi['gal']))
    ax1.set_xlim([1., 6.])
    ax1.set_ylim([16, 26])
    ax1.set_xlabel('Radius in pixels', fontsize=10)
    ax1.set_ylabel('Magnitude i', fontsize=10)
    ax1.tick_params(labelsize=10)
    ax1.legend(loc="lower left", fontsize=10)

    ax2.hist(radius['star'][magi['star'] < 23], bins=80, range=[1.7, 3], color='b')
    ax2.hist(radius['gal'][magi['gal'] < 23], bins=80, range=[1.7, 3], color='r')
    ax1.set_xlabel('Radius in pixels', fontsize=10)

    P.tight_layout()

    e1source, e2source = compute_elipticities(moments['star']['xx'],
                                              moments['star']['yy'],
                                              moments['star']['xy'])
    e1psf, e2psf = compute_elipticities(moments['psfs']['xx'],
                                        moments['psfs']['yy'],
                                        moments['psfs']['xy'])

    idx = magi['star'] < 22.5
    fig, (ax0, ax1) = P.subplots(ncols=2)
    ax0.scatter(e1source[idx], e2source[idx], s=1, color='b')
    ax0.set_xlabel('e1(source)', fontsize=10)
    ax0.set_ylabel('e2(source)', fontsize=10)
    ax0.set_xlim([-0.10, 0.10])
    ax0.set_ylim([-0.10, 0.10])
    ax1.tick_params(labelsize=10)
    ax1.scatter(e1source[idx] - e1psf[idx],
                e2source[idx] - e2psf[idx], s=1, color='b')
    ax1.set_xlim([-0.10, 0.10])
    ax1.set_ylim([-0.10, 0.10])
    ax1.set_xlabel('e1(source) - e1(psf)', fontsize=10)
    ax1.set_ylabel('e2(source) - e2(psf)', fontsize=10)
    ax1.tick_params(labelsize=10)

    P.tight_layout()

    P.show()
