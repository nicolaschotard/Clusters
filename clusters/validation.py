"""Data validation utilisites and plots."""

import numpy as N
import pylab as P
import seaborn
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
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


def stellarLocus(d, mag_type="modelfit_CModel_mag_extcorr",
                 ifilt="i_new", cat='deepCoadd_forced_src'):
    """
    Check colors by plotting stellar loci and comparing with analytical fits.

    First a few color-color (and one mag-color) plots are plotted based on the
    input magnitudes. Since analytical fits are based on SDSS data, the given
    magnitudes are then converted to SDSS mags. Fits are overplotted with the
    derived SDSS magnitudes, and then residuals are calculated and plotted.
    The analytical plots are plotted as an intermediary as well.

    Three plots are saved. Nothing is returned.
    """
    # Get number of filters in table. Need ugriz bands for the plots here.
    filters, nfilters = get_filter_list(d[cat])
    assert 'i' in filters and 'g' in filters and 'r' in filters \
        and 'u' in filters and 'z' in filters, "filter list must contains ugriz"

    # Get the id
    oid = 'objectId' if cat == 'deepCoadd_forced_src' else 'id'

    # Define selection filter
    filt = define_selection_filter(d, cat)

    # Separate the stars from the galaxies
    star, gal = separate_star_gal(d, cat, oid, nfilters, filt=filt)

    # care about only stars for now. get their color magnitudes.
    mrs = star[star['filter'] == 'r'][mag_type]
    mis = star[star['filter'] == 'i'][mag_type]
    mgs = star[star['filter'] == 'g'][mag_type]
    mus = star[star['filter'] == 'u'][mag_type]
    mzs = star[star['filter'] == 'z'][mag_type]

    ########################################################################################
    # Plot color plots corresponding to some in Convey et al. 2007 (arXiv:0707.4473v2)
    # Cut on  the magnitude of the stars
    idx_star = mrs < 23

    nrow, ncol = 4, 2
    size, legendfontsize, labelsize, ticklabelsize = 10, 20, 25, 22
    fig, ax = P.subplots(nrows=nrow, ncols=ncol)
    P.subplots_adjust(hspace=0.4)
    # u-g vs g-r
    r, c = 0, 0
    ax[r, c].scatter(mgs[idx_star] - mrs[idx_star],
                     mus[idx_star] - mgs[idx_star], s=size, color='b')
    ax[r, c].set_title("%s for %d stars" % (mag_type, len(mgs[idx_star])),
                       fontsize=legendfontsize)
    ax[r, c].set_xlabel("g - r", fontsize=labelsize)
    ax[r, c].set_ylabel("u - g", fontsize=labelsize)
    ax[r, c].set_xlim([-0.5, 2.0])
    ax[r, c].set_ylim([-0.5, 5.])
    ax[r, c].tick_params(labelsize=labelsize)

    # g-r vs r-i
    r, c = 0, 1
    ax[r, c].scatter(mrs[idx_star] - mis[idx_star], mgs[idx_star] - mrs[idx_star],
                     s=size, color='b')
    ax[r, c].set_title("%s for %d stars" % (mag_type, len(mgs[idx_star])),
                       fontsize=legendfontsize)
    ax[r, c].set_xlim([-0.5, 3.0])
    ax[r, c].set_ylim([-0.5, 2.0])
    ax[r, c].tick_params(labelsize=labelsize)
    ax[r, c].set_ylabel("g - r", fontsize=ticklabelsize)
    ax[r, c].set_xlabel("r - i", fontsize=ticklabelsize)

    # r-i vs i-z
    r, c = 1, 0
    ax[r, c].scatter(mis[idx_star] - mzs[idx_star], mrs[idx_star] - mis[idx_star],
                     s=size, color='b')
    ax[r, c].set_xlim([-0.5, 2.0])
    ax[r, c].set_ylim([-0.5, 3.0])
    ax[r, c].tick_params(labelsize=labelsize)
    ax[r, c].set_xlabel("i - z", fontsize=ticklabelsize)
    ax[r, c].set_ylabel("r - i", fontsize=ticklabelsize)

    # g-r vs g-i
    r, c = 1, 1
    ax[r, c].scatter(mgs[idx_star] - mis[idx_star], mgs[idx_star] - mrs[idx_star],
                     s=size, color='b')
    ax[r, c].set_xlabel("g - i", fontsize=ticklabelsize)
    ax[r, c].set_ylabel("g - r", fontsize=ticklabelsize)
    ax[r, c].set_xlim([-1., 4.0])
    ax[r, c].set_ylim([-0.5, 2.])
    ax[r, c].tick_params(labelsize=labelsize)

    # u-g vs g-i
    r, c = 2, 0
    ax[r, c].scatter(mgs[idx_star] - mis[idx_star], mus[idx_star] - mgs[idx_star],
                     s=size, color='b')
    ax[r, c].set_xlabel("g - i", fontsize=ticklabelsize)
    ax[r, c].set_ylabel("u - g", fontsize=ticklabelsize)
    ax[r, c].set_xlim([-1., 4.0])
    ax[r, c].set_ylim([-0.5, 5.])
    ax[r, c].tick_params(labelsize=labelsize)

    # r-i vs g-i
    r, c = 2, 1
    ax[r, c].scatter(mgs[idx_star] - mis[idx_star], mrs[idx_star] - mis[idx_star],
                     s=size, color='b')
    ax[r, c].set_xlim([-1., 4.0])
    ax[r, c].set_ylim([-0.5, 2.5])
    ax[r, c].tick_params(labelsize=labelsize)
    ax[r, c].set_xlabel("g - i", fontsize=ticklabelsize)
    ax[r, c].set_ylabel("r - i", fontsize=ticklabelsize)

    # i-z vs g-i
    r, c = 3, 0
    ax[r, c].scatter(mgs[idx_star] - mis[idx_star], mis[idx_star] - mzs[idx_star],
                     s=size, color='b')
    ax[r, c].set_xlim([-1., 4.0])
    ax[r, c].set_ylim([-0.5, 2.5])
    ax[r, c].tick_params(labelsize=labelsize)
    ax[r, c].set_xlabel("g - i", fontsize=ticklabelsize)
    ax[r, c].set_ylabel("i - z", fontsize=ticklabelsize)

    # i vs g-i
    r, c = 3, 1
    ax[r, c].scatter(mgs[idx_star] - mis[idx_star], mis[idx_star], s=size, color='b')
    ax[r, c].set_xlim([-2., 5.0])
    ax[r, c].set_ylim([10., 25.])
    ax[r, c].invert_yaxis()
    ax[r, c].tick_params(labelsize=labelsize)
    ax[r, c].set_xlabel("g - i", fontsize=ticklabelsize)
    ax[r, c].set_ylabel("i", fontsize=ticklabelsize)

    fig.set_size_inches(20, 17)
    filename = 'stellarLocusPlot_' + mag_type + '_' + ifilt + '_' + cat + '.png'
    P.savefig(filename, format='png', bbox_inches='tight')
    print 'Saved plot: ', filename
    P.show()

    ########################################################################################
    ########################################################################################
    # Covert magnitudes to SDSS magnitudes for a direct comparison with Convey et al. 2007.
    # Color corrections CFHT --> SDSS
    def g_sdss_gr(g_mega, r_mega):
        return g_mega + 0.195 * (g_mega - r_mega)

    def r_sdss_gr(r_mega, g_mega):
        return r_mega + 0.011 * (g_mega - r_mega)

    def i_sdss_ri(i_mega, r_mega):
        return i_mega + 0.079 * (r_mega - i_mega)

    def i2_sdss_ri(i2_mega, r_mega):
        return i2_mega + 0.001 * (r_mega - i2_mega)

    def u_sdss_ug(u_mega, g_mega):
        return u_mega + 0.181 * (u_mega - g_mega)

    def g_sdss_gi(g_mega, i_mega):
        return g_mega + 0.103 * (g_mega - i_mega)

    def i_sdss_gi(i_mega, g_mega):
        return i_mega + 0.044 * (g_mega - i_mega)

    def i2_sdss_gi(i2_mega, g_mega):
        return i2_mega - 0.003 * (g_mega - i2_mega)

    def z_sdss_iz(i_mega, z_mega):
        return z_mega - 0.099 * (i_mega - z_mega)

    # covert to SDSS colors
    gsdss = g_sdss_gr(mgs, mrs)
    rsdss = r_sdss_gr(mrs, mgs)
    usdss = u_sdss_ug(mus, mgs)
    zsdss = z_sdss_iz(mis, mzs)
    if ifilt == "i_new":
        print "Will use new i2 (MP9702) filter"
        isdss = i2_sdss_ri(mis, mrs)
        # isdss2 = i2_sdss_gi(mis, mgs)
    else:
        print "Will use old i (MP9701) filter"
        isdss = i_sdss_ri(mis, mrs)
        # isdss2 = i_sdss_gi(mis, mgs)

    ########################################################################################
    # Covey et al. have analytical fits for four color-color plots.
    # Coefficients of 5-degree polynomials in (g-i) + residuals for stellar locus parametrization
    # according to Covey et al. 2007 - arXiv:0707.4473
    # The color reference system is SDSS

    uminusg = N.asarray([1.0636113, -1.6267818, 4.9389572, -3.2809081,
                         0.8725109, -0.0828035, 0.6994, 0.0149])
    gminusr = N.asarray([0.4290263, -0.8852323, 2.0740616, -1.1091553,
                         0.2397461, -0.0183195, 0.2702, 0.0178])
    rminusi = N.asarray([-0.4113500, 1.8229991, -1.9989772, 1.0662075,
                         -0.2284455, 0.0172212, 0.2680, 0.0164])
    iminusz = N.asarray([-0.2270331, 0.7794558, -0.7350749, 0.3727802,
                         -0.0735412, 0.0049808, 0.1261, 0.0071])

    def poly_5(p, x):
        return p[0] + p[1] * x + p[2] * x**2 + p[3] * x**3 + p[4] * x**4 + p[5] * x**5

    # plot these fits to see whats what.
    P.clf()
    model = N.linspace(0.05, 4.4, 80)
    P.plot(model, poly_5(uminusg[:6], model), color='r', label='u-g', lw=2)
    P.plot(model, poly_5(gminusr[:6], model), color='g', label='g-r', lw=2)
    P.plot(model, poly_5(rminusi[:6], model), color='b', label='r-i', lw=2)
    P.plot(model, poly_5(iminusz[:6], model), color='k', label='i-z', lw=2)
    P.title('Covey et al. 2007: Analytical Fits')
    P.xlabel("g - i", fontsize=ticklabelsize)
    P.legend(loc="upper left", fontsize=15)
    fig.set_size_inches(20, 17)
    P.show()

    ################################################################################################
    # Now plot the derived SDSS magnitudes and compare them with the analytical fits.
    # Plot the analog the above plot with sdss mags, now overplotted with
    # Covey et al.'s analytical fit

    # Note: Covey et al.'s fits are a function of (g-i). For dependency on other colors, (g-i)
    # is first fitted wrt to that color and filtered for smoothing. See the first subplot for
    # the documented steps.
    idx_star = (mrs < 23.) & (mgs < 23.)
    nrow, ncol = 4, 2
    size, legendfontsize, labelsize, ticklabelsize = 10, 20, 25, 22
    fig, ax = P.subplots(nrows=nrow, ncols=ncol)
    P.subplots_adjust(hspace=0.4)
    # u-g vs g-r
    r, c = 0, 0
    ax[r, c].scatter(gsdss[idx_star] - rsdss[idx_star], usdss[idx_star] - gsdss[idx_star], s=size,
                     color='b', label="Dervied SDSS mags for %d stars" % (len(mgs[idx_star])))
    model_gminusr = N.linspace(-0.3, 1.7, 80)      # range of g-r we want to plot here
    # analytical expression is valid over 0.05<(g-i)<4.4 only
    tofit_gminusi = N.linspace(0.05, 4.4, 1000)  
    tofit_gminusr = poly_5(gminusr[:6], tofit_gminusi)
    # Fit g-i as a function of g-r
    fit_gminusi_gminusr= interp1d(tofit_gminusr, tofit_gminusi, kind= 'nearest', bounds_error=False)
    # Use the fit and the range of (g-r) to get (u-g) as a function of (g-r).
    # Filter is used for smooth curve; not always necessary but gets rid of kinky fits in some cases
    ax[r, c].plot(model_gminusr, savgol_filter(poly_5(uminusg[:6],
                                                      fit_gminusi_gminusr(model_gminusr)), 7, 1),
                  color='r', label='Covey 2007 model', lw=2)
    ax[r, c].set_xlabel("g - r", fontsize=labelsize)
    ax[r, c].set_ylabel("u - g", fontsize=labelsize)
    ax[r, c].set_xlim([-0.5, 2.0])
    ax[r, c].set_ylim([0., 5.])
    ax[r, c].tick_params(labelsize=labelsize)
    ax[r, c].legend(loc="upper left", fontsize=legendfontsize, labelspacing=0.05)

    # g-r vs r-i
    r, c = 0, 1
    ax[r, c].scatter(rsdss[idx_star] - isdss[idx_star], gsdss[idx_star] - rsdss[idx_star], s=size,
                     color='b', label="Dervied SDSS mags for %d stars" % (len(mgs[idx_star])))
    model_rminusi = N.linspace(-0.3, 2.7, 80)
    tofit_rminusi = poly_5(rminusi[:6], tofit_gminusi)
    fit_gminusi_rminusi = interp1d(tofit_rminusi,tofit_gminusi, kind= 'nearest', bounds_error=False)   
    ax[r, c].plot(model_rminusi, savgol_filter(poly_5(gminusr[:6],
                                                      fit_gminusi_rminusi(model_rminusi)), 7, 1),
                  color='r', label='Covey 2007 model', lw=2)
    ax[r, c].set_xlim([-0.5, 3.0])
    ax[r, c].set_ylim([-0.5, 2.0])
    ax[r, c].tick_params(labelsize=labelsize)
    ax[r, c].set_ylabel("g - r", fontsize=ticklabelsize)
    ax[r, c].set_xlabel("r - i", fontsize=ticklabelsize)
    ax[r, c].legend(loc="lower right", fontsize=legendfontsize, labelspacing=0.05)

    # r-i vs i-z
    r, c = 1, 0
    ax[r, c].scatter(isdss[idx_star] - zsdss[idx_star], rsdss[idx_star] - isdss[idx_star], s=size,
                     color='b', label="Dervied SDSS mags for %d stars" % (len(mgs[idx_star])))
    model_iminusz = N.linspace(-0.3, 1.7, 80)
    tofit_iminusz = poly_5(iminusz[:6], tofit_gminusi)
    fit_gminusi_iminusz = interp1d(tofit_iminusz, tofit_gminusi, kind='nearest', bounds_error=False)
    ax[r, c].plot(model_iminusz, savgol_filter(poly_5(rminusi[:6],
                                                      fit_gminusi_iminusz(model_iminusz)), 7, 1),
                  color='r', label='Covey 2007 model', lw=2)
    ax[r, c].set_xlim([-0.5, 2.0])
    ax[r, c].set_ylim([-0.5, 3.0])
    ax[r, c].tick_params(labelsize=labelsize)
    ax[r, c].set_xlabel("i - z", fontsize=ticklabelsize)
    ax[r, c].set_ylabel("r - i", fontsize=ticklabelsize)
    ax[r, c].legend(loc="upper left", fontsize=legendfontsize, labelspacing=0.05)

    # g-r vs g-i
    r, c = 1, 1
    ax[r, c].scatter(gsdss[idx_star] - isdss[idx_star], gsdss[idx_star] - rsdss[idx_star], s=size,
                     color='b', label="Dervied SDSS mags for %d stars" % (len(mgs[idx_star])))
    model = N.linspace(0.05, 3.8, 80)  # range has to be between 0.05, 4.4
    ax[r, c].plot(model, poly_5(gminusr[:6], model), color='r', label='Covey 2007 model', lw=2)
    ax[r, c].set_xlabel("g - i", fontsize=ticklabelsize)
    ax[r, c].set_ylabel("g - r", fontsize=ticklabelsize)
    ax[r, c].set_xlim([0., 4.0])
    ax[r, c].set_ylim([-0.5, 3.])
    ax[r, c].tick_params(labelsize=labelsize)
    ax[r, c].legend(loc="upper left", fontsize=legendfontsize,labelspacing=0.05)

    # u-g vs g-i
    r, c = 2, 0
    ax[r, c].scatter(gsdss[idx_star] - isdss[idx_star], usdss[idx_star] - gsdss[idx_star], s=size,
                     color='b', label="Dervied SDSS mags for %d stars" % (len(mgs[idx_star])))
    model = N.linspace(0.05, 3.8, 80)  # range has to be between 0.05, 4.4
    ax[r, c].plot(model, poly_5(uminusg[:6], model), color='r', label='Covey 2007 model', lw=2)

    ax[r, c].set_xlabel("g - i", fontsize=ticklabelsize)
    ax[r, c].set_ylabel("u - g", fontsize=ticklabelsize)
    ax[r, c].set_xlim([0., 4.0])
    ax[r, c].set_ylim([0., 5.])
    ax[r, c].tick_params(labelsize=labelsize)
    ax[r, c].legend(loc="upper left", fontsize=legendfontsize, labelspacing=0.05,)

    # r-i vs g-i
    r, c = 2, 1
    ax[r, c].scatter(gsdss[idx_star] - isdss[idx_star], rsdss[idx_star] - isdss[idx_star], s=size,
                     color='b', label="Dervied SDSS mags for %d stars" % (len(mgs[idx_star])))
    model = N.linspace(0.1, 3.8, 80)    # range has to be between 0.05, 4.4
    ax[r, c].plot(model, poly_5(rminusi[:6], model), color='r', label='Covey 2007 model', lw=2)
    ax[r, c].set_xlim([-0.5, 4.0])
    ax[r, c].set_ylim([-0.5, 2.5])
    ax[r, c].tick_params(labelsize=labelsize)
    ax[r, c].set_xlabel("g - i", fontsize=ticklabelsize)
    ax[r, c].set_ylabel("r - i", fontsize=ticklabelsize)
    ax[r, c].legend(loc="upper left", fontsize=legendfontsize, labelspacing=0.05,)

    # i-z vs g-i
    r, c = 3, 0
    ax[r, c].scatter(gsdss[idx_star] - isdss[idx_star], isdss[idx_star] - zsdss[idx_star], s=size,
                     color='b', label="Dervied SDSS mags for %d stars" % (len(mgs[idx_star])))
    model = N.linspace(0.1, 3.8, 80)    # range has to be between 0.05, 4.4
    ax[r, c].plot(model, poly_5(iminusz[:6], model), color='r', label='Covey 2007 model', lw=2)
    ax[r, c].set_xlim([-0.5, 4.0])
    ax[r, c].set_ylim([-0.5, 2.5])
    ax[r, c].tick_params(labelsize=labelsize)
    ax[r, c].set_xlabel("g - i", fontsize=ticklabelsize)
    ax[r, c].set_ylabel("i - z", fontsize=ticklabelsize)
    ax[r, c].legend(loc="upper left", fontsize= legendfontsize, labelspacing= 0.05,)

    # i vs g-i
    r, c = 3, 1
    ax[r, c].scatter(gsdss[idx_star] - isdss[idx_star], isdss[idx_star], s=size, color='b',
                     label="Dervied SDSS mags for %d stars" % (len(mgs[idx_star])))
    ax[r, c].set_xlim([-2., 5.0])
    ax[r, c].set_ylim([10., 25.])
    ax[r, c].invert_yaxis()
    ax[r, c].tick_params(labelsize=labelsize)
    ax[r, c].set_xlabel("g - i", fontsize=ticklabelsize)
    ax[r, c].set_ylabel("i", fontsize=ticklabelsize)
    ax[r, c].legend(loc="upper left", fontsize=legendfontsize, labelspacing=0.05)

    fig.set_size_inches(20, 17)
    filename = 'derivedSDSSstellarLocusPlot_' + mag_type + \
               '_' + ifilt + '_' + cat + '.png'
    P.savefig(filename, format='png', bbox_inches='tight')
    print 'Saved plot: ', filename
    P.show()

    ################################################################################################
    # Plot the residuals; code follows the calcualtion in stellarLocus.
    # For the cases where a direct analytical expression is not available,
    # the fits calcualted above are used.
    nrow, ncol = 4, 2
    size, legendfontsize, labelsize, ticklabelsize = 10, 20, 25, 22
    fig, ax = P.subplots(nrows=nrow, ncols=ncol)
    P.subplots_adjust(hspace=0.4)

    numstep = 40
    # u-g vs g-r
    r, c = 0, 0
    xmin, xmax = -0.3, 2.0
    colval, step = N.linspace(xmin, xmax, num=numstep, retstep=True)
    res, err = [], []
    for val in colval:
        idc = ((gsdss[idx_star] - rsdss[idx_star]) > (val - step / 2)) & \
              ((gsdss[idx_star] - rsdss[idx_star]) < (val + step / 2))
        res.append(poly_5(uminusg[:6], fit_gminusi_gminusr(val)) -
                   N.mean((usdss[idx_star] - gsdss[idx_star])[idc]))
        err.append(N.std((usdss[idx_star] - gsdss[idx_star])[idc]) /
                   N.sqrt(len((usdss[idx_star] - gsdss[idx_star])[idc])))

    ax[r, c].errorbar(colval, res, yerr=err, ls='None', fmt='s',
                      capsize=25, color='b', lw=1)
    ax[r, c].set_title('%s - Residuals' % mag_type, fontsize=legendfontsize)
    ax[r, c].set_xlim([xmin - 0.3, xmax + 0.3])
    ax[r, c].set_xlabel("g - r", fontsize=ticklabelsize)
    ax[r, c].set_ylabel("Model - data (u - g)", fontsize=ticklabelsize)
    ax[r, c].tick_params(labelsize=labelsize)

    # g-r vs r-i
    r, c = 0, 1
    xmin, xmax = -0.3, 2.0
    colval, step = N.linspace(xmin, xmax, num=numstep, retstep=True)
    res, err = [], []
    for val in colval:
        idc = ((rsdss[idx_star] - isdss[idx_star]) > (val - step / 2)) & \
              ((rsdss[idx_star] - isdss[idx_star]) < (val + step / 2))
        res.append(poly_5(gminusr[:6], fit_gminusi_rminusi(val)) - N.mean((gsdss[idx_star] -
                                                                           rsdss[idx_star])[idc]))
        err.append(N.std((gsdss[idx_star] - rsdss[idx_star])[idc]) /
                   N.sqrt(len((gsdss[idx_star] - rsdss[idx_star])[idc])))

    ax[r, c].errorbar(colval, res, yerr=err, ls='None', fmt='s',
                      capsize=25, color='b', lw=1)
    ax[r, c].set_title('%s - Residuals' % mag_type, fontsize=legendfontsize)
    ax[r, c].set_xlim([xmin - 0.3, xmax + 0.3])
    ax[r, c].set_xlabel("r - i", fontsize=ticklabelsize)
    ax[r, c].set_ylabel("Model - data (g - r)", fontsize=ticklabelsize)
    ax[r, c].tick_params(labelsize=labelsize)

    # r-i vs i-z
    r, c = 1, 0
    xmin, xmax = -0.3, 2.0
    colval, step = N.linspace(xmin, xmax, num=numstep, retstep=True)
    res, err = [], []
    for val in colval:
        idc = ((isdss[idx_star] - zsdss[idx_star]) > (val - step / 2)) & \
              ((isdss[idx_star] - zsdss[idx_star]) < (val + step / 2))
        res.append(poly_5(rminusi[:6], fit_gminusi_iminusz(val)) -
                   N.mean((gsdss[idx_star] - rsdss[idx_star])[idc]))
        err.append(N.std((rsdss[idx_star] - isdss[idx_star])[idc]) /
                   N.sqrt(len((rsdss[idx_star] - isdss[idx_star])[idc])))

    ax[r, c].errorbar(colval, res, yerr=err, ls='None', fmt='s',
                      capsize=25, color='b', lw=1)
    ax[r, c].set_xlim([xmin - 0.2, xmax + 0.2])
    ax[r, c].set_xlabel("i - z", fontsize=ticklabelsize)
    ax[r, c].set_ylabel("Model - data (r - i)", fontsize=ticklabelsize)
    ax[r, c].tick_params(labelsize=labelsize)

    # g-r vs g-i
    r, c = 1, 1
    xmin, xmax = 0.05, 4.0
    colval, step = N.linspace(xmin, xmax, num=numstep, retstep=True)
    res, err = [], []

    for val in colval:
        idc = ((gsdss[idx_star] - isdss[idx_star]) > (val - step / 2)) & \
              ((gsdss[idx_star] - isdss[idx_star]) < (val + step / 2))
        res.append(poly_5(gminusr[:6], val) - N.mean((gsdss[idx_star] -
                                                      rsdss[idx_star])[idc]))
        err.append(N.std((gsdss[idx_star] - rsdss[idx_star])[idc]) /
                   N.sqrt(len((gsdss[idx_star] - rsdss[idx_star])[idc])))
    ax[r, c].errorbar(colval, res, yerr=err, ls='None', fmt='s', capsize=25,
                      color='b', lw=1, label='%s - Residuals' % mag_type)
    ax[r, c].set_xlim([xmin - 0.3, xmax + 0.3])
    ax[r, c].set_ylim([-0.1, 0.1])
    ax[r, c].set_xlabel("g - i", fontsize=ticklabelsize)
    ax[r, c].set_ylabel("Model - data (g - r)", fontsize=ticklabelsize)
    ax[r, c].tick_params(labelsize=labelsize)

    # u-g vs g-i
    r, c = 2, 0
    xmin, xmax = 0.05, 4.0
    colval, step = N.linspace(xmin, xmax, num=numstep, retstep=True)
    res, err = [], []
    for val in colval:
        idc = ((gsdss[idx_star] - isdss[idx_star]) > (val - step / 2)) & \
              ((gsdss[idx_star] - isdss[idx_star]) < (val + step / 2))
        res.append(poly_5(uminusg[:6], val) - N.mean((usdss[idx_star] -
                                                      gsdss[idx_star])[idc]))
        err.append(N.std((usdss[idx_star] - gsdss[idx_star])[idc]) /
                   N.sqrt(len((usdss[idx_star] - gsdss[idx_star])[idc])))
    ax[r, c].errorbar(colval, res, yerr=err, ls='None', fmt='s',
                      capsize=25, color='b', lw=1)
    ax[r, c].set_xlim([xmin - 0.3, xmax + 0.3])
    ax[r, c].set_xlabel("g - i", fontsize=ticklabelsize)
    ax[r, c].set_ylabel("Model - data (u - g)", fontsize=ticklabelsize)
    ax[r, c].tick_params(labelsize=labelsize)

    # r-i vs g-i
    r, c = 2, 1
    xmin, xmax = 0.05, 4.0
    colval, step = N.linspace(xmin, xmax, num=numstep, retstep=True)
    res, err = [], []
    for val in colval:
        idc = ((gsdss[idx_star] - isdss[idx_star]) > (val - step / 2)) & \
              ((gsdss[idx_star] - isdss[idx_star]) < (val + step / 2))
        res.append(poly_5(rminusi[:6], val) - N.mean((rsdss[idx_star] -
                                                      isdss[idx_star])[idc]))
        err.append(N.std((rsdss[idx_star] - isdss[idx_star])[idc]) /
                   N.sqrt(len((rsdss[idx_star] - isdss[idx_star])[idc])))

    ax[r, c].errorbar(colval, res, yerr=err, ls='None', fmt='s',
                      capsize=25, color='b', lw=1)
    ax[r, c].set_xlim([xmin - 0.3, xmax + 0.3])
    ax[r, c].set_xlabel("g - i", fontsize=ticklabelsize)
    ax[r, c].set_ylabel("Model - data (r - i)", fontsize=ticklabelsize)
    ax[r, c].tick_params(labelsize=labelsize)

    # i-z vs g-i
    r, c = 3, 0
    xmin, xmax = 0.05, 4.0
    colval, step = N.linspace(xmin, xmax, num=numstep, retstep=True)
    res, err = [], []
    for val in colval:
        idc = ((gsdss[idx_star] - isdss[idx_star]) > (val - step / 2)) & \
              ((gsdss[idx_star] - isdss[idx_star]) < (val + step / 2))
        res.append(poly_5(iminusz[:6], val) - N.mean((isdss[idx_star] -
                                                      zsdss[idx_star])[idc]))
        err.append(N.std((isdss[idx_star] - zsdss[idx_star])[idc]) /
                   N.sqrt(len((isdss[idx_star] - zsdss[idx_star])[idc])))

    ax[r, c].errorbar(colval, res, yerr=err, ls='None', fmt='s',
                      capsize=25, color='b', lw=1)
    ax[r, c].set_xlim([xmin - 0.3, xmax + 0.3])
    ax[r, c].set_xlabel("g - i", fontsize=ticklabelsize)
    ax[r, c].set_ylabel("Model - data (i - z)", fontsize=ticklabelsize)
    ax[r, c].tick_params(labelsize=labelsize)

    ax[3, 1].axis('off')
    fig.set_size_inches(20, 17)
    filename = 'residualsDerivedSDSSvsAnalyticalFit_' + mag_type + '_' + \
               ifilt + '_' + cat + '.png'
    P.savefig(filename, format='png', bbox_inches='tight')
    print 'Saved plot: ', filename
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

    # & (d[cat]['filter'] == 'i')

    # Separate the stars from the galaxies
    star, gal = separate_star_gal(d, cat, oid, nfilters, filt=filt)
    # filti = d[cat]['filter'] == 'i'
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
    ax2.set_xlabel('Radius in pixels', fontsize=10)
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
