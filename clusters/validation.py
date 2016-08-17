import numpy as N
import pylab as P
import seaborn
from Clusters import data

"""
Data validation utilisites and plots
"""

def load_cluster(cluster="MACSJ2243.3-0935", ifilt="i_new"):
    """Load the data for a given cluster."""
    # load astropy tables from hdf5 file
    d = data.read_data(cluster + "_all.hdf5")
    # read extinction law parameters
    d2 = data.read_data(cluster + "_all_extinction.hdf5", path="extinction")

    # correct maggnitude for extinction
    data.correct_for_extinction(d["forced"], d2, ifilt=ifilt)
    data.correct_for_extinction(d["meas"], d2, ifilt=ifilt)

    return d

def stellarLocus(d, mag_type="modelfit_CModel_mag_extcorr", ifilt="i_new", cat="forced"):
    """Plot stellar locus."""
    # get number of filters in table.
    # For stellar locus we need at least the g, r ani filters
    nfilters = len(d['meas'].group_by('id').groups[0])

    if cat == 'forced':
        oid = 'objectId'
    else:
        oid = 'id'

    # define selection filter
    filt = d['meas']['base_ClassificationExtendedness_flag'] == 0
    filt &= d['meas']['detect_isPrimary'] == 1
    filt &= d[cat]['modelfit_CModel_flag'] == 0
    filt &= d[cat]['modelfit_CModel_flux'] > 0

    filt &= (d[cat]['modelfit_CModel_flux'] /
             d[cat]['modelfit_CModel_fluxSigma']) > 10

    filtS = d['meas']['base_ClassificationExtendedness_value'] < 0.5
    filtG = d['meas']['base_ClassificationExtendedness_value'] > 0.5

    s = d[cat][filt&filtS].group_by(oid)
    f = (s.groups.indices[1:] - s.groups.indices[:-1]) == nfilters
    star = s.groups[f]
    g = d[cat][filt&filtG].group_by(oid)
    f = (g.groups.indices[1:] - g.groups.indices[:-1]) == nfilters
    gal = g.groups[f]

    # get color magnitudes for stars and galaxies
    mrS = star[star['filter']=='r'][mag_type]
    miS = star[star['filter']=='i'][mag_type]
    mgS = star[star['filter']=='g'][mag_type]
    mrG = gal[gal['filter']=='r'][mag_type]
    miG = gal[gal['filter']=='i'][mag_type]
    mgG = gal[gal['filter']=='g'][mag_type]

    # plot r-i Vs g-r for stars and galaxies

    # cut on magnitude
    idxS = mrS < 23
    idxG = mrG < 22

    fig, (ax1) = P.subplots(ncols=1)
    ax1.scatter(mgS[idxS] - mrS[idxS], mrS[idxS] - miS[idxS], s=1,
                color='b', label="stars %d"%len(mgS[idxS]))
    ax1.scatter(mgG[idxG] - mrG[idxG], mrG[idxG] - miG[idxG], s=1,
                color='r', label="galaxies %d"%len(mgG[idxG]))
    ax1.set_xlim([-0.5, 2.0])
    ax1.set_ylim([-0.5, 2.5])
    ax1.tick_params(labelsize=20)
    ax1.set_xlabel("g - r", fontsize=20)
    ax1.set_ylabel("r - i", fontsize=20)
    ax1.legend(loc="upper left")

    P.show()

    # Coefficients of 5-degree polynomials in (g-i) + residuals for stellar locus parametrization
    # according to Covey et al. 2007 - arXiv:0707.4473
    # The color reference system is SDSS

    uMinusg = N.asarray([1.0636113, -1.6267818, 4.9389572, -3.2809081,
                         0.8725109, -0.0828035, 0.6994, 0.0149])
    gMinusr = N.asarray([0.4290263, -0.8852323, 2.0740616, -1.1091553,
                         0.2397461, -0.0183195, 0.2702, 0.0178])
    rMinusi = N.asarray([-0.4113500, 1.8229991, -1.9989772, 1.0662075,
                         -0.2284455, 0.0172212, 0.2680, 0.0164])
    iMinusz = N.asarray([-0.2270331, 0.7794558, -0.7350749, 0.3727802,
                         -0.0735412, 0.0049808, 0.1261, 0.0071])
    poly_5 = (lambda p, x: p[0] + p[1]*x + p[2]*x**2 + p[3]*x**3 + p[4]*x**4 + p[5]*x**5)

    # Color corrections CFHT --> SDSS
    u_SDSS_ug = (lambda u_Mega, g_Mega : u_Mega +0.181*(u_Mega - g_Mega))
    g_SDSS_gr = (lambda g_Mega, r_Mega : g_Mega +0.195*(g_Mega - r_Mega))
    g_SDSS_gi = (lambda g_Mega, i_Mega : g_Mega +0.103*(g_Mega - i_Mega))
    r_SDSS_gr = (lambda r_Mega, g_Mega : r_Mega +0.011*(g_Mega - r_Mega))
    i_SDSS_ri = (lambda i_Mega, r_Mega : i_Mega +0.079*(r_Mega - i_Mega))
    i_SDSS_gi = (lambda i_Mega, g_Mega : i_Mega +0.044*(g_Mega - i_Mega))
    i2_SDSS_ri = (lambda i2_Mega, r_Mega : i2_Mega +0.001*(r_Mega - i2_Mega))
    i2_SDSS_gi = (lambda i2_Mega, g_Mega : i2_Mega -0.003*(g_Mega - i2_Mega))
    z_SDSS_iz = (lambda z_Mega, i_Mega : z_Mega -0.099*(i_Mega - z_Mega))

    gSDSS = g_SDSS_gr(mgS, mrS)
    rSDSS = r_SDSS_gr(mrS, mgS)
    if ifilt == "i_new":
        print "Will use new i2 (MP9702) filter"
        iSDSS = i2_SDSS_ri(miS, mrS)
        iSDSS2 = i2_SDSS_gi(miS, mgS)
    else:
        print "Will use old i (MP9701) filter"
        iSDSS = i_SDSS_ri(miS, mrS)
        iSDSS2 = i_SDSS_gi(miS, mgS)

    # plot stellar locus for stars and compare to Covey et al. 2007 model

    idxS = (mrS < 23) & (mgS < 23)
    fig, (ax1, ax2) = P.subplots(ncols=2)

    ax1.scatter(gSDSS[idxS]-iSDSS[idxS], rSDSS[idxS]-iSDSS[idxS],
                s=1, color='b', label='%s '%mag_type)
    ax1.set_xlim([-0.5, 4.])
    ax1.set_ylim([-0.3, 2.5])
    nbins = 80
    gMi_model = N.linspace(0.4, 3.8, nbins)
    ax1.plot(gMi_model, poly_5(rMinusi[:6], gMi_model), color='r', label='Covey 2007 model', lw=2)
    ax1.set_xlabel("g - i", fontsize=10)
    ax1.set_ylabel("r - i", fontsize=10)
    ax1.tick_params(labelsize=10)
    ax1.legend(loc="upper left", fontsize=10)

    ax2.scatter(gSDSS[idxS]-iSDSS[idxS], gSDSS[idxS]-rSDSS[idxS],
                s=1, color='b', label='%s' % mag_type)
    ax2.set_xlim([-0.5, 4.])
    ax2.set_ylim([-0.3, 2.0])
    gMr_model = N.linspace(0.4, 3.8, nbins)
    ax2.plot(gMi_model, poly_5(gMinusr[:6], gMr_model), color='r', label='Covey 2007 model', lw=2)
    ax2.set_xlabel("g - i", fontsize=10)
    ax2.set_ylabel("g - r", fontsize=10)
    ax2.tick_params(labelsize=10)
    ax2.legend(loc="upper left", fontsize=10)

    # Plot residuals
    fig, (ax1, ax2) = P.subplots(ncols=2)

    xMin = 0.5
    xMax = 3.0
    numStep = 40
    res = []
    err = []
    colVal, step = N.linspace(xMin,xMax,num=numStep,retstep=True)
    for c in colVal :
        idc = ((gSDSS[idxS]-iSDSS[idxS]) > (c-step/2)) & ((gSDSS[idxS]-iSDSS[idxS]) < (c+step/2))
        res.append(poly_5(rMinusi[:6], c) - N.mean((rSDSS[idxS]-iSDSS[idxS])[idc]))
        err.append(N.std((rSDSS[idxS]-iSDSS[idxS])[idc]) /
                   N.sqrt(len((rSDSS[idxS]-iSDSS[idxS])[idc])))
    ax1.errorbar(colVal, res, yerr=err, ls='None', fmt='s',
                 capsize=25, color='b', lw=1, label='%s - Residuals'%mag_type)
    ax1.set_xlim([xMin-0.2, xMax+0.2])
    ax1.set_ylim([-0.1, 0.1])
    ax1.set_xlabel("g - i", fontsize=10)
    ax1.set_ylabel("Model - data (r-i)", fontsize=10)
    ax1.tick_params(labelsize=10)
    ax1.legend(loc="upper left", fontsize=10)


    xMin = 0.5
    xMax = 3.0
    numStep = 40
    res = []
    err = []
    colVal, step = N.linspace(xMin,xMax,num=numStep,retstep=True)
    for c in colVal:
        idc = ((gSDSS[idxS]-iSDSS[idxS]) > (c-step/2)) & ((gSDSS[idxS]-iSDSS[idxS]) < (c+step/2))
        res.append(poly_5(gMinusr[:6], c) - N.mean((gSDSS[idxS]-rSDSS[idxS])[idc]))
        err.append(N.std((gSDSS[idxS]-rSDSS[idxS])[idc]) /
                   N.sqrt(len((gSDSS[idxS]-rSDSS[idxS])[idc])))
    ax2.errorbar(colVal, res, yerr=err, ls='None', fmt='s', capsize=25,
                 color='b', lw=1, label='%s - Residuals'%mag_type)
    ax2.set_xlim([xMin-0.2, xMax+0.2])
    ax2.set_ylim([-0.1, 0.1])
    ax2.set_xlabel("g - i", fontsize=10)
    ax2.set_ylabel("Model - data (g-r)", fontsize=10)
    ax2.tick_params(labelsize=10)
    ax2.legend(loc="upper left", fontsize=10)
    P.tight_layout()
    P.show()

def starElipticities(d) :

    """Compute star elipticities from second momments and check if
    psf correction is valid - Also check magnitude Vs radius"""

    nfilters = len(d['meas'].group_by('id').groups[0])
    cat = 'meas'
    oid = 'id'

    filt = d['meas']['base_ClassificationExtendedness_flag'] == 0
    filt &= d['meas']['detect_isPrimary'] == 1
    filt &= d[cat]['modelfit_CModel_flag'] == 0
    filt &= d[cat]['modelfit_CModel_flux'] > 0

    filt &= (d[cat]['modelfit_CModel_flux'] / \
              d[cat]['modelfit_CModel_fluxSigma']) > 10

    filtS = d['meas']['base_ClassificationExtendedness_value'] < 0.5
    filtG = d['meas']['base_ClassificationExtendedness_value'] > 0.5

    s = d[cat][filt&filtS].group_by(oid)
    f = (s.groups.indices[1:] - s.groups.indices[:-1]) == nfilters
    star = s.groups[f]
    g = d[cat][filt&filtG].group_by(oid)
    f = (g.groups.indices[1:] - g.groups.indices[:-1]) == nfilters
    gal = g.groups[f]

    shapeHSMSource_xx_s = star[star['filter'] == 'i']['ext_shapeHSM_HsmSourceMoments_xx']
    shapeHSMSource_yy_s = star[star['filter'] == 'i']['ext_shapeHSM_HsmSourceMoments_yy']
    shapeHSMSource_xy_s = star[star['filter'] == 'i']['ext_shapeHSM_HsmSourceMoments_xy']
    shapeHSMPsf_xx_s = star[star['filter'] == 'i']['ext_shapeHSM_HsmPsfMoments_xx']
    shapeHSMPsf_yy_s = star[star['filter'] == 'i']['ext_shapeHSM_HsmPsfMoments_yy']
    shapeHSMPsf_xy_s = star[star['filter'] == 'i']['ext_shapeHSM_HsmPsfMoments_xy']
    magI_s = star[star['filter']=='i']['modelfit_CModel_mag']
    radius_s = r = N.sqrt(shapeHSMSource_xx_s + shapeHSMSource_yy_s)

    shapeHSMSource_xx_g = gal[gal['filter'] == 'i']['ext_shapeHSM_HsmSourceMoments_xx']
    shapeHSMSource_yy_g = gal[gal['filter'] == 'i']['ext_shapeHSM_HsmSourceMoments_yy']
    shapeHSMSource_xy_g = gal[gal['filter'] == 'i']['ext_shapeHSM_HsmSourceMoments_xy']
    magI_g = gal[gal['filter'] == 'i']['modelfit_CModel_mag']
    radius_g = r = N.sqrt(shapeHSMSource_xx_g + shapeHSMSource_yy_g)

    # Plot magnitude as a function of the source radius computed from second momments
    fig, (ax1, ax2) = P.subplots(ncols=2)
    ax1.scatter(radius_s, magI_s, s=1, color='b', label='Stars %d'%len(magI_s))
    ax1.scatter(radius_g, magI_g, s=1, color='r', label='Galaxies %d'%len(magI_g))
    ax1.set_xlim([1., 6.])
    ax1.set_ylim([16, 26])
    ax1.set_xlabel('Radius in pixels', fontsize=10)
    ax1.set_ylabel('Magnitude i', fontsize=10)
    ax1.tick_params(labelsize=10)
    ax1.legend(loc="lower left", fontsize=10)

    ax2.hist(radius_s[magI_s < 23], bins=80, range=[1.7, 3], color='b')
    ax2.hist(radius_g[magI_g < 23], bins=80, range=[1.7, 3], color='r')
    ax1.set_xlabel('Radius in pixels', fontsize=10)

    P.tight_layout()
    P.show()

    denomSource = shapeHSMSource_xx_s +2.*N.sqrt(shapeHSMSource_xx_s*shapeHSMSource_yy_s - N.square(shapeHSMSource_xy_s))
    e1source = (shapeHSMSource_xx_s - shapeHSMSource_yy_s) / denomSource
    e2source = 2.0*shapeHSMSource_xy_s / denomSource

    denomPsf = shapeHSMPsf_xx_s +2.*N.sqrt(shapeHSMPsf_xx_s*shapeHSMPsf_yy_s - N.square(shapeHSMPsf_xy_s))
    e1psf = (shapeHSMPsf_xx_s - shapeHSMPsf_yy_s) / denomSource
    e2Psf = 2.0*shapeHSMPsf_xy_s / denomSource

    idx = magI_s < 22.5
    fig, (ax0, ax1) = P.subplots(ncols=2)
    ax0.scatter(e1source[idx], e2source[idx], s=1, color='b')
    ax0.set_xlabel('e1(source)', fontsize=10)
    ax0.set_ylabel('e2(source)', fontsize=10)
    ax0.set_xlim([-0.10, 0.10])
    ax0.set_ylim([-0.10, 0.10])
    ax1.tick_params(labelsize=10)
    ax1.scatter(e1source[idx]-e1psf[idx], e2source[idx]-e2Psf[idx], s=1, color='b')
    ax1.set_xlim([-0.10, 0.10])
    ax1.set_ylim([-0.10, 0.10])
    ax1.set_xlabel('e1(source) - e1(psf)', fontsize=10)
    ax1.set_ylabel('e2(source) - e2(psf)', fontsize=10)
    ax1.tick_params(labelsize=10)

    P.tight_layout()
    P.show()
