import numpy as np
import pzmassfitter as pmf
import pzmassfitter.nfwutils as nfwutils
import pzmassfitter.nfwmodeltools as nfwmodeltools
import astropy.table as table
import tempfile, shutil
import yaml
import clusters.main
import cPickle
import astropy.io.fits as pyfits
import pzmassfitter.ldac as ldac

def setupNearPerfectData(m200 = 1e15):

    ngals = 10000

    seqnr = np.arange(ngals)
    
    zcluster = 0.3
    z_source0 = 0.8
    sigma_z = 0.05
    z_true = z_source0 + sigma_z*np.random.standard_normal(ngals)
    zrange = np.arange(0.5, 1.1, 0.05)
    z_pdf = np.exp(-0.5*((zrange - z_source0)/sigma_z)**2)/np.sqrt(2*np.pi*sigma_z**2)
    assert np.abs(np.trapz(z_pdf, zrange) - 1) < 1e-6

    z_best = z_source0*np.ones(ngals)
    pdfgrid = np.vstack(ngals*[z_pdf])


    beta_true = nfwutils.global_cosmology.beta_s(z_true, zcluster)
    print 'beta', beta_true
    
    x_mpc = np.random.uniform(-4, 4, size=ngals)
    y_mpc = np.random.uniform(-4, 4, size=ngals)
    r_mpc = np.sqrt(x_mpc**2 + y_mpc**2)

    Dl = nfwutils.global_cosmology.angulardist(zcluster)

    x_deg = (x_mpc/Dl)*(180./np.pi)
    y_deg = (y_mpc/Dl)*(180./np.pi)

    rho_c = nfwutils.global_cosmology.rho_crit(zcluster)


    c200 = 4.
    r200 = (3*abs(m200)/(4*200*np.pi*rho_c))**(1./3.)
    rscale = r200 / c200

    
    rho_c_over_sigma_c = 1.5 * nfwutils.global_cosmology.angulardist(zcluster) * nfwutils.global_cosmology.beta([1e6], zcluster)[0] * nfwutils.global_cosmology.hubble2(zcluster) / nfwutils.global_cosmology.v_c**2
    gamma_inf = nfwmodeltools.NFWShear(r_mpc, c200, rscale, rho_c_over_sigma_c, 200.)
    kappa_inf = nfwmodeltools.NFWKappa(r_mpc, c200, rscale, rho_c_over_sigma_c, 200.)

    g_t_true = beta_true*gamma_inf/(1-beta_true*kappa_inf)
    g_x_true = np.zeros_like(g_t_true)

    print 'shear', g_t_true

    shapenoise = 0.005
    g_t = g_t_true + shapenoise*np.random.standard_normal(ngals)

    posangle = np.arctan2(y_mpc, x_mpc)

    cos2phi = np.cos(2*posangle)
    sin2phi = np.sin(2*posangle)

    e1 = -g_t*cos2phi
    e2 = -g_t*sin2phi

    deepCoadd_meas = table.Table([seqnr, x_deg, y_deg, e1, e2], names=('id', 'coord_ra_deg', 'coord_dec_deg', 'ext_shapeHSM_HsmShapeRegauss_e1', 'ext_shapeHSM_HsmShapeRegauss_e2'))

    pdz_values = table.Table([seqnr, z_best, pdfgrid], names=('objectId', 'Z_BEST', 'pdz'))
    pdz_bins = table.Table([zrange,], names=('zbins',))


    tmpdir = tempfile.mkdtemp()

    shearcatfile = '{}/shearfile.hdf5'.format(tmpdir)
    deepCoadd_meas.write(shearcatfile, path='deepCoadd_meas')

    pdzfile = '{}/pdz.hdf5'.format(tmpdir)
    pdz_values.write(pdzfile, path='pdz_values')
    pdz_bins.write(pdzfile, path='pdz_bins', append=True)

    config=dict(cluster = 'nearperfect',
                ra = 0.,
                dec = 0.,
                redshift = zcluster,
                filter = ["u", "g", "r", "i", "i2", "z"])

    
    configfile = '{}/config.yaml'.format(tmpdir)
    with open(configfile, 'w') as output:
        output.write(yaml.dump(config))

    return tmpdir, shearcatfile, pdzfile, configfile



######


def cleanuptest(tmpdir):

    shutil.rmtree(tmpdir)

######


def test_pzmassfitter():

    try:

        m200 = 1e15

        tmpdir, shearcatfile, pdzfile, configfile = setupNearPerfectData(m200)

        argv = '--zcode none --testing --output {tmpdir}/mass.out {configfile} {shearcatfile} {pdzfile}'.format(tmpdir=tmpdir,
                                                                                                                configfile=configfile,
                                                                                                                shearcatfile=shearcatfile,
                                                                                                                pdzfile=pdzfile).split()

        clusters.main.mass(argv)


        scanfile = '{tmpdir}/mass.out.m200.scan.fits'.format(tmpdir=tmpdir)
        scan = ldac.openObjectFile(scanfile)
        maxlike_mass = scan['Mass'][scan['prob'] == np.max(scan['prob'])]
        print m200, maxlike_mass, maxlike_mass/m200
        

        assert np.abs(maxlike_mass - m200)/m200 < 0.02

    except AssertionError:
        print tmpdir
        raise

    
    
    cleanuptest(tmpdir)




if __name__ == '__main__':

    test_pzmassfitter()
    
