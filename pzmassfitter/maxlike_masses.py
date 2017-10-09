"""Run a ML pdz fit for an nfw model."""

from __future__ import with_statement
try:
    import cPickle as pickle  # python 2
except ImportError:
    import pickle  # python 3
import pymc
import numpy as np
import astropy.io.fits as pyfits
from . import ldac
from . import nfwutils
from . import util
from . import nfwmodeltools as tools
from . import pymc_mymcmc_adapter as pma


usage = '''
   maxlike_masses.py  filehandler_module shapedistro_module <options>

maxlike_masses.py may be called from the command line or from
python. On the command line, the first two arguments must be the names
of python modules that 1.) handle reading in file IO and 2.) contain
information on the shape distribution. Module names should not have
the .py on them.

ex: maxlike_masses.py maxlike_sim_filehandler nfwmodel_normshapedistro -o outfile -n 50000 ...

Use the -h option, after specifying the filehandler and shapedistro to see all available options.
'''


######################
# ML Reconstruction
######################


class ModelInitException(Exception):
    pass


massscale = 1e14


class LensingModel(object):

    def __init__(self):

        self.cuts = [self.modelCut]

    def addCLOps(self, parser):

        parser.add_option('--deltaz95low', dest='deltaz95low',
                          help='Lower limit on the width of the PDZ',
                          type='float', default=-1)
        parser.add_option('--deltaz95high', dest='deltaz95high',
                          help='Upper limit on the with of the PDZ',
                          type='float', default=2.5)
        parser.add_option('--zbhigh', dest='zbhigh',
                          help='Upper limit on photoz point estimate',
                          type='float', default=1.25)
        parser.add_option('--zcut', dest='zcut',
                          help='To what deltaZ behind the cluster should galaxies be exlucded?',
                          type='float', default=0.1)
        parser.add_option('--masslow', dest='masslow',
                          help='Mass prior low cutoff',
                          type='float', default=1e13)
        parser.add_option('--masshigh', dest='masshigh',
                          help='Mass prior high cutoff',
                          type='float', default=1e16)
        parser.add_option('--ztypecut', dest='ztypecut',
                          help='Turn on type dependent redshift cuts',
                          default=False, action='store_true')
        parser.add_option('--radlow', dest='radlow',
                          help='Low radius cutoff wrt cluster center, Mpc',
                          default=0.75, type='float')
        parser.add_option('--radhigh', dest='radhigh',
                          help='High radius cutoff wrt cluster center, Mpc',
                          default=3.0, type='float')
        parser.add_option('--concentration', dest='concentration',
                          help='Assumed concentration in the fit',
                          default=None, type='float')
        parser.add_option('--logprior', dest='logprior',
                          help='Turn on log10 mass prior',
                          default=False, action='store_true')
        
    #######################################################

    def createOptions(self, deltaz95low=-1, deltaz95high=2.5, zbhigh=5,  # zbhigh=1.25 default
                      zcut=0.1, masslow=1e13, masshigh=1e16,
                      ztypecut=False, radlow=0.75, radhigh=3.0,  # radlow=0.75 default
                      concentration=None, delta=200.,
                      options=None, args=None, logprior=False):

        if options is None:
            options = util.VarContainer()

        options.deltaz95low = deltaz95low
        options.deltaz95high = deltaz95high
        options.zbhigh = zbhigh
        options.zcut = zcut
        options.masslow = masslow
        options.masshigh = masshigh
        options.ztypecut = ztypecut
        options.radlow = radlow
        options.radhigh = radhigh
        options.concentration = concentration
        options.delta = delta
        options.logprior = logprior

        return options, None

    #######################################################

    def modelCut(self, manager, minMPC=0.5, maxMPC=3.):

        options = manager.options
        inputcat = manager.inputcat

#        if 'r500' in manager:
 #           manager.comment('Using r500')
        minMPC = options.radlow  # *manager.r500
        maxMPC = options.radhigh  # *manager.r500

        goodObjs = np.logical_and(
            np.logical_and(np.logical_and(manager.inputcat['r_mpc'] > minMPC,
                                          manager.inputcat['r_mpc'] < maxMPC),
                           np.logical_and(manager.inputcat['z_b'] > 0,
                                          manager.inputcat['z_b'] < options.zbhigh)),
            np.abs(manager.inputcat['ghats']) < 5)

# Definition of pdz is changing, need to change this.
#        pdz = manager.pz # used to be manager.pdz (name change in astropytable_filehandler)
#        pdzrange = manager.pdzrange
#        delta95Z = np.zeros(len(pdz))
#
#        for i in range(len(delta95Z)):
#
#            cumpdz = pdz[i].cumsum() / pdz[i].cumsum()[-1]
#
#            delta95Z[i] = pdzrange[cumpdz >= 0.95][0] - pdzrange[cumpdz >= 0.05][0]
#
#        deltaZcut = np.logical_and(options.deltaz95low <= delta95Z,
#                                   delta95Z < options.deltaz95high)
##
####
        if options.zcut is None:

            zcut = np.ones(len(manager.inputcat)) == 1

        else:

            zcut = manager.inputcat['z_b'] > (manager.zcluster + options.zcut)

        ztypecut = np.ones(len(inputcat)) == 1
        if options.ztypecut:
            zt = inputcat['z_t']
            zb = inputcat['z_b']

            type1 = np.logical_and(zt >= 1, zt < 2)
            ztypecut[type1] = np.logical_and(type1, zb < 1.15)

            type2 = np.logical_and(zt >= 2, zt < 3)
            ztypecut[type2] = np.logical_and(type2, zb < 1.3)

            type3 = np.logical_and(zt >= 3, zt < 4)
            ztypecut[type3] = np.logical_and(type3,
                                             np.logical_or(zb <= 1,
                                                           np.logical_and(1.15 < zb, zb < 1.3)))

            type4 = np.logical_and(zt >= 4, zt < 5)
            ztypecut[type4] = np.logical_and(type4, np.logical_or(
                zb < 0.95, np.logical_and(1.15 < zb, zb < 1.3)))

            type5 = np.logical_and(zt >= 5, zt < 6)
            ztypecut[type5] = False

        # basic_cuts = reduce(np.logical_and, [goodObjs, deltaZcut, zcut, ztypecut])
        basic_cuts = reduce(np.logical_and, [goodObjs, zcut, ztypecut])

        return basic_cuts

    ########################################################################################

    def makeModelPrior(self, manager, parts):

        options = manager.options

        if options.concentration is None:
            parts.log10concentration = pymc.TruncatedNormal('log10concentration', 0.6, 1. / 0.116**2,
                                                            np.log10(1.), np.log10(10.))  # tau!

            @pymc.deterministic
            def cdelta(log10concentration=parts.log10concentration):
                return 10**log10concentration
            parts.cdelta = cdelta
        else:
            parts.cdelta = options.concentration

        manager.massdelta = options.delta
        parts.massdelta = options.delta

        if options.logprior:
            # Uniform sampling of log(m)
            parts.log10mdelta = pymc.Uniform('log10mdelta', np.log10(
                options.masslow), np.log10(options.masshigh))

            @pymc.deterministic
            def mdelta(log10mdelta=parts.log10mdelta):
                return 10**log10mdelta
            parts.mdelta = mdelta

        else:
         # Uniform sampling of m
            parts.scaledmdelta = pymc.Uniform(
                'scaledmdelta', options.masslow / massscale, options.masshigh / massscale)

            @pymc.deterministic
            def mdelta(scaledmdelta=parts.scaledmdelta):
                return massscale * scaledmdelta
            parts.mdelta = mdelta

        #############################

    def makeShapePrior(self, datamanager, parts):
        # This is just a stand-in. Subclass for specific examples.

        inputcat = datamanager.inputcat

        parts.shearcal_m = np.zeros(len(inputcat))
        parts.shearcal_c = np.zeros(len(inputcat))

        parts.sigma = 0.005

    ##############################################################
    # Likelihood
    ###########

    def makeLikelihood(self, datamanager, parts):

        inputcat = datamanager.inputcat

        pz = datamanager.pz

        parts.r_mpc = np.ascontiguousarray(
            inputcat['r_mpc'].astype(np.float64))
        parts.ghats = np.ascontiguousarray(
            inputcat['ghats'].astype(np.float64))
        parts.pz = np.ascontiguousarray(pz.astype(np.float64))

        parts.zs = np.ascontiguousarray(
            np.array(datamanager.pdzrange).astype(np.float64))

        parts.betas = np.ascontiguousarray(nfwutils.global_cosmology.beta_s(
            parts.zs, parts.zcluster).astype(np.float64))
        parts.nzbins = len(parts.betas)

        parts.rho_c = nfwutils.global_cosmology.rho_crit(parts.zcluster)
        parts.rho_c_over_sigma_c = 1.5 * nfwutils.global_cosmology.angulardist(parts.zcluster) * \
            nfwutils.global_cosmology.beta([1e6], parts.zcluster)[0] * \
            nfwutils.global_cosmology.hubble2(parts.zcluster) / \
            nfwutils.global_cosmology.v_c**2

        parts.data = None
        for i in range(20):
            try:
                @pymc.stochastic(observed=True, name='data_%d' % i)
                def data(value=parts.ghats,
                         mdelta=parts.mdelta,
                         cdelta=parts.cdelta,
                         r_mpc=parts.r_mpc,
                         zs=parts.zs,
                         betas=parts.betas,
                         pz=parts.pz,
                         shearcal_m=parts.shearcal_m,
                         shearcal_c=parts.shearcal_c,
                         sigma=parts.sigma,
                         rho_c=parts.rho_c,
                         rho_c_over_sigma_c=parts.rho_c_over_sigma_c,
                         massdelta=parts.massdelta):

                    return tools.gauss_like(mdelta,
                                            cdelta,
                                            r_mpc,
                                            value,
                                            zs,
                                            betas,
                                            pz,
                                            shearcal_m,
                                            shearcal_c,
                                            sigma,
                                            rho_c,
                                            rho_c_over_sigma_c,
                                            massdelta)

                parts.data = data
                break
            except pymc.ZeroProbability:
                pass

        if parts.data is None:
            raise ModelInitException

        #######################

    def makeModelParts(self, datamanager, parts=None):

        if parts is None:
            parts = util.VarContainer()

        parts.zcluster = datamanager.zcluster

        for i in range(10):
            try:
                self.makeShapePrior(datamanager, parts)
                self.makeModelPrior(datamanager, parts)
                self.makeLikelihood(datamanager, parts)
                return parts
            except pymc.ZeroProbability:
                pass

        raise ModelInitException

    #############

    def createModel(self, datamanager):
        datamanager.ngalaxies = len(datamanager.inputcat)
        parts = self.makeModelParts(datamanager)
        return pymc.Model(parts)


#########################################################################
#########################################################################


class ScanModelToFile(object):

    def addCLOps(self, parser):
        pass

    ################

    def createOptions(self,
                      outputFile,
                      options=None, args=None):

        if options is None:
            options = util.VarContainer()

        options.outputFile = outputFile
        return options, args

    ##########

    def run(self, manager):

        # SCANNING MASS(<1.5MPC)

        mass = np.arange(5e13, 1e16, 5e12)
        model = manager.model

        scan = np.zeros_like(mass)
        for i, m in enumerate(mass):
            try:
                model.scaledmdelta.value = m / massscale
                scan[i] = model.logp
            except pymc.ZeroProbability:
                scan[i] = pymc.PyMCObjects.d_neg_inf

        cols = [pyfits.Column(name='Mass', format='E', array=mass),
                pyfits.Column(name='prob', format='E', array=scan)]
        manager.cat = ldac.LDACCat(
            pyfits.BinTableHDU.from_columns(pyfits.ColDefs(cols)))
        manager.cat.hdu.header.set('EXTNAME', 'OBJECTS')

        manager.cat.saveas('{}.m{}.scan.fits'.format(manager.options.outputFile,
                                                     int(manager.model.massdelta)))

#        self.calcMasses(manager)

    ##########

    def calcMasses(self, manager):

        masses = manager.cat['Mass']
        scan = manager.cat['prob']

        pdf = np.exp(scan - max(scan))
        cdf = np.cumsum(pdf)
        cdf = cdf / cdf[-1]

        buffered_cdf = np.zeros(len(cdf) + 2)
        buffered_cdf[-1] = 1.
        buffered_cdf[1:-1] = cdf

        nsamples = 10000
        manager.masses = np.zeros(nsamples)
        for i in range(nsamples):

            cdf_pick = np.random.uniform()
            inbin = np.logical_and(
                np.roll(buffered_cdf, 1) <= cdf_pick, buffered_cdf > cdf_pick)
            manager.masses[i] = masses[inbin[1:-1]][0]

    ##########

    def dump(self, manager):
        pass

#        manager.cat.saveas(manager.options.outputFile, clobber=True)
#
#
#
#        outputFile = manager.options.outputFile
#        pma.dumpMasses(manager.masses,'%s.mass15mpc' % outputFile)
#

    ##########

    def finalize(self, manager):
        pass

################################################


class SampleModelToFile(object):

    def run(self, manager):

        model = manager.model
        nsamples = manager.options.nsamples
        outputFile = manager.options.outputFile
        burn = manager.options.burn

        mcmc_manager = util.VarContainer()
        mcmc_options = util.VarContainer()

        mcmc_manager.options = mcmc_options

        mcmc_options.singlecore = True
        mcmc_options.adapt_every = 100
        mcmc_options.adapt_after = 100
        mcmc_options.nsamples = nsamples
        mcmc_manager.model = model

        runner = pma.MyMCMemRunner()
        runner.run(mcmc_manager)
        runner.finalize(mcmc_manager)

        manager.chain = mcmc_manager.chain

    def addCLOps(self, parser):

        raise NotImplementedError

    def createOptions(self, outputFile, nsamples=2000, burn=500, options=None, args=None):

        if options is None:
            options = util.VarContainer()

        options.outputFile = outputFile
        options.nsamples = nsamples
        options.burn = burn
        return options, args

    def dump(self, manager):

        outputFile = manager.options.outputFile

        with open('%s.chain.pkl' % outputFile, 'wb') as output:
            pickle.dump(manager.chain, output)

        pma.dumpMasses(np.array(manager.chain['mdelta'][manager.options.burn:]),
                       '%s.m%d' % (outputFile, manager.massdelta))

    def finalize(self, manager):
        pass
