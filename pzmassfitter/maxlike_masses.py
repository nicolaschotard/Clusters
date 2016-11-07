#!/usr/bin/env python
#######################
# Run a ML pdz fit for an nfw model
########################

from __future__ import with_statement
import cPickle
import numpy as np, pymc
import astropy.io.fits as pyfits
import nfwmodeltools as tools, varcontainer
import nfwutils, shearprofile as sp, ldac

##########################

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


class ModelInitException(Exception): pass    


##########################

class LensingModel(object):

    def __init__(self):

        self.cuts = [self.modelCut]

    

    def makeLikelihood(self, datamanager, parts):

        inputcat = datamanager.inputcat

        pdz = datamanager.pdz


        bin_selectors = self.bin_selectors(inputcat)


        parts.r_mpc = [ np.ascontiguousarray(inputcat['r_mpc'][x].astype(np.float64)) for x in bin_selectors ]
        parts.ghats = [ np.ascontiguousarray(inputcat['ghats'][x].astype(np.float64)) for x in bin_selectors ]
        parts.pdz = [ np.ascontiguousarray(pdz[x].astype(np.float64)) for x in bin_selectors ]


        parts.nshapebins = len(bin_selectors)
        parts.pdzrange = datamanager.pdzrange

        
        parts.betas = np.ascontiguousarray(nfwutils.beta_s(parts.pdzrange, parts.zcluster).astype(np.float64))
        parts.nzbins = len(parts.betas)


        parts.data = np.empty(parts.nshapebins, dtype=object)


        for i, cur_ghats, cur_shape_param, cur_r_mpc, cur_pdz in zip(np.arange(parts.nshapebins), 
                                                                     parts.ghats, 
                                                                     parts.shape_params, 
                                                                     parts.r_mpc, 
                                                                     parts.pdz):

            @pymc.stochastic(observed=True, name='data_%d' % i)
            def data(value = cur_ghats, r_scale = parts.r_scale, shape_params = cur_shape_param,
                     r_mpc = cur_r_mpc, pdz = cur_pdz, betas = parts.betas, 
                     concentration = parts.concentration,
                     zcluster = parts.zcluster):


                return self.likelihood_func(r_scale, r_mpc, 
                                       value, betas, 
                                       pdz, shape_params,
                                       concentration, zcluster)





            parts.data[i] = data





    #######################################################

    def addCLOps(self, parser):

        parser.add_option('--deltaz95low', dest='deltaz95low',
                          help = 'Lower limit on the width of the PDZ',
                          type='float', default = -1)
        parser.add_option('--deltaz95high', dest='deltaz95high',
                          help = 'Upper limit on the with of the PDZ',
                          type='float', default = 2.5)
        parser.add_option('--zbhigh', dest='zbhigh',
                          help = 'Upper limit on photoz point estimate',
                          type='float', default = 1.25)
        parser.add_option('--zcut', dest='zcut',
                          help = 'To what deltaZ behind the cluster should galaxies be exlucded?',
                          type='float', default = 0.1)
        parser.add_option('--masslow', dest='masslow',
                          help = 'Mass prior low cutoff',
                          type = 'float', default = 1e13)
        parser.add_option('--masshigh', dest='masshigh',
                          help = 'Mass prior high cutoff',
                          type = 'float', default = 1e16)
        parser.add_option('--ztypecut', dest='ztypecut',
                          help='Turn on type dependent redshift cuts',
                          default = False, action='store_true')
        parser.add_option('--radlow', dest='radlow',
                          help = 'Low radius cutoff wrt cluster center, Mpc',
                          default = 0.75, type = 'float')
        parser.add_option('--radhigh', dest='radhigh',
                          help = 'High radius cutoff wrt cluster center, Mpc',
                          default = 3.0, type='float')
        parser.add_option('--concentration', dest='concentration',
                          help = 'Assumed concentration in the fit',
                          default = None, type='float')
        parser.add_option('--logprior', dest='logprior',
                          help = 'Turn on log10 mass prior',
                          default = False, action = 'store_true')

                          


    #######################################################



    def createOptions(self, deltaz95low = -1, deltaz95high = 2.5, zbhigh = 1.25,
                      zcut = 0.1, masslow = 1e13, masshigh = 1e16,
                      ztypecut = False, radlow = 0.75, radhigh = 3.0, 
                      concentration = None, logprior = False,
                      options = None, args = None):

        if options is None:
            options = varcontainer.VarContainer()

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
        options.logprior = logprior

        return options, None




    #######################################################


    def modelCut(self, manager, minMPC=0.5, maxMPC=3.):

        options = manager.options
        inputcat = manager.inputcat

#        if 'r500' in manager:
 #           manager.comment('Using r500')
        minMPC = options.radlow   #*manager.r500
        maxMPC = options.radhigh  #*manager.r500

        goodObjs = np.logical_and( \
                   np.logical_and(np.logical_and(manager.inputcat['r_mpc'] > minMPC, 
                                                 manager.inputcat['r_mpc'] < maxMPC),
                                  np.logical_and(manager.inputcat['z_b'] > 0,
                                                 manager.inputcat['z_b'] < options.zbhigh)),
                   np.abs(manager.inputcat['ghats']) < 5)

        



        pdz = manager.pdz
        pdzrange = manager.pdzrange
        delta95Z = np.zeros(len(pdz))

        for i in range(len(delta95Z)):

            cumpdz = pdz[i].cumsum() / pdz[i].cumsum()[-1]

            delta95Z[i] = pdzrange[cumpdz >= 0.95][0] - pdzrange[cumpdz >= 0.05][0]

        deltaZcut = np.logical_and(options.deltaz95low <= delta95Z,
                                   delta95Z < options.deltaz95high)

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
            ztypecut[type3] = np.logical_and(type3, np.logical_or(zb <= 1, np.logical_and( 1.15 < zb, zb < 1.3)))
            
            type4 = np.logical_and(zt >= 4, zt < 5)
            ztypecut[type4] = np.logical_and(type4, np.logical_or(zb < 0.95, np.logical_and(1.15 < zb, zb < 1.3)))

            type5 = np.logical_and(zt >= 5, zt < 6)
            ztypecut[type5] = False

        
            
        basic_cuts = reduce(np.logical_and, [goodObjs, deltaZcut, zcut, ztypecut])


        return basic_cuts





    ########################################################################################


    def makeModelPrior(self, manager, parts):

        options = manager.options


        if options.concentration is None:
            parts.log10concentration = pymc.TruncatedNormal('log10concentration', 0.6, 1./0.116**2, 
                                                      np.log10(1.), np.log10(10.))   #tau!
            @pymc.deterministic
            def concentration(log10concentration = parts.log10concentration):
                return 10**log10concentration
            parts.concentration = concentration
        else:
            parts.concentration = options.concentration


        if options.logprior:
            
            parts.logmass_15mpc = pymc.Uniform('logmass_15mpc', np.log10(options.masslow), 
                                               np.log10(options.masshigh))

            @pymc.deterministic
            def mass_15mpc(logmass = parts.logmass_15mpc):
                return 10**logmass

            parts.mass_15mpc = mass_15mpc

        else:


            parts.mass_15mpc = pymc.Uniform('mass_15mpc', options.masslow, options.masshigh)

        @pymc.deterministic
        def r_scale(mass = parts.mass_15mpc, 
                    concentration = parts.concentration, 
                    zcluster = parts.zcluster):
            
            try:
                rs = nfwutils.RsMassInsideR(mass, concentration, zcluster, 1.5)
            except ValueError:
                raise pymc.ZeroProbability

            return rs

        parts.r_scale = r_scale



    ########################################################################################
    # Shape PRIORS
    ########################

    def makeFixedPrior(self, manager, parts):

        parts.shape_params = self.shapedistro_params

    ######################


    def makeSampledPrior(self, manager, parts):


        parts.shape_sample_index = [pymc.DiscreteUniform('shape_index_%d' % i, 0, len(x)-1) \
                                        for i, x in enumerate(self.shapedistro_params)]

        parts.shape_params = np.empty(self.nshapebins, dtype=object)
        for i, index, samples in zip(range(self.nshapebins), 
                                     parts.shape_sample_index, 
                                     self.shapedistro_params):

            @pymc.deterministic(name = 'shape_params_%d' % i)
            def shape_param_func(index = index, samples = samples):
                return np.ascontiguousarray(samples[:,index])

            parts.shape_params[i] = shape_param_func


    #######################

    def makeNormPrior(self, model, parts):

        parts.shape_params = np.empty(self.nshapebins, dtype=object)
        for i, (mu, cov) in enumerate(self.shapedistro_params):
            parts.shape_params[i] = pymc.MvNormalCov('shape_params_%d' % i, mu, cov)

    ##############################################################



    def makeModelParts(self, datamanager, parts = None):

        if parts is None:
            parts = varcontainer.VarContainer()

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

    ##########

    def run(self, manager):

        #SCANNING MASS(<1.5MPC)

        mass = np.arange(5e13, 1e16, 5e12)
        model = manager.model

        scan = np.zeros_like(mass)
        for i, m in enumerate(mass):
            try:
                model.mass_15mpc.value = m
                scan[i] = model.logp
            except pymc.ZeroProbability:
                scan[i] =  pymc.PyMCObjects.d_neg_inf

        


        cols = [ pyfits.Column(name = 'Mass', format = 'E', array = mass),
                 pyfits.Column(name = 'prob', format = 'E', array = scan)]
        manager.cat = ldac.LDACCat(pyfits.BinTableHDU.from_columns(pyfits.ColDefs(cols)))
        manager.cat.hdu.header.update('EXTNAME', 'OBJECTS')


        self.calcMasses(manager)


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
            inbin = np.logical_and(np.roll(buffered_cdf, 1) <= cdf_pick, buffered_cdf > cdf_pick)
            manager.masses[i] = masses[inbin[1:-1]][0]
            



        



    ##########

    def dump(self, manager):

        manager.cat.saveas(manager.options.outputFile, clobber=True)


        
        outputFile = manager.options.outputFile
        dumpMasses(manager.masses,'%s.mass15mpc' % outputFile)





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


        manager.mcmc = pymc.MCMC(input = model, db='pickle', dbname=outputFile)


        try:
            manager.shapedistro.sampler_callback(mcmc)
        except:
            #doesn't have it, doesn't matter
            pass


        manager.mcmc.sample(nsamples)

        if isinstance(manager.mcmc.concentration, float):
            manager.masses = np.array([ nfwutils.massInsideR(rs,
                                                             manager.mcmc.concentration,
                                                             manager.mcmc.zcluster,
                                                             manager.r500) \
                                            for rs in manager.mcmc.trace('r_scale')[burn:] ])
        else:
            manager.masses = np.array([ nfwutils.massInsideR(rs,
                                                             c,
                                                             manager.mcmc.zcluster,
                                                             manager.r500) \
                                            for rs, c in zip(manager.mcmc.trace('r_scale')[burn:],
                                                             manager.mcmc.trace('concentration')[burn:])])
            
        



    ############

    def addCLOps(self, parser):


        parser.add_option('-s', '--nsamples', dest='nsamples',
                          help='Number of MCMC samples to draw or scan model', default=None, type='int')
        parser.add_option('--burn', dest='burn',
                          help='Number of MCMC samples to discard before calculated mass statistics', 
                          default=10000, type=int)


    ##############

    def dump(self, manager):

        masses = manager.masses

        outputFile = manager.options.outputFile

        dumpMasses(masses, '%s.mx500' % outputFile)

        dumpMasses(manager.mcmc.trace('mass_15mpc')[manager.options.burn:],
                   '%s.mass15mpc' % outputFile)


    ##############

    def finalize(self, manager):


        manager.mcmc.db.close()


    ######################



def dumpMasses(masses, outputFile):


    with open('%s.mass.pkl' % outputFile, 'wb') as output:
        cPickle.dump(masses, output)
        
        
    mean = np.mean(masses)
    stddev = np.std(masses)
    quantiles = pymc.utils.quantiles(masses, qlist=[2.5, 15.8, 25, 50, 75, 84.1, 97.5])
    hpd68 = pymc.utils.hpd(masses, 0.32)
    hpd95 = pymc.utils.hpd(masses, 0.05)
    ml, (m, p) = sp.ConfidenceRegion(masses)
    lml, (lm, lp) = sp.ConfidenceRegion(np.log10(masses))


    with open('%s.mass.summary.txt' % outputFile, 'w') as output:
        output.write('mean\t%e\n' % mean)
        output.write('stddev\t%e\n' % stddev)
        output.write('Q2.5\t%e\n' % quantiles[2.5])
        output.write('Q25\t%e\n' % quantiles[25])
        output.write('Q50\t%e\n' % quantiles[50])
        output.write('Q75\t%e\n' % quantiles[75])
        output.write('Q97.5\t%e\n' % quantiles[97.5])
        output.write('HPD68\t%e\t%e\n' % (hpd68[0], hpd68[1]))
        output.write('HPD95\t%e\t%e\n' % (hpd95[0], hpd95[1]))
        output.write('MaxLike\t%e\t%e\t%e\n' % (ml, m, p))
        output.write('Log10 Maxlike\t%e\t%e\t%e\n' % (lml, lm, lp))
        output.close()

    with open('%s.mass.summary.pkl' % outputFile, 'wb') as output:
        stats = {'mean' : mean,
                 'stddev' : stddev,
                 'quantiles' : quantiles,
                 'hpd68' : hpd68,
                 'hpd95' : hpd95,
                 'maxlike' : (ml, (m, p)),
                 'log10maxlike' : (lml, (lm, lp))}
        cPickle.dump(stats, output)
        output.close()

    print 'mean\t%e' % mean
    print 'stddev\t%e' % stddev
    print 'Q2.5\t%e' % quantiles[2.5]
    print 'Q25\t%e' % quantiles[25]
    print 'Q50\t%e' % quantiles[50]
    print 'Q75\t%e' % quantiles[75]
    print 'Q97.5\t%e' % quantiles[97.5]
    print 'HPD68\t%e\t%e' % (hpd68[0], hpd68[1])
    print 'HPD95\t%e\t%e' % (hpd95[0], hpd95[1])
    print 'MaxLike\t%e\t%e\t%e\n' % (ml, m, p)
    print 'Log10 Maxlike\t%e\t%e\t%e\n' % (lml, lm, lp)



