##########
# For mass recon using Voigt model
############

import cPickle, numpy as np, pymc
import nfwmodeltools, maxlike_masses as maxlike

##############

__default_samplefile__ = '/nfs/slac/g/ki/ki06/anja/SUBARU/shapedistro/voigt_posterior.normapprox.pkl'


class VoigtnormShapedistro(object):

    def __init__(self, samplefile = __default_samplefile__):

        self.samplefile = samplefile

        self.likelihood_func = nfwmodeltools.voigt_like

        self.makeShapePrior = maxlike.makeNormPrior

        self.loadDistroFile()

    #####

    def loadDistroFile(self):

        input = open(self.samplefile, 'rb')

        self.shape_params = cPickle.load(input)

        input.close()

#####


    def bin_selectors(self, cat):

        bins = [cat['snratio'] < 4, 
                np.logical_and(cat['snratio'] >= 4, cat['snratio'] < 5),
                np.logical_and(cat['snratio'] >= 5, cat['snratio'] < 6),
                np.logical_and(cat['snratio'] >= 6, cat['snratio'] < 8),
                np.logical_and(cat['snratio'] >= 8, cat['snratio'] < 10),
                np.logical_and(cat['snratio'] >= 10, cat['snratio'] < 20),
                cat['snratio'] > 20]

        return bins

###

    def sampler_callback(self, mcmc):

        for shapeprior, covar in zip(mcmc.shape_params, [ x[1] for x in shape_params ]):
            mcmc.use_step_method(pymc.AdaptiveMetropolis, shapeprior, cov = covar, shrink_if_necessary=True)
