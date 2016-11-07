########
# Contains All info for running fits assuming
# a simple Normal shapedistro
#########

import numpy as np
import maxlike_masses as maxlike, nfwmodeltools

class NormShapedistro(maxlike.LensingModel):

    ########

    def addCLOps(self, parser):
        
        super(NormShapedistro, self).addCLOps(parser)

        parser.add_option('--sigma', dest = 'sigma',
                          help = 'With of assumed Gaussian Distribution',
                          default = 0.25, type = 'float')

    ########

    def createOptions(self, sigma = 0.25, *args, **keywords):
        
        options, args = super(NormShapedistro, self).createOptions(*args, **keywords)

        options.sigma = sigma

        return options, args

    ########

        

    def __init__(self):

        maxlike.LensingModel.__init__(self)

        self.likelihood_func = nfwmodeltools.dumbnormal_like

    #######




    ####

    def bin_selectors(self, cat):
        return [np.ones(len(cat)) == 1]


    ###

    def makeShapePrior(self, manager, parts):
        
        self.shapedistro_params = [manager.options.sigma * np.ones(1)]
        
        super(NormShapedistro, self).makeFixedPrior(manager, parts)

