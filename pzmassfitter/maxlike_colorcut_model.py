###############################
# Implement overhead for running a single colorcut cluster
# in the maxlike framework, with mcmc sampling
###############################

import maxlike_masses as mm
import numpy as np
import nfwmodeltools as tools
import maxlike_benstep_voigt as bentvoigt
import histogramProb


#############################




class ColorcutModel(mm.LensingModel):

    def __init__(self):

        self.cuts = [self.modelCut]

    #########

    def modelCut(self, manager):

        options = manager.options
        inputcat = manager.inputcat

        minMPC = options.radlow   
        maxMPC = options.radhigh  

        goodObjs = np.logical_and(manager.inputcat['r_mpc'] > minMPC, 
                                  manager.inputcat['r_mpc'] < maxMPC)

        
        return goodObjs


    ############################

    def addCLOps(self, parser):

        parser.add_option('--masslow', dest='masslow',
                          help = 'Mass prior low cutoff',
                          type = 'float', default = 1e13)
        parser.add_option('--masshigh', dest='masshigh',
                          help = 'Mass prior high cutoff',
                          type = 'float', default = 1e16)
        parser.add_option('--radlow', dest='radlow',
                          help = 'Low radius cutoff wrt cluster center, Mpc',
                          default = 0.75, type = 'float')
        parser.add_option('--radhigh', dest='radhigh',
                          help = 'High radius cutoff wrt cluster center, Mpc',
                          default = 3.0, type='float')
        parser.add_option('--nshells', dest='nshells',
                          help = 'Number of radial bins in shear profile',
                          default = 12, type='int')



    ############################

    def createOptions(self, options = None, args = None, masslow = 1e13, masshigh = 1e16, radlow = 0.75, radhigh = 3.0,
                      nshells = 12):

        if options is None:
            options = varcontainer.VarContainer()

        options.masslow = masslow
        options.masshigh = masshigh        
        options.radlow = radlow
        options.radhigh = radhigh
        options.nshells = nshells


        return options, args

    #############################

    def makeShapePrior(self, datamanager, parts):


        parts.step_m_delta = pymc.MvNormalCov('step_m_delta', [0., 0., 0.], bentvoigt.bentvoigt3cov)
        parts.step_c_delta = pymc.Normal('step_c_delta', 0., 1./(0.0004**2))

        
        inputcat = datamanager.inputcat

        pivot, m_slope, m_b, m_cov, c = bentvoigt.bentvoigt3PsfDependence('interp',
                                                                          datamanager.psfsize)

        @pymc.deterministic(trace = False, name='shearcal_m')
        def shearcal_m(size = inputcat['size'], 
                       pivot = pivot, m_slope = m_slope, m_b = m_b, 
                       m_delta = parts.step_m_delta):


            pivot = pivot + m_delta[0]
            m_b = m_b + m_delta[1]
            m_slope = m_slope + m_delta[2]


            m = np.zeros_like(size)
            m[size >= pivot] = m_b
            m[size < pivot] = m_slope*(size[size < pivot] - pivot) + m_b

            return np.ascontiguousarray(m.astype(np.float64))

        parts.shearcal_m = shearcal_m


        @pymc.deterministic(trace = False, name='shearcal_c')
        def shearcal_c(size = inputcat['size'], c = c, cdelta = parts.step_c_delta):
            c = (c + cdelta)*np.ones_like(size)

            return np.ascontiguousarray(c.astype(np.float64))

        parts.shearcal_c = shearcal_c


    ##############################

    def makeModelPrior(datamanager, parts):

        super(ColorcutModel, self).makeModelPrior(datamanager, parts)

        parts.beta = histogramProb.HistoDistribution(data = datamanager.beta_samples, name = 'beta', bins = (50,50))

        parts.contamfrac500 = pymc.Normal('contamfrac500', datamanager.contamfrac[0], 1./datamanager.contamfrac[1]**2)


    ###############################

    def assignBins(self, r_mpc, options):

        numobjs = len(r_mpc)

        radial_order = np.argsort(r_mpc)

        assignment = np.ones(numobjs, dtype=np.int)

        ngals_pershell = numobjs / options.nshells

        curBin = 0
        for i in range(0, numobjs, ngals_pershell):

            maxTake = np.min(numobjs, i+ngals_pershell)
            assignment[radial_order[i:maxTake]] = curBin

            curBin += 1

        return assignment, curBin

    #######

    def calcBinRmpc(self, obj_r_mpc, obj_weights, binAssignments, nbins):

        bin_rmpc = np.zeros(nbins)

        for i in range(nbins):
            inBin = binAssignments == i
            bin_rmpc[i] = np.sum(obj_r_mpc[inBin]*obj_weights[inBin])/np.sum(obj_weights[inBin])

        return bin_rmpc

    ########

    def calcBinSigmas(self, obj_ghats, obj_weights, binAssignments, nbins):

        nbootstraps = 1000

        bin_sigmas = np.zeros(nbins)
        bootstrapmeans = np.zeros(nbootstraps)        

        for i in range(nbins):
            inBin = binAssignments == i


            gInBin = obj_ghats[inBin]
            weightsInBin = obj_weights[inBin]

            for j in range(nbootstraps):
                bootstap = np.random.randint(0, len(gInBin), len(gInBin))
                bootstrapmeans[j] = np.sum(gInBin[bootstrap]*weightsInBin[bootstrap])/np.sum(weightsInBin[bootstrap])

            bin_sigmas[i] = np.std(bootstrapmeans)

        return bin_sigmas
            


    ########

    def makeLikelihood(manager, parts):

        inputcat = datamanager.inputcat
        options = datamanager.options

        parts.obj_r_mpc = np.ascontiguousarray(inputcat['r_mpc'].astype(np.float64))
        parts.obj_ghats = np.ascontiguousarray(inputcat['ghats'].astype(np.float64))
        parts.obj_weights = np.ascontiguousarray(inputcat['sigma_g'].astype(np.float64))

        parts.binAssignments, parts.nBins = self.assignBins(parts.obj_r_mpc, options)

        parts.bin_r_mpc = self.calcBinRmpc(parts.obj_r_mpc, parts.obj_weights, parts.binAssignments, parts.nbins)
        parts.bin_sigmas = self.calcBinSigmas(parts.obj_ghats, parts.obj_weights, parts.binAssignments, parts.nbins)




        @pymc.deterministic(trace=False):
        def bin_ghats(binAssignments = parts.binAssignments,
                     nbins = parts.nBins,
                    obj_ghats = parts.obj_ghats, 
                    obj_weights = parts.obj_weights,
                    m = parts.shearcal_m,
                    c = parts.shearcal_c):

            gcorr = (obj_ghats - c)/(1+m)
            bin_ghats = np.zeros(nbins)
            for i in range(nbins):
                inBin = binAssignments == i
                bin_ghats[i] = np.sum(gcorr[inBin]*obj_weights[inBin])/np.sum(obj_weights[inBin])

            return bin_ghats

        parts.bin_ghats = bin_ghats


        @pymc.deterministic(trace=False):
        def bin_boostcor(f500 = parts.contamfrac500
                     r_mpc = parts.bin_r_mpc,
                     r500 = datamanager.r500):
            
            x = r_mpc / r500
            
            return 1./(1-f500*np.exp(1-x))

        parts.bin_boostcor = bin_boostcor

        

        parts.data = None
        for i in range(10):

            try:

                @pymc.stochastic(observed=True, name='data_%d' % i)
                def data(value = 0.,
                         r_mpc = parts.bin_r_mpc,
                         r_scale = parts.r_scale,
                         ghats = parts.bin_ghats,
                         sigma_ghat = parts.bin_sigmas,
                         concentration = parts.concentration,
                         zcluster = parts.zcluster,
                         beta = parts.beta,
                         boostcor = parts.bin_boostcor):

                    
                    avebeta, avebeta2 = beta
                    
                    return nfwtools.bentcolorcut_like(r_mpc,r_scale, ghats, sigma_ghat,concentration,
                                                      zcluster, avebeta, avebeta2, boostcor)






                parts.data = data


                break
            except pymc.ZeroProbability:
                pass

        if parts.data is None:
            raise ModelInitException




    #######
