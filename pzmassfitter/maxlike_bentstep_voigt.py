#########################
# Implements a size dependent STEP correction in one bin
# correction is flat at large sizes, and linear at small sizes
########################


import maxlike_masses as mm, nfwmodeltools as nfwtools
import numpy as np, pymc, nfwutils

#######################


class BentVoigtShapedistro(mm.LensingModel):

    def addCLOps(self, parser):

        super(BentVoigtShapedistro, self).addCLOps(parser)

        parser.add_option('--steppsf', dest='steppsf',
                          help = 'Which STEP PSF to use, A,C,D,interp?',
                          default = 'interp')

    def createOptions(self, steppsf = 'interp', *args, **keywords):

        options, args = super(BentVoigtShapedistro, self).createOptions(*args, **keywords)

        options.steppsf = steppsf

        return options, args

    def psfDependence(self, steppsf, psfSize):

        if (steppsf == 'interp' and psfSize >= 1.87) or steppsf == 'c':
            m_slope = 0.254
            m_b = 0.050
            c = 0.00089
        elif (steppsf == 'interp' and psfSize <= 1.53) or steppsf == 'a':
            m_slope = 0.212
            m_b = -0.011
            c = 0.00042
        elif steppsf == 'interp' or steppsf == 'linear':
            m_slope = 0.12*(psfSize - 1.53) + 0.212
            m_b = 0.18*(psfSize - 1.53) - 0.011
            c = 0.0014*(psfSize - 1.53) + 0.00042
        elif steppsf == 'd':
            m_slope = 0.255
            m_b = 0.0051
            c = 0.0015

        if steppsf == 'd':
            m_cov = np.array([[ 0.00068239,  0.00119369], #m_b, m_slope
                              [ 0.00119369,  0.00306567]])
        else:
            m_cov = np.array([[ 0.0008,  0.0014],   #m_b, m_slope
                              [ 0.0014,  0.0033]])


        return m_slope, m_b, m_cov, c
        

    def makeShapePrior(self, data, parts):

        inputcat = data.inputcat

#        psfSize = data.psfsize # rh, in pixels

        #disabled since we don't have a STEP calibration for the regauss pipeline in DMSTACK
#        m_slope, m_b, m_cov, c = self.psfDependence(data.options.steppsf, psfSize)
#
#        
#        parts.step_m_prior = pymc.MvNormalCov('step_m_prior', [m_b, m_slope], m_cov)
#
#        @pymc.deterministic(trace = False)
#        def shearcal_m(size = inputcat['size'], mprior = parts.step_m_prior):
#
#            m_b = mprior[0]
#            m_slope = mprior[1]
#                           
#
#            m = np.zeros_like(size)
#            m[size >= 2.0] = m_b
#            m[size < 2.0] = m_slope*(size[size < 2.0] - 2.0) +m_b
#            
#            return np.ascontiguousarray(m.astype(np.float64))
#
#        parts.shearcal_m = shearcal_m
#
#        parts.step_c_prior = pymc.Normal('step_c_prior', c, 1./(0.0004**2))
#
#        @pymc.deterministic(trace = False)
#        def shearcal_c(size = inputcat['size'], cprior = parts.step_c_prior):
#            c = cprior*np.ones_like(size)
#            
#            return np.ascontiguousarray(c.astype(np.float64))
#
#        parts.shearcal_c = shearcal_c
#

        parts.shearcal_m = np.ones(len(inputcat))
        parts.shearcal_c = np.zeros(len(inputcat))


# turned off sampling of the wings for now
#        parts.sigma = pymc.Uniform('sigma', 0.15, 0.5) #sigma
#        parts.gamma = pymc.Uniform('gamma', 0.003, 0.1) #gamma

        parts.sigma = 0.2
        parts.gamma = 1.4e-2


    ####################

    def sampler_callback(self, mcmc):

        mcmc.use_step_method(pymc.AdaptiveMetropolis, 
                             [mcmc.sigma, mcmc.gamma], 
                             shrink_if_necessary = True)

    ####################

    def makeLikelihood(self, datamanager, parts):

        inputcat = datamanager.inputcat

        pz = datamanager.pz

        parts.r_mpc = np.ascontiguousarray(inputcat['r_mpc'].astype(np.float64))
        parts.ghats = np.ascontiguousarray(inputcat['ghats'].astype(np.float64))
        parts.pz = np.ascontiguousarray(pz.astype(np.float64))



        parts.zs = np.ascontiguousarray(np.array(datamanager.pdzrange).astype(np.float64))


        parts.betas = np.ascontiguousarray(nfwutils.global_cosmology.beta_s(parts.zs, parts.zcluster).astype(np.float64))
        parts.nzbins = len(parts.betas)

        parts.rho_c = nfwutils.global_cosmology.rho_crit(parts.zcluster)
        parts.rho_c_over_sigma_c = 1.5 * nfwutils.global_cosmology.angulardist(parts.zcluster) * nfwutils.global_cosmology.beta([1e6], parts.zcluster)[0] * nfwutils.global_cosmology.hubble2(parts.zcluster) / nfwutils.global_cosmology.v_c**2




        parts.data = None
        for i in range(10):

            try:

                @pymc.stochastic(observed=True, name='data_%d' % i)
                def data(value = parts.ghats,
                         mdelta = parts.mdelta,
                         cdelta = parts.cdelta,
                         r_mpc = parts.r_mpc,
                         zs = parts.zs,
                         betas = parts.betas,
                         pz = parts.pz, 
                         shearcal_m = parts.shearcal_m,
                         shearcal_c = parts.shearcal_c,
                         sigma = parts.sigma,
                         gamma = parts.gamma,
                         rho_c = parts.rho_c,
                         rho_c_over_sigma_c = parts.rho_c_over_sigma_c,
                         massdelta = parts.massdelta):                    
                    
                    return nfwtools.bentvoigt_like(mdelta,
                                                   cdelta,
                                                   r_mpc,
                                                   value,
                                                   zs, 
                                                   betas, 
                                                   pz, 
                                                   shearcal_m,
                                                   shearcal_c,
                                                   sigma,
                                                   gamma,
                                                   rho_c,
                                                   rho_c_over_sigma_c,
                                                   massdelta)





                parts.data = data


                break
            except pymc.ZeroProbability:
                pass

        if parts.data is None:
            raise ModelInitException



#########################################

bentvoigt3cov = np.array([[ 0.05106948,  0.00287574, -0.0167086 ],
                       [ 0.00287574,  0.00059742, -0.00013534],
                       [-0.0167086 , -0.00013534,  0.01099542]])


def bentvoigt3PsfDependence(steppsf, psfSize):


    if (steppsf == 'interp' and psfSize >= 1.87) or steppsf == 'c':
        pivot = 1.93
        m_b = 0.012
        m_slope = 0.22
        c = 0.0006

    elif (steppsf == 'interp' and psfSize <= 1.53) or steppsf == 'a':
        pivot = 1.97
        m_b = -.028
        m_slope = 0.20
        c = -9.6e-06

    elif steppsf == 'd':
        pivot = 2.08
        m_b = 9.5e-5
        m_slope = 0.23
        c = 0.00107



    elif steppsf == 'interp' or steppsf == 'linear':
        delta = psfSize - 1.53  #reference to psfA
        pivot = delta*(-0.118) + 1.97
        m_b = delta*(0.118) - 0.028
        m_slope = 0.059*delta + 0.2
        c = 0.0018*delta - 9.6e-06

# pivot, m_b, m_slope for dec 7 from PSF Afg
    m_cov = bentvoigt3cov

            

    return pivot, m_slope, m_b, m_cov, c




###########################################


class BentVoigt3Shapedistro(BentVoigtShapedistro):

            

        

    def makeShapePrior(self, data, parts):

        inputcat = data.inputcat

        psfSize = data.psfsize

        pivot, m_slope, m_b, m_cov, c = bentvoigt3PsfDependence(data.options.steppsf, psfSize)

        parts.step_m_prior = pymc.MvNormalCov('step_m_prior', [pivot, m_b, m_slope], m_cov)

        @pymc.deterministic(trace = False)
        def shearcal_m(size = inputcat['size'], mprior = parts.step_m_prior):

            pivot = mprior[0]
            m_b = mprior[1]
            m_slope = mprior[2]

                           

            m = np.zeros_like(size)
            m[size >= pivot] = m_b
            m[size < pivot] = m_slope*(size[size < pivot] - pivot) + m_b
            
            return np.ascontiguousarray(m.astype(np.float64))

        parts.shearcal_m = shearcal_m

        parts.step_c_prior = pymc.Normal('step_c_prior', c, 1./(0.0004**2))

        @pymc.deterministic(trace = False)
        def shearcal_c(size = inputcat['size'], cprior = parts.step_c_prior):
            c = cprior*np.ones_like(size)
            
            return np.ascontiguousarray(c.astype(np.float64))

        parts.shearcal_c = shearcal_c


        
        parts.sigma = pymc.Uniform('sigma', 0.15, 0.5) #sigma
        parts.gamma = pymc.Uniform('gamma', 0.003, 0.1) #gamma


