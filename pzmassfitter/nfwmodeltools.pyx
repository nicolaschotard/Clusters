##########################
# Implements an NFW model for investigation 
# of redshift and contamination effects
###########################
# Compiling info: gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I /u/ki/dapple/include/python2.7/ -I /u/ki/dapple/lib/python2.7/site-packages/numpy/core/include/ -o nfwmodeltools.so nfwmodeltools.c voigt.c
# Compiling at AIFA: gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I /users/dapple/anaconda/pkgs/numpy-1.7.1-py27_2/lib/python2.7/site-packages/numpy/core/include/ -I /users/dapple/anaconda/include/python2.7/ -o nfwmodeltools.so nfwmodeltools.c


# cython: profile=False

import numpy as np
cimport numpy as np
cimport cython

import scipy.optimize

import nfwutils

cdef extern from "math.h":
    double exp(double)
    double log(double)
    double atanh(double)
    double sqrt(double)
    double atan(double)

cdef extern from "voigt.h":
    double voigt(double, double, double)



########################

__cvs_id__ = "$Id: nfwmodeltools.pyx,v 1.5 2011-02-09 01:59:14 dapple Exp $"

#########################

logsqrt2pi = np.log(np.sqrt(2*np.pi))

DTYPE = np.double
ctypedef np.double_t DTYPE_T

#########################


############################
# NFW Profile
############################

cdef double deltaC(double c, double delta = 200.):
    return (delta/3.) * c**3 / (log(1+c) - c/(1+c))

##############


@cython.boundscheck(False)
@cython.wraparound(False)
def NFWShear(np.ndarray[np.double_t, ndim=1, mode='c'] r, 
              double concentration, 
              double rs,
              double rho_c_over_sigma_c, 
             double delta = 200.):
    
    cdef double delta_c = deltaC(concentration, delta = delta)
    cdef double amp = rs*delta_c*rho_c_over_sigma_c

    cdef double x,a,b,c
    cdef Py_ssize_t i, npos
    npos = r.shape[0]
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] g = np.zeros(r.shape[0], dtype=np.float64)

    for i from npos > i >= 0:
        
        x = r[i]/rs

        if x < 1:
            
            a = atanh(sqrt((1-x)/(1+x)))
            b = sqrt(1-x**2)
            c = (x**2) - 1
    
            g[i] = 8*a/(b*x**2) + 4*log(x/2)/x**2 - 2/c + 4*a/(b*c)

        elif x > 1:

            a = atan(sqrt((x-1)/(1+x)))
            b = sqrt(x**2-1)
    
            g[i] = 8*a/(b*x**2) + 4*log(x/2)/x**2 - 2/b**2 + 4*a/b**3

        else:

            g[i] = 10./3 + 4*log(.5)

    return amp*g

###################

@cython.boundscheck(False)
@cython.wraparound(False)
def NFWKappa(np.ndarray[np.double_t, ndim=1, mode='c'] r, 
              double concentration, 
              double rs,
              double rho_c_over_sigma_c,
             double delta = 200.):

    
    cdef double delta_c = deltaC(concentration, delta = delta)
    cdef double amp = 2*rs*delta_c*rho_c_over_sigma_c

    cdef Py_ssize_t i, npos
    npos = r.shape[0]
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] kappa = np.zeros(npos, dtype=np.float64)

    cdef double x, a,b,c

    for i from npos > i >= 0:

        x = r[i]/rs

        if x < 1:
            
            a = atanh(sqrt((1-x)/(1+x)))
            b = sqrt(1-x**2)
            c = 1./(x**2 - 1)
            kappa[i] = c*(1 - 2.*a/b)

        elif x > 1:
            a = atan(sqrt((x-1)/(1+x)))
            b = sqrt(x**2-1)
            c = 1./(x**2 - 1)
            kappa[i] = c*(1 - 2.*a/b)

        else:
            kappa[i] = 1./3.

    return kappa*amp


###############################

@cython.boundscheck(False)
@cython.wraparound(False)
def aveEnclosedKappa(np.ndarray[np.double_t, ndim=1, mode='c'] r, 
                     double concentration, 
                     double rs,
                     double rho_c_over_sigma_c,
                     double delta = 200.):

    cdef double delta_c = deltaC(concentration, delta = delta)
    cdef double amp = 4*rs*delta_c*rho_c_over_sigma_c

    cdef Py_ssize_t i, npos
    npos = r.shape[0]
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] avekappa = np.zeros(npos, dtype=np.float64)

    cdef double x, a,b,c

    for i from npos > i >= 0:

        x = r[i]/rs

        if x < 1:
            
            a = atanh(sqrt((1-x)/(1+x)))
            b = sqrt(1-x**2)
            c = log(x/2.)
            avekappa[i] = (2*a/b + c)/(x**2)

        elif x > 1:
            a = atan(sqrt((x-1)/(1+x)))
            b = sqrt(x**2-1)
            c = log(x/2.)
            avekappa[i] = (2*a/b + c)/(x**2)

        else:
            avekappa[i] = 1 + log(0.5)

    return avekappa*amp






###############################

@cython.boundscheck(False)
@cython.wraparound(False)
def rdelta2rs(double rdelta, 
              double c200, 
              double delta):

    cdef double delta_c = deltaC(c200)
    
    # x = r_delta / rs
    def f(x):
        
        return 3*delta_c*(log(1+x) - (x/(1+x)))/x**3 - delta

    
    x0 = scipy.optimize.brenth(f, 0.1, 20)

    cdef double rs = rdelta / x0

    return rs

#####################


@cython.boundscheck(False)
@cython.wraparound(False)
def rscaleConstM(double mdelta,
                 double c,
                 double rho_c, 
                 double delta):

    cdef double rdelta = (3*mdelta/(4*delta*np.pi*rho_c))**(1./3.)

    cdef double rs = rdelta2rs(rdelta, c, delta)

    return rs


#######################

@cython.boundscheck(False)
@cython.wraparound(False)
def massInsideR(double rs, 
                double c, 
                double R, 
                double rho_c):

    cdef double x = R/rs
    cdef double delta_c = deltaC(c)

    cdef double massInsideR = (log(1+x) - (x/(1+x)))*4*np.pi*delta_c*rho_c*rs**3


    return massInsideR



######################
# ML Reconstruction Tools
######################
#########



@cython.boundscheck(False)
@cython.wraparound(False)
def bentvoigt_like(double mdelta,
                   double cdelta,
                   np.ndarray[DTYPE_T, ndim=1, mode='c'] r_mpc not None, 
                   np.ndarray[DTYPE_T, ndim=1, mode='c'] ghats not None,
                   np.ndarray[DTYPE_T, ndim=1, mode='c'] zs not None, 
                   np.ndarray[DTYPE_T, ndim=1, mode='c'] betas not None, 
                   np.ndarray[DTYPE_T, ndim=2, mode='c'] pz not None, 
                   np.ndarray[DTYPE_T, ndim=1, mode='c'] m not None,
                   np.ndarray[DTYPE_T, ndim=1, mode='c'] c not None,
                   double sigma,
                   double gamma,
                   double rho_c,
                   double rho_c_over_sigma_c,
                   double massdelta):


    cdef Py_ssize_t nobjs = pz.shape[0]
    cdef Py_ssize_t npz = pz.shape[1]



    cdef np.ndarray[DTYPE_T, ndim=1, mode='c'] gamma_inf
    cdef np.ndarray[DTYPE_T, ndim=1, mode='c'] kappa_inf
    cdef double rscale

    if mdelta == 0:
        gamma_inf = np.zeros(nobjs)
        kappa_inf = np.zeros(nobjs)
    else:

        rdelta = (3*abs(mdelta)/(4*massdelta*np.pi*rho_c))**(1./3.)
        rscale = rdelta / cdelta

        gamma_inf = NFWShear(r_mpc, cdelta, rscale, rho_c_over_sigma_c, delta = massdelta)
        kappa_inf = NFWKappa(r_mpc, cdelta, rscale, rho_c_over_sigma_c, delta = massdelta)

    if mdelta < 0.:
        gamma_inf = -gamma_inf



    cdef double delta = 0.
    cdef double curPZ = 0.
    cdef double beta = 0.
    cdef double g = 0.

    cdef Py_ssize_t i, j, s
    cdef DTYPE_T galProb = 0.
    cdef DTYPE_T logProb = 0.
    cdef np.ndarray[DTYPE_T, ndim=1, mode='c'] integrand = np.zeros(npz)
        
    for i from nobjs > i >= 0:

        #zero out the array for the new object
        for j from npz > j >= 0:
            integrand[j] = 0.


        for j from npz > j >= 0:

            curPZ = pz[i,j]
            if curPZ > 1e-6:

                beta = betas[j]

                g = (beta*gamma_inf[i] / (1 - beta*kappa_inf[i]))


                delta = ghats[i] - (1+m[i])*g - c[i]

                integrand[j] = curPZ*voigt(delta, sigma, gamma)

        galProb = np.trapz(integrand, zs)
        logProb = logProb + log(galProb)

        
    return logProb


#############################################################


@cython.boundscheck(False)
@cython.wraparound(False)
def bentcolorcut_like(double r_scale, 
               np.ndarray[DTYPE_T, ndim=1, mode='c'] r_mpc not None, 
               np.ndarray[DTYPE_T, ndim=1, mode='c'] ghats not None,
               np.ndarray[DTYPE_T, ndim=1, mode='c'] sigma_ghat not None,
               double concentration, 
               double zcluster,
               double beta_s,
               double beta_s2,
               np.ndarray[DTYPE_T, ndim=1, mode='c'] boostcor not None):

    cdef Py_ssize_t nbins = r_mpc.shape[0]

    cdef np.ndarray[DTYPE_T, ndim=1, mode='c'] gamma_inf = NFWShear(r_mpc, concentration, r_scale, zcluster)
    cdef np.ndarray[DTYPE_T, ndim=1, mode='c'] kappa_inf = NFWKappa(r_mpc, concentration, r_scale, zcluster)


    cdef double g = 0.
    cdef double binlogProb = 0.
    cdef logProb = 0.

    cdef Py_ssize_t i
        
    for i from nbins > i >= 0:

        g = beta_s*gamma_inf / (1 - ((beta_s2/beta_s)*kappa_inf) )

        binlogProb = -0.5*((boostcor[i]*ghats[i] - g)/sigma_ghat[i])**2

        logProb = logProb + binlogProb

    return logProb

##################################################

@cython.boundscheck(False)
@cython.wraparound(False)
def bentcolorcut_like(np.ndarray[DTYPE_T, ndim=1, mode='c'] bin_r_mpc not None, 
                      double r_scale, 
                      np.ndarray[DTYPE_T, ndim=1, mode='c'] ghats not None, 
                      np.ndarray[DTYPE_T, ndim=1, mode='c'] bin_sigma not None,
                      double concentration,
                      double zcluster,
                      double avebeta,
                      double avebeta2,
                      np.ndarray[DTYPE_T, ndim=1, mode='c'] bin_contamboost not None):




    cdef Py_ssize_t nobjs = ghats.shape[0]
    cdef Py_ssize_t nbins = bin_r_mpc.shape[0]

    cdef np.ndarray[DTYPE_T, ndim=1, mode='c'] gamma_inf = NFWShear(bin_r_mpc, concentration, r_scale, zcluster)
    cdef np.ndarray[DTYPE_T, ndim=1, mode='c'] kappa_inf = NFWKappa(bin_r_mpc, concentration, r_scale, zcluster)


    cdef int curbin = 0
    cdef double gtilde = 0.
    cdef double delta = 0.
    cdef double modelg = 0.

    cdef Py_ssize_t i, j, s
    cdef DTYPE_T logProb = 0.

        
    #calculate logprob
    for i from nbins > i >= 0:

        gtilde = bin_contamboost[i] * ghats[i] 
       
        modelg = (avebeta*gamma_inf[i] / (1 - (avebeta2/avebeta)*kappa_inf[i]))

        delta = gtilde - modelg

        logProb = logProb -.5*(delta/bin_sigma[i])**2  - logsqrt2pi - np.log(bin_sigma[i])

        
    return logProb


