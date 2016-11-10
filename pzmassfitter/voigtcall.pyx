##########################
# Allows python calls to voigt function
#########################
# Compiling info: gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I /u/ki/dapple/include/python2.7/ -I /u/ki/dapple/lib/python2.7/site-packages/numpy/core/include/ -o voigtcall.so voigtcall.c voigt.c

# cython: profile=True

import numpy as np
cimport numpy as np
cimport cython

import nfwutils

cdef extern from "math.h":
    double exp(double)
    double log(double)

cdef extern from "voigt.h":
    double voigt(double, double, double)



@cython.boundscheck(False)
@cython.wraparound(False)
def voigtcall(np.ndarray[np.double_t, ndim=1, mode='c'] delta not None,
              double sigma,
	      double gamma):
	      
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] results = np.zeros_like(delta)
    cdef Py_ssize_t nobjs = delta.shape[0]
    cdef Py_ssize_t i

    for i from nobjs > i >= 0:
        results[i] = voigt(delta[i], sigma, gamma)

    return results


##########################################


@cython.boundscheck(False)
@cython.wraparound(False)
def interpvoigtcall(np.ndarray[np.double_t, ndim=1, mode='c'] delta not None,
              double sigma,
	      double gamma,
              np.ndarray[np.double_t, ndim=1, mode='c'] samplepoints = np.arange(-5,5,0.005)):

    cdef np.ndarray[np.double_t, ndim=1, mode='c'] voigtcalls = voigtcall(samplepoints, sigma, gamma)
    
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] interpolation = np.interp(delta, samplepoints, voigtcalls, -1, -1)

    cdef Py_ssize_t nobjs = delta.shape[0]
    cdef Py_ssize_t i
    cdef double cursample

    for i from nobjs > i >= 0:

        if interpolation[i] == -1:
            interpolation[i] = voigt(delta[i], sigma, gamma)


    return interpolation

############################################
###########################################
###########################################


@cython.boundscheck(False)
@cython.wraparound(False)
def gausscall(np.ndarray[np.double_t, ndim=1, mode='c'] delta not None,
              double sigma):
	      
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] results = np.zeros_like(delta)
    cdef Py_ssize_t nobjs = delta.shape[0]
    cdef Py_ssize_t i

    for i from nobjs > i >= 0:
        results[i] = exp(-.5*(delta[i]/sigma)**2)

    return results





############################################



@cython.boundscheck(False)
@cython.wraparound(False)
def interpgausscall(np.ndarray[np.double_t, ndim=1, mode='c'] delta not None,
              double sigma,
              np.ndarray[np.double_t, ndim=1, mode='c'] samplepoints = np.arange(-5,5,0.005)):

    cdef np.ndarray[np.double_t, ndim=1, mode='c'] gausscalls = gausscall(samplepoints, sigma)
    
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] interpolation = np.interp(delta, samplepoints, gausscalls, -1, -1)

    cdef Py_ssize_t nobjs = delta.shape[0]
    cdef Py_ssize_t i
    cdef double cursample

    for i from nobjs > i >= 0:

        if interpolation[i] == -1:
            interpolation[i] = exp(-.5*(delta[i]/sigma)**2)


    return interpolation




