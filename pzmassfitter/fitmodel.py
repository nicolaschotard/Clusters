###################################
# Utilities for fitting models
# 
# Based on solution by abeardmore found on http://code.google.com/p/pyminuit/issues/detail?id=6
#
# Modified and extended by Douglas Applegate
###################################

import numpy
import minuit
import math, inspect
import scipy.stats as stats

###############################

__cvs_id__ = "$Id: fitmodel.py,v 1.2 2010-07-02 23:08:47 dapple Exp $"

###############################

###############################
# Statistical Distribution Look-up functions
###############################

def chisq_exceeds_prob(chisq, dof):
    '''
    Probability that chisq exceeds value, given degrees of freedom dof
    '''
    return stats.chi2.sf(chisq, dof)

###

def f_test_exceeds_prob(chisq_old, dof_old, chisq_new, dof_new):
    '''
    Probability that the improvement in a fit by adding extra parameters is random
    '''

    deltaDOF = dof_old - dof_new

    F = (chisq_old - chisq_new)/(deltaDOF*chisq_new/dof_new)

    return stats.f.sf(F, deltaDOF, dof_new)

###############################
# Common Models
###############################

def ConstantModel(x, a0):
    
    return a0

#######

def LinearModel(x, a0, a1):

    return a0 + a1*x

########

def QuadraticModel(x, a0, a1, a2):

    return a0 + a1*x + a2*x**2

########

class PolynomialModel(object):
    """
    Creates a polynomial model of the form
    a0 + a1*x + a2*x**2 + ... 
    where the order parameter controls which orders are included
    """

    def __init__(self, order):
        '''
        order is a list of positive integers specifying polynomial order to include
        0: constant, 1: linear, 2: quadratic, etc.
        Does not include lower order terms implicitly (ie specify [0,1,2], etc
        '''

        self.order = order

        self.basis = {}
        for o in order:
            param = 'a%d' % o
            def base(x, a, order=o):
                return a*(x**order)
            self.basis[param] = base

        self.params=self.basis.keys()

    def __call__(self, x, *params, **keyword_params):

        for key, val in zip(self.params, params):
            keyword_params[key] = val

        sum = 0.
        for key, val in keyword_params.iteritems():
            sum += self.basis[key](x, val)

        return sum

###########

def PowerLawModel(x, alpha, beta):

    return alpha*x**beta

###########

def GaussianModel(x, A, mu, sigma):

    z = (x - mu) / sigma
    return A*numpy.exp(-0.5*z**2)

###############################
# Statistical Fuctions for Minimization
###############################

def ChiSqStat(ydata, yerr, ymodel):
    """
    Returns the chi-square given arrays of ydata, yerr, and ymodel values.
    """
    chisquared = ((ydata - ymodel)/yerr)**2 
    stat = chisquared.sum()
    return stat

####################

def CStat(ydata, yerr, ymodel):
    """
    Returns the cstat a la xspec given arrays of data and model values.
    This is a -2.0 log likelihood statistic.
    """

    lmodel = numpy.zeros(ymodel.size)

    lmodel[ymodel <= 0.0] = -32.

    lmodel[ymodel > 0.0] = numpy.log(ymodel[ymodel > 0.0])

    ldata = numpy.zeros(ydata.size)

    ldata[ydata <= 0.0] = -32.0

    ldata[ydata > 0.0] = numpy.log(ydata[ydata > 0.0])

    # fitstat = ymodel - ydata  + ydata * (ldata - lmodel)

    fitstat = ymodel + ydata  * ((ldata - lmodel) - 1.0)

    stat = 2.0* fitstat.sum()

    return stat



###############################
# Fitting Class -- Use to perform minimizations
###############################

class FitModel:
    """
    Fits a generic model (provided by the class Model to data (numpy arrays
    xdata and ydata), with a fit statistic provided by StatFunc.
    """
    def __init__(self, xdata, ydata, yerr, model, 
                 statfunc = ChiSqStat, guess = []):

        self.xdata = numpy.array(xdata, dtype=numpy.float64)
        self.ydata = numpy.array(ydata, dtype=numpy.float64)
        self.yerr = numpy.array(yerr, dtype=numpy.float64)

        self.model = model
        
        self.statfunc = statfunc

        self.guess = guess

        self.fcn = FCN(self.xdata, self.ydata, self.yerr, model, statfunc)

        self.m = minuit.Minuit( self.fcn )

        self.params = self.m.parameters

        if self.guess == []:
            self.guess = numpy.ones(len(self.params))

        for param, value in zip(self.params, self.guess):
            self.m.values[param] = value
            self.m.errors[param] = math.fabs(value) * 0.05

        self.m.strategy = 1
        self.m.tol = 1.0

        self.have_fit = False

    def fixed(self, fparams):
        """
        Fix or unfix the parameters specified in the dictionary fparams, which
        contain True or False values.
        """
        for key in fparams.keys():
            self.m.fixed[key] = fparams[key]

    def limits(self, lparams):
        """
        Set limits given by the parameters in the dictionary lparams.
        """
        for key in lparams.keys():
            self.m.limits[key] = lparams[key]

    def fit(self, printmode = 0):
        """
        Call migrad to fit the model to the data.
        Set printmode = 1 to monitor the progress of the fitting.
        """

        self.m.printMode = printmode

        self.par_vals = {}

        self.ymodel = None

        try :

            self.m.migrad()

            print "fval = %g, nfcn %d" % (self.m.fval, self.m.ncalls)

            self.m.migrad()

            print "fval = %g, nfcn %d" % (self.m.fval, self.m.ncalls)

            print "Fit parameters : "
            print self.m.values

            self.par_vals = self.m.values

            # calculate the best fit model
            self.ymodel = self.model( self.xdata, **self.m.values )

            self.statval = self.m.fval

            self.have_fit = True

        except minuit.MinuitError :

            # reset have_fit if migrad fails
            self.have_fit = False



    def uncert(self, nsigma = 1.0):
        """
        Calculate the parameter uncertainties at the nsigma**2
        confidence level. E.g. for one parameter of interest
        nsigma = 1.0   for 68%
                 1.645 for 90%
                 2.0   for 95.45%
                 3.0   for 99.73%
        """

        if not(self.have_fit) :
            print "Warning: uncert requires a valid fit."
            return

        # in case minos fails
        self.m.hesse()

        print "Hesse errors : "
        print self.m.errors

        self.par_err = {}

        for key in self.m.values.keys():

            if (self.m.fixed[key] == True):
                continue

            try:

                self.m.minos(key, -nsigma)
                self.m.minos(key,  nsigma)


                error = (self.m.merrors[key, -nsigma], 
                         self.m.merrors[key,  nsigma])

            except minuit.MinuitError :

                print "Caught MinuitError: Minos failed. using Hesse error."
                print "Only really valid for a well behaved fitting FCN !"
                error = self.m.errors[key] * nsigma
    

            self.par_err[key] = error

            
        print "Parameter errors :"
        print self.par_err



    def corr_matrix(self):
        """
        Display the fit parameter correlation matrix."
        """
        if not(self.have_fit) :
            print "Warning: uncert requires a valid fit."
            return

        print "Correlation matrix :"
        print numpy.array(self.m.matrix(correlation=True))

    



#####################################
# Utilities
####################################

    
def FCN(x,y,yerr, model, statfunc):
    """
    Calculates the fitting FCN for pyMinuit(2) given the data (xdata & ydata)
    and model (class Model, with a tuple of initial parameters, params),
    using the class StatFunc to calculate the statistic.
    """

    #assumes model is a function with first arg being X values
    if inspect.isfunction(model):
        params = inspect.getargspec(model)[0][1:]
    elif hasattr(model, '__call__'):
        args = inspect.getargspec(model.__call__)[0]
        if len(args) < 3:
            paramAttr = inspect.getargspec(model.__call__)[1]
            params = getattr(model, paramAttr)
        else:
            params = args[2:]

    paramstring = ','.join(params)

    class_template = '''class fitclass(object):
    def __init__(self, x, y, yerr, model, statfunc):
        self.x = x
        self.y = y
        self.yerr = yerr
        self.model = model
        self.statfunc = statfunc

    def __call__(self, %s):
        return self.statfunc(self.y, self.y, self.model(self.x, %s))
''' % (paramstring, paramstring)

    
    exec class_template


    return fitclass(x,y,yerr,model,statfunc)


