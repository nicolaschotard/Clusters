# Copyright (C) 2011, 2012 Adam Mantz
"""A module for Markov Chain Monte Carlo. See help() and example()."""

import csv
import glob
import sys
import cPickle as pickle
import numpy as np
try:
    from mpi4py import MPI
except ImportError:
    pass
try:
    import numdifftools
    have_numdifftools = True
except ImportError:
    have_numdifftools = False
try:
    import scipy.optimize
    have_scipy = True
except ImportError:
    have_scipy = False

parallel_filename_base = 'mymc'
parallel_filename_ext = '.chk'


def help():
    print """
The module should be very flexible, but is designed with these things foremost in mind:
  1. use with expensive likelihood calculations which probably have a host of hard-to-modify
    code associated with them.
  2. making is straightforward to break the parameter space into subspaces which can be sampled
    using different proposal methods and at different rates. For example, if changing some parameters
    requires very expensive calulations in the likelihood, the other, faster parameters can be sampled
    at a higher rate. Or, some parameters may lend themselves to Gibbs sampling, while others may not,
    and these can be block updated independently.
  3. keeping the overhead low to facilitate large numbers of parameters. Some of this has been lost in
    the port from C++, but, for example, the package provides automatic tuning of the proposal
    covariance for block updating without needing to store traces of the parameters in memory.

Real-valued parameters are assumed, but the framework can be used with other types of parameters,
with suitable overloading of classes.

A byproduct of item (1) is that the user is expected to handle all aspects of the calculation of the
posterior. The module doesn't implement assignment of canned, standard priors, or automatic discovery
of shortcuts (like Gibbs sampling). The idea is that the user is in the best position to know how the 
details of the likelihood and priors should be implemented. (Note that, as is, the code very naturally
handles log-posterior values of -infinity, corresponding to zero probability, e.g. when a proposed
parameter value is outside the allowed range. The only point to be careful of is that the sampler
needs to start at a point with non-zero probability.)

It's probably worth reading the docstrings for:
 - Parameter
 - ParameterSpace
 - Updater, CartesianSequentialUpdater, CartesianPermutationUpdater, MultiDimSequentialUpdater, 
   MultiDimPermutationUpdater, emceeUpdater
 - Slice, Metropolis 
 - randNormalExp, randChiExp
 - textBackend, stdoutBackend, dictBackend
 - Engine


Here is a quick overview of the class structure:
 Parameter objects should be self explanatory.
 ParameterSpace objects represent sets of Parameters for organizational purposes.
 Updater objects provide the mechanism to sample a ParameterSpace. They come in three varieties:
   Cartesian updaters perform updates to each Parameter in their ParameterSpace individually.
   MultiDim updaters perform block updates to all parameters in their ParameterSpace at the same time.
   emceeUpdater is an interface to the emcee package, and is somewhat different than described below;
    see its docstring.
  Each of these comes in Sequential and Permutation flavors, corresponding to sampling each direction
    in the ParameterSpace in fixed or random order. There is also a Rotation version of the MultiDim
    updater, which proposes along completely random directions in the multi-dimension parameter space,
    with the step length scaled appropriately for the direction chosen, given the covariance matrix.
  Updaters can adapt their proposal distributions. Cartesian updaters tune the typical scale of
    proposed steps for each parameter, while MultiDim updaters tune the set of basis directions used
    to explore the ParameterSpace, as well as the scale along each direction.
  Updater states (e.g. the proposal basis after adapting) can be saved to and restored from files
    using the save() and restore() methods. The relevant data are stored as a dictionary using cPickle;
    this should be safer and preferable to pickling the Updater object directly. Cartesian and
    MultiDim updaters are mutually compatible as far as this functionality is concerned, inasmuch as 
    that's possible.
 Step objects implement the specific algorithm used to propose a step along a given direction. Slice
    and Metropolis algorithms are implemented. The distribution of proposal lengths by Metropolis Step
    objects is customizable.
 Backend objects handle the storage of Parameter values as the chain progresses.
 Engine objects hold a list of Updater objects, each of which is called in a single iteration of the
    chain.

Finally, the ChiSquareLikelihood class may simplify the setup of least-squares style problems. Similar
    helper classes could clearly be added.


Communication between parallel chains can significantly speed up convergence. In parallel mode,
    adaptive Updaters use information from all running chains to tune their proposals, rather than
    only from their own chain. The Gelman-Rubin convergence criterion (ratio of inter- to intra-chain
    variances) for each free parameter is also calculated. Parallelization is implemented in two ways;
    see ?Updater for instructions on using each.
  1. Via MPI (using mpi4py). MPI adaptations are synchronous: when a chain reaches a communication
    point, it stops until all chains have caught up. All Updaters in a given chain should use the
    same communicator (at least, there's no advantage in doing otherwise). (This could be made safely
    asynchronous using shared memory, although I'm not sure that would be an improvement.)
  2. Via the filesystem. When a chain adapts, it will write its covariance information to a file. It
    will then read in any information from other chains that is present in similar files, and
    incorporate it when tuning. This process is asynchronous; chains will not wait for one another,
    they will simply adapt using whatever information has been shared at the time. The global
    variables parallel_filename_base and parallel_filename_ext can be used to customize the prefix
    and suffix of the files written.


The definition of an MCMC iteration in this implementation can be a little confusing. As far as an
MultiDim updater object's 'count' attribute (controling the timing of adaptations, for example) is
 concerned, an iteration corresponds to one proposal along some direction. For Cartesian updaters,
 each iteration corresponds to a separate proposal for all of the Parameters in the Updater's
 ParameterSpace. (This disparity is unfortunate, but ensures that the count attribute is the correct
 value to use when, e.g., calculating chain vairnaces on the fly.) In either case, the Updater's rate
 attribute determines how many of these iterations occur each time the updater is called by the
 Engine (default=1). However, as far as the Engine is concerned, an iteration corresponds to a loop
 over all loaded Updaters. So if there are two MultiDim updaters, u1 and u2, with rates set to
 1 and 2, and the Engine's updater list is [u2, u1, u2] (repetition is perfectly allowed, and
 sometimes desirable), then each iteration of the Engine actually corresponds to
 1*u1.rate + 2*u2.rate = 5 proposals. However, if u1 is instead a Cartesian updater controlling
 3 Parameters, then each engine iteration implies 7 proposals. The chain is stored to the Backends
 after each Engine iteration. The upshot of all this is that there's lots of flexibility to sample
 parameters at different rates and/or thin the chain as it's being run, at the expense of a
 little complexity.

See also mymc.example()


Crude non-MCMC functionality:
1. ParameterSpace.optimize provides deterministic minimization (using scipy.optimize.fmin_powell)
 to get an decent initial fit. This is not very efficient.
2. Updater.set_covariance_from_hessian uses finite differencing to estimate an appropriate
 cartesian width in each direction. This will fail if the state is not in a local minimum,
 or just because.
"""


class Parameter(object):
    """
    Class to handle a single free parameter. Can/often should be overriden.

    The only critical attributes are:
     1. width: initial guess for step lengths. Adaptive CartesianUpdaters change this value.
     2. (): return the current parameter value.
     3. set( ): set the parameter to a new value.
    """
    def __init__(self, value=0.0, width=1.0, name=''):
        self.value = value
        self.width = width
        self.name = name

    def __call__(self):
        return self.value

    def __str__(self):
        return self.name + ': ' + str(self.value) + ' ' + str(self.width)

    def set(self, value):
        self.value = value


class DerivedParameter(Parameter):
    """
    Prototype of a Parameter whose value may not be sampled because it is a deterministic function
    of other parameters. The user is responsible for setting its value. This should be done directly
    rather than by using set(), for safety reasons.
    """
    def __init__(self, value=0.0, name=''):
        Parameter.__init__(self, value, None, name)

    def set(self, value):
        raise Exception('Attempt to set() value of a DerivedParameter.')


class postgetter(object):
    # needs space and struct to be defined
    # actually returns -2*loglike
    def __init__(self):
        self.verbose = False
        self.last = 1e300

    def __call__(self, x):
        for i, p in enumerate(self.space):
            self.space[i].set(x[i])
        chisq = -2.0 * self.space.log_posterior(self.struct)
        if np.isinf(chisq):
            chisq = 1e300
        if self.verbose and chisq < self.last:
            self.last = chisq
            print chisq, x
        return chisq

class ParameterSpace(list):
    """
    Class to define sets of parameters (parameter spaces); inherits list.

    To sample the parameter space, attribute log_posterior must be set to a function of one argument
    that evaluates the *complete* posterior likelihood, including priors and parameters not in
    this ParameterSpace.
    """
    def __init__(self, parameterList=None, log_posterior=None):
        if parameterList is None:
            parameterList = []
        list.__init__(self, parameterList)
        self.log_posterior = log_posterior

    def __str__(self):
        st = ''
        for p in self:
            st = st + '\n' + str(p)
        return st

    def optimize(self, struct, xtol=0.01, ftol=0.01, maxiter=10):
        if not have_scipy:
            print "ParameterSpace.optimize requires the scipy package -- aborting."
            return None
        g = postgetter()
        g.verbose = True
        g.space = self
        g.struct = struct
        origin = [p() for p in self]
        try:
            m = scipy.optimize.fmin_powell(g, origin, full_output=True,
                                           xtol=xtol, ftol=ftol, maxiter=maxiter)
            ret = -0.5 * m[1] # log-likelihood for best point
        except:
            print "ParameterSpace.optimize: warning -- some kind of error in scipy.optimize.fmin_powell."
            for i, p in enumerate(self):
                p.set(origin[i])
            ret = None
        for i, p in enumerate(self):
            p.set(m[0][i])
        #print m[2]
        return ret


class Updater(object):
    """
    Abstract base class for updaters. Do not instantiate directly.
    Constructor arguments:
     1* ParameterSpace to update.
     2* Step (proposal method) to use. This can be safely changed later using the set_step( ) method.
     3  Number of steps between proposal adaptations (zero to not adapt).
        For parallel chains, the Gelman-Rubin convergence criterion is calculated at the same interval.
     4  Number of steps before the first proposal adaptation.
     5  A function (with one arbitrary argument) to call after each adaptation.
     6  For no parallelization, set to None.
        For filesystem parallelization, set to a unique, scalar identifier.
        For MPI parallelization, set to an mpi4py.MPI.Comm object (e.g. MPI.COMM_WORLD)
        See module docstring for more details.
    """
    def __init__(self, space, step, every, start, on_adapt, parallel):
        self.space = space
        self.set_step(step)
        self.engine = None
        self.index = None
        self.adapt = (every > 0)
        self.adapt_every = every
        self.adapt_start = max(1, start)
        self.onAdapt = on_adapt
        self.count = 0
        if self.adapt:
            self.R = None
        self.rate = 1

    def restore(self, filename):
        f = open(filename, 'rb')
        s = pickle.load(f)
        self.restoreBits(s)
        f.close()

    def save(self, filename):
        f = open(filename, 'wb')
        pickle.dump(self.saveBits(), f)
        f.close()

    def set_step(self, step):
        if step is None:
            return
        # todo: prevent direct assignment bypassing this
        self.step = step
        self.step.updater = self        

class CartesianUpdater(Updater):
    """
    Abstract base class for updaters that proposal one parameter at a time.

    Do not instantiate directly.
    """
    def __init__(self, space, step, adapt_every, adapt_starting, on_adapt, parallel):
        Updater.__init__(self, space, step, adapt_every, adapt_starting, on_adapt, parallel)
        if self.adapt:
            self.means = np.zeros(len(self.space))
            self.variances = np.zeros(len(self.space))
            if parallel is not None:
                mpi = False
                try:
                    if isinstance(parallel, MPI.Comm):
                        mpi = True
                        self.gatherAdapt = self.gatherMPI
                        self.comm = parallel
                except NameError:
                    pass
                if not mpi:
                    self.pid = str(parallel)
                    self.uind = '_' + str(self.index)
                    self.gatherAdapt = self.gatherFilesys
            else:
                self.gatherAdapt = self.gatherSerial
        self.current_direction = 0
        self.origin = 0.0

    def __call__(self, struct):
        if self.adapt and self.count >= self.adapt_start and self.count % self.adapt_every == 0:
            self.do_adapt(struct)
        self.count += 1
        for j in range(len(self.space)):
            self.choose_direction(j)
            self.origin = self.space[self.current_direction]()
            self.step(struct)
            self.accumulate()

    def accumulate(self):
        if self.adapt:
            # Golub, Chan and Levesque one-pass mean and variance algorithm
            d = self.space[self.current_direction]() - self.means[self.current_direction]
            self.means[self.current_direction] += d / self.count
            # this is actually (n-1) times the variance (below)
            self.variances[self.current_direction] += (self.count-1.0) / self.count * d**2
    def do_adapt(self, struct):
        stdevs = self.gatherAdapt()
        for i, p in enumerate(self.space):
            if stdevs[i] != 0.0:
                p.width = stdevs[i]
        if not self.onAdapt is None:
            self.onAdapt(struct)

    def gatherFilesys(self):
        filename = parallel_filename_base + self.pid + self.uind + parallel_filename_ext
        self.save(filename)
        total = 0
        moment1 = np.zeros(len(self.space))
        moment2 = np.zeros(len(self.space))
        grandMeans = np.zeros(len(self.space))
        grandMeanVar = np.zeros(len(self.space))
        grandVarMean = np.zeros(len(self.space))
        j = 0
        for filename in glob.iglob(parallel_filename_base + '*' + self.uind + parallel_filename_ext):
            try:
                f = open(filename, 'rb')
                s = pickle.load(f)
                f.close()
                total += s['count'] # becomes Ntot
                moment1 += s['count'] * s['means'] # becomes Ntot*<x>
                moment2 += s['count'] * (s['variances'] / \
                                         (s['count']-1.0) + s['means']**2) # becomes Ntot*<x^2>
                d = s['means'] - grandMeans
                grandMeans += d / (j + 1.0)
                grandMeanVar += j / (j + 1.0) * d**2
                d = s['variances'] / (s['count']-1.0) - grandVarMean
                grandVarMean += d / (j + 1.0)
                j += 1
            except IOError:
                print "Warning: IO error while reading " + filename + " to update covariance (process " + self.pid + ", updater " + self.uind + ")."
        if j > 1:
            B = self.count / (j - 1.0) * grandMeanVar
            W = grandVarMean / j
#            self.R = np.sqrt((self.count-1.0)/self.count + B/(self.count*W))
        return np.sqrt((moment2 - moment1**2 / total) / (total - 1.0))

    def gatherMPI(self):
        alls = self.saveBits()
        alls = self.comm.allgather(alls)
        total = 0
        moment1 = np.zeros(len(self.space))
        moment2 = np.zeros(len(self.space))
        grandMeans = np.zeros(len(self.space))
        grandMeanVar = np.zeros(len(self.space))
        grandVarMean = np.zeros(len(self.space))
        for j, s in enumerate(alls):
            total += s['count'] # becomes Ntot
            moment1 += s['count'] * s['means'] # becomes Ntot*<x>
            moment2 += s['count'] * (s['variances'] / \
                                     (s['count']-1.0) + s['means']**2) # becomes Ntot*<x^2>
            d = s['means'] - grandMeans
            grandMeans += d / (j + 1.0)
            grandMeanVar += j / (j + 1.0) * d**2
            d = s['variances'] / (s['count']-1.0) - grandVarMean
            grandVarMean += d / (j + 1.0)
        if len(alls) > 1:
            B = self.count / (len(alls) - 1.0) * grandMeanVar
            W = grandVarMean / len(alls)
            self.R = np.sqrt((self.count-1.0) / self.count + B / (self.count*W))
        return np.sqrt((moment2 - moment1**2 / total) / (total - 1.0))

    def gatherSerial(self):
        return np.sqrt(self.variances / (self.count-1.0))

    def move(self, x):
        p = self.space[self.current_direction]
        p.set(self.origin + x * p.width)

    def restoreBits(self, s):
        self.count = s['count']
        if s['type'] == 'Cartesian':
            for i, p in enumerate(self.space):
                p.width = s['widths'][i]
            self.means = s['means']
            self.variances = s['variances']
        elif  s['type'] == 'MultiDim':
            for i, p in enumerate(self.space):
                p.width = np.sqrt(s['covariances'][i, i])
            self.means = s['means']
            self.variances = s['covariances'].diagonal()
        else:
            raise Exception('CartesianUpdater.restoreBits: error restoring updater state -- unknown updater type')

    def saveBits(self):
        if self.adapt:
            return {'type': 'Cartesian', 'count': self.count, 'means': self.means,
                    'variances': self.variances, 'widths': [p.width for p in self.space]}
        else:
            return None

    def scatter(self, struct, ntries=10):
        c = self.current_direction
        origin = [p() for p in self.space]
        for i in range(ntries):
            for self.current_direction in range(len(self.space)):
                self.origin = origin[self.current_direction]
                self.move(np.random.randn())
            self.engine.current_logP = self.space.log_posterior(struct)
            if self.engine.current_logP != -np.inf:
                self.current_direction = c
                return True
        for self.current_direction in range(len(self.space)):
            self.origin = origin[self.current_direction]
            self.move(0.0)
        self.engine.current_logP = self.space.log_posterior(struct)
        self.current_direction = c
        return False

    def set_covariance(self, cov):
        ok = True
        for i, p in enumerate(self.space):
            if cov[i, i] > 0:
                p.width = np.sqrt(cov[i][i])
            else:
                ok = False
        return ok

    def set_covariance_from_hessian(self, struct, h=0.1):
        ok = True
        g = postgetter()
        g.space = self.space
        g.struct = struct
        if self.engine.current_logP is None:
            self.engine.current_logP = self.space.log_posterior(struct)
        chisq1 = -2.0 * self.engine.current_logP
        origin = [p() for p in self.space]
        for i, p in enumerate(self.space):
            trial = origin
            trial[i] = (1.0-h) * origin[i]
            chisq0 = g(trial)
            trial[i] = (1.0 + h) * origin[i]
            chisq2 = g(trial)
            d2 = (chisq2 - 2.0*chisq1 + chisq0) / (h * origin[i])**2
            if d2 > 0.0:
                p.width = 1.0 / np.sqrt(d2)
            else:
                ok = False
        for i, p in enumerate(self.space):
            p.set(origin[i])
        return ok
        # if not have_numdifftools:
        #     print "Error: numdifftools package is required to calculate Hessian matrix"
        #     return False
        # origin = [p() for p in self.space]
        # g = postgetter()
        # g.space = self.space
        # g.struct = struct
        # try:
        #     Hfun = numdifftools.Hessdiag(g)
        #     m = Hfun(origin)
        #     good = True
        # except:
        #     print "CartesianUpdater.set_covariance_from_hessian: warning -- aborting due to Hessian evaluation failure"
        #     good = False
        # for i, p in enumerate(self.space):
        #     p.set(origin[i])
        # if not good:
        #     return False
        # for i, p in enumerate(self.space):
        #     if m[i] > 0.0:
        #         p.width = 1.0 / np.sqrt(m[i])
        #     else:
        #         good = False
        # return good


class SequentialUpdater(Updater):
    """
    Abstract class for updaters that propose each parameter in order.
    """
    def choose_direction(self, j):
        self.current_direction = j

class PermutationUpdater(Updater):
    """
    Abstract class for updaters that propose parameters in random order.
    """
    def __init__(self):
        self.permutation = np.arange(len(self.space))
        if self.count % len(self.space) != 0:
            np.random.shuffle(self.permutation)
    def choose_direction(self, j):
        if j == 0:
            np.random.shuffle(self.permutation)
        self.current_direction = self.permutation[j]

class CartesianSequentialUpdater(CartesianUpdater, SequentialUpdater):
    """
    Updater class to propose parameters individually, in sequence. See Updater.__init__.
    """
    def __init__(self, parameter_space, step, adapt_every=0,
                 adapt_starting=0, on_adapt=None, parallel=None):
        CartesianUpdater.__init__(self, parameter_space, step, adapt_every,
                                  adapt_starting, on_adapt, parallel)

class CartesianPermutationUpdater(CartesianUpdater, PermutationUpdater):
    """
    Updater class to propose parameters individually, in random order. See Updater.__init__.
    """
    def __init__(self, parameter_space, step, adapt_every=0,
                 adapt_starting=0, on_adapt=None, parallel=None):
        CartesianUpdater.__init__(self, parameter_space, step, adapt_every,
                                  adapt_starting, on_adapt, parallel)
        PermutationUpdater.__init__(self)



class MultiDimUpdater(Updater):
    """
    Abstract base class for block updates. Do not instantiate directly.
    """
    def __init__(self, space, step, adapt_every, adapt_starting, on_adapt, parallel):
        Updater.__init__(self, space, step, adapt_every, adapt_starting, on_adapt, parallel)
        if self.adapt:
            self.rescale = 1.0
            self.means = np.zeros(len(self.space))
            self.d = np.zeros(len(self.space))
            self.covariances = np.zeros((len(self.space), len(self.space)))
            if parallel is not None:
                mpi = False
                try:
                    if isinstance(parallel, MPI.Comm):
                        mpi = True
                        self.gatherAdapt = self.gatherMPI
                        self.comm = parallel
                except NameError:
                    pass
                if not mpi:
                    self.pid = str(parallel)
                    self.uind = '_' + str(self.index)
                    self.gatherAdapt = self.gatherFilesys
            else:
                self.gatherAdapt = self.gatherSerial
        self.current_direction = np.zeros(len(self.space))
        self.origin = np.zeros(len(self.space))
        self.widths = [p.width for p in self.space]
        self.width = 0.0
        self.basis = np.eye(len(self.space), len(self.space))
    def __call__(self, struct):
        if self.adapt and self.count >= self.adapt_start and self.count % self.adapt_every == 0:
            self.do_adapt(struct)
        self.choose_direction()
        self.origin = [p() for p in self.space]
        self.step(struct)
        self.accumulate()
    def accumulate(self):
        self.count += 1
        if self.adapt:
            for i, p in enumerate(self.space):
                self.d[i] = p() - self.means[i];
                self.means[i] += self.d[i] / self.count;
                for j in range(i + 1):
                    self.covariances[i,j] += (self.count-1.0) / self.count * self.d[i] * self.d[j]
                    # don't need anything in the upper triangle
    def do_adapt(self, struct):
        cov = self.gatherAdapt()
        if self.set_covariance(cov) and not self.onAdapt is None:
            self.onAdapt(struct)
    def gatherFilesys(self):
        filename = parallel_filename_base + self.pid + self.uind + parallel_filename_ext
        self.save(filename)
        total = 0
        moment1 = np.zeros(len(self.space))
        moment2 = np.zeros((len(self.space), len(self.space)))
        grandMeans = np.zeros(len(self.space))
        grandMeanVar = np.zeros(len(self.space))
        grandVarMean = np.zeros(len(self.space))
        j = 0
        for filename in glob.iglob(parallel_filename_base + '*' + self.uind + parallel_filename_ext):
            try:
                f = open(filename, 'rb')
                s = pickle.load(f)
                f.close()
                total += s['count'] # becomes Ntot
                moment1 += s['count'] * s['means'] # becomes Ntot*<x>
                moment2 += s['count'] * (s['covariances'] / (s['count']-1.0) + \
                                         np.outer(s['means'], s['means'])) # becomes Ntot*<xy>
                d = s['means'] - grandMeans
                grandMeans += d / (j + 1.0)
                grandMeanVar += j / (j + 1.0) * d**2
                d = s['covariances'].diagonal() / (s['count']-1.0) - grandVarMean
                grandVarMean += d / (j + 1.0)
                j += 1
            except:
                print "Warning: IO error while reading " + filename + \
                    " to update covariance (process " + self.pid + ", updater " + self.uind + ")."
        if j > 1:
            B = self.count / (j - 1.0) * grandMeanVar
            W = grandVarMean / j
            self.R = np.sqrt((self.count-1.0) / self.count + B / (self.count*W))
        return (moment2 - np.outer(moment1 / total, moment1)) / (total - 1.0)

    def gatherMPI(self):
        alls = self.saveBits()
        alls = self.comm.allgather(alls)
        total = 0
        moment1 = np.zeros(len(self.space))
        moment2 = np.zeros((len(self.space), len(self.space)))
        grandMeans = np.zeros(len(self.space))
        grandMeanVar = np.zeros(len(self.space))
        grandVarMean = np.zeros(len(self.space))
        for j, s in enumerate(alls):
            total += s['count'] # becomes Ntot
            moment1 += s['count'] * s['means'] # becomes Ntot*<x>
            moment2 += s['count'] * (s['covariances'] / (s['count']-1.0) + \
                                     np.outer(s['means'], s['means'])) # becomes Ntot*<xy>
            d = s['means'] - grandMeans
            grandMeans += d / (j + 1.0)
            grandMeanVar += j / (j + 1.0) * d**2
            d = s['covariances'].diagonal() / (s['count']-1.0) - grandVarMean
            grandVarMean += d / (j + 1.0)
        if len(alls) > 1:
            B = self.count / (len(alls) - 1.0) * grandMeanVar
            W = grandVarMean / len(alls)
            self.R = np.sqrt((self.count-1.0) / self.count + B / (self.count*W))
        return (moment2 - np.outer(moment1 / total, moment1)) / (total - 1.0)

    def gatherSerial(self):
        return self.covariances / (self.count-1.0)

    def move(self, x):
        for i, p in enumerate(self.space):
            p.set(self.origin[i] + x * self.current_direction[i] * self.width)

    def restoreBits(self, s):
        self.count = s['count']
        if s['type'] == 'Cartesian':
            self.means = s['means']
            self.covariances = np.eye(len(self.space), len(self.space)) * s['variances']
            self.widths = s['widths']
            self.basis = np.eye(len(self.space), len(self.space))
        elif  s['type'] == 'MultiDim':
            self.means = s['means']
            self.covariances = s['covariances']
            self.widths = s['widths']
            self.basis = s['basis']
        else:
            raise Exception('MultiDimUpdater.restoreBits: error restoring updater state -- unknown updater type')

    def saveBits(self):
        if self.adapt:
            return {'type': 'MultiDim', 'count': self.count, 'means': self.means,
                    'covariances': self.covariances, 'widths': self.widths, 'basis': self.basis}
        else:
            return None

    def scatter(self, struct, ntries=10):
        for i in range(ntries):
            self.origin = [p() for p in self.space]
            for j in range(len(self.space)):
                self.current_direction = self.basis[:,j]
                self.width = self.widths[j]
                self.move(np.random.randn())
            self.engine.current_logP = self.space.log_posterior(struct)
            if self.engine.current_logP != -np.inf:
                return True
        for j in range(len(self.space)):
            self.current_direction = self.basis[:,j]
            self.move(0.0)
            self.engine.current_logP = self.space.log_posterior(struct)
        return False

    def set_covariance(self, cov):
        try:
            evals, self.basis = np.linalg.eigh(cov)
            self.widths = [self.rescale*np.sqrt(abs(v)) for v in evals]

            if not (np.array(self.widths) != 0.).all():
            
                print "ERROR: (np.array(self.widths) != 0.).all() != 0. Aborting"

                return False
            return False
        except np.linalg.LinAlgError:
            print "MultiDimUpdater.set_covariance: warning -- aborting due to covariance diagonalization failure"
            return False

    def set_covariance_from_hessian(self, struct, h=0.1):
        ok = True
        g = postgetter()
        g.space = self.space
        g.struct = struct
        if self.engine.current_logP is None:
            self.engine.current_logP =  self.space.log_posterior(struct)
        chisq1 = -2.0 * self.engine.current_logP
        self.origin = [p() for p in self.space]
        self.basis = np.eye(len(self.space), len(self.space))
        for i, p in enumerate(self.space):
            trial = self.origin
            trial[i] = (1.0-h) * self.origin[i]
            chisq0 = g(trial)
            trial[i] = (1.0 + h) * self.origin[i]
            chisq2 = g(trial)
            d2 = (chisq2 - 2.0*chisq1 + chisq0) / (h * self.origin[i])**2
            if d2 > 0.0:
                self.widths[i] = 1.0 / np.sqrt(d2)
            else:
                ok = False
        for i, p in enumerate(self.space):
            p.set(self.origin[i])
        return ok
        # if not have_numdifftools:
        #     print "Error: numdifftools package is required to calculate Hessian matrix"
        #     return False
        # self.origin = [p() for p in self.space]
        # g = postgetter()
        # g.verbose = True
        # g.space = self.space
        # g.struct = struct
        # try:
        #     Hfun = numdifftools.Hessian(g, numTerms=0, stepRatio=1.01) #stepNom=self.widths, 
        #     m = Hfun(self.origin)
        #     good = True
        # except:
        #     print "MultiDimUpdater.set_covariance_from_hessian: warning -- aborting due to Hessian evaluation failure"
        #     good = False
        # for i, p in enumerate(self.space):
        #     p.set(self.origin[i])
        # if not good:
        #     return False
        # try:
        #     cov = np.linalg.inv(m)
        #     return self.set_covariances(cov)
        # except np.linalg.LinAlgError:
        #     print "MultiDimUpdater.set_covariance_from_hessian: warning -- aborting due to Hessian inversion failure"
        #     return False


class MDSequentialUpdater(Updater):
    """
    Abstract class for sequential block updates.
    """
    def choose_direction(self):
        j = self.count % len(self.space)
        self.current_direction = self.basis[:,j]
        self.width = self.widths[j]


class MDPermutationUpdater(Updater):
    """
    Abstract class for block updates in random order.
    """
    def __init__(self):
        self.permutation = np.arange(len(self.space))
        if self.count % len(self.space) != 0:
            np.random.shuffle(self.permutation)

    def choose_direction(self):
        j = self.count % len(self.space)
        if j == 0:
            np.random.shuffle(self.permutation)
        self.current_direction = self.basis[:, self.permutation[j]]
        self.width = self.widths[self.permutation[j]]


class MDRotationUpdater(Updater):
    """
    Abstract class for block updates in random directions.
    """
    def __init__(self):
        self.j = 0
        self.q0 = 0
        self.q1 = 1
        self.cosr = 1
        self.sinr = 0

    def choose_direction(self):
        # Choose a random pair of basic vectors and rotate randomly in that plane.
        # On the next call, try the orthogonal direction.
        if self.j == 0:
            self.q0 = np.random.randint(0, len(self.space))
            self.q1 = np.random.randint(0, len(self.space)-1)
            if self.q0 == self.q1:
                self.q1 = len(self.space) - 1
            r = np.random.random() * 2.0 * np.pi
            self.cosr = np.cos(r)
            self.sinr = np.sin(r)
            self.current_direction = self.cosr * self.widths[self.q0] * self.basis[:, self.q0] - \
                                     self.sinr * self.widths[self.q1] * self.basis[:, self.q1]
        else:
            self.current_direction = self.sinr * self.widths[self.q0] * self.basis[:, self.q0] + \
                                     self.cosr * self.widths[self.q1] * self.basis[:, self.q1]
        try:
            self.width = np.sqrt(sum(self.current_direction**2))
            self.current_direction /= self.width

        except np.FloatingPointError as fpe:
            print 'DEBUG'
            print self.widths
            print self.widths[self.q0], self.widths[self.q1]
            print self.sinr, self.cosr
            self.basis
            self.basis[:,self.q0], self.basis[:,self.q1]
            print self.current_direction
            raise fpe

        self.j = (self.j + 1) % 2
    # NB: the below seems more elegant, but appears to have a fatal bug.
#     def choose_direction(self):
#         self.current_direction = np.random.randn(len(self.space))
#         self.current_direction /= np.sqrt(sum(self.current_direction**2))
#         self.width = np.sqrt(1.0 / sum((self.current_direction/self.widths)**2))
#         self.current_direction = np.dot(self.basis, self.current_direction)


class MultiDimSequentialUpdater(MultiDimUpdater, MDSequentialUpdater):
    """
    Updater class to propose block-update parameters in sequence. See Updater.__init__.
    """
    def __init__(self, parameter_space, step, adapt_every=0,
                 adapt_starting=0, on_adapt=None, parallel=None):
        MultiDimUpdater.__init__(self, parameter_space, step, adapt_every,
                                 adapt_starting, on_adapt, parallel)


class MultiDimPermutationUpdater(MultiDimUpdater, MDPermutationUpdater):
    """
    Updater class to block-update parameters in random order. See Updater.__init__.
    """
    def __init__(self, parameter_space, step, adapt_every=0,
                 adapt_starting=0, on_adapt=None, parallel=None):
        MultiDimUpdater.__init__(self, parameter_space, step, adapt_every,
                                 adapt_starting, on_adapt, parallel)
        MDPermutationUpdater.__init__(self)


class MultiDimRotationUpdater(MultiDimUpdater, MDRotationUpdater):
    """
    Updater class to block-update parameters in random directions. See Updater.__init__.
    """
    def __init__(self, parameter_space, step, adapt_every=0,
                 adapt_starting=0, on_adapt=None, parallel=None):
        MultiDimUpdater.__init__(self, parameter_space, step, adapt_every,
                                 adapt_starting, on_adapt, parallel)
        MDRotationUpdater.__init__(self)


try:
    import emcee
    def emcee_lnprobfn(x, up):
        for j, p in enumerate(up.space):
            p.set(x[j])
        return up.space.log_posterior(up.structptr)
            
    class emceeUpdater(Updater):
        """
    Updater that uses the EMCEE Hammer (http://danfm.ca/emcee/ and arxiv:1202.3665).
    Special constructor arguments:
     1* ParameterSpace to update.
     2  Number of walkers, must be even and >= len(ParameterSpace)
     3  Threads/pool arguments as defined by the emcee package for multithread parallelization.
    The initial states are set using the parameter values and .width attributes on initialization.
    Note that MPI parallelization is incompatible with this updater - MPI processes will not communicate.
    Note also that it must be the ONLY updater, containing all parameters exactly once.
    After each call to emceeUpdater, the parameter values are set to those of one of the walkers 
        (a different one each time). The usual backend output should therefore be useful, although 
        arguably it would be better to output the values for all the walkers every N steps. 
        The onStep argument to the Engine could be used for this, or a specialized 
        ParameterSpace/backend combination could be defined.
        """
        def __init__(self, space, nwalkers=None, threads=1, pool=None):
            Updater.__init__(self, space, None, 0, 0, None, None)
            if nwalkers is None:
                nwalkers = 2 * len(space)
            self.sampler = emcee.EnsembleSampler(nwalkers, len(space), emcee_lnprobfn, args=[self],
                                                 threads=threads, pool=pool)
            self.prob = None
            self.rstate = None
            self.pos = emcee.EnsembleSampler.sampleBall([p() for p in space],
                                                        [p.width for p in space], nwalkers)
            self.nwalkers = nwalkers

        def __call__(self, struct):
            self.structptr = struct
            for self.pos, self.prob, self.rstate in self.sampler.sample(self.pos, self.prob,
                                                                        self.rstate, 1):
                pass
            for j, p in enumerate(self.space):
                p.set(self.pos[self.count % self.nwalkers, j])
            self.count += 1

        def restoreBits(self, s):
            if s['type'] == 'Emcee':
                self.sampler = s['sampler']
                self.pos = s['pos']
                self.prob = s['prob']
            else:
                raise Exception('emceeUpdater.restoreBits: incompatible updater type')
        def saveBits(self):
            return {'sampler':self.sampler, 'pos':self.pos, 'prob':self.prob}
except ImportError:
    pass


class Step(object):
    """
    Abstract base class for proposal methods. Do not instantiate directly.
    """
    def __init__(self):
        self.updater = None


class Slice(Step):
    """
    Class implementing the slice proposal method (http://arxiv.org/abs/physics/0009028).
    Constructor arguments:
     1. factor by which to increase the initial slice width.
     2. maximum number of iterations in stepping-out and stepping-in loops before giving up.
     3. whether to suppress warnings if said loops reach the maximum number of iterations.
     4. whether to print a ridiculous amount of information (possibly useful for debugging posterior functions).
    """
    def __init__(self, width_factor=2.4, maxiter=100, quiet=True, obnoxious=False):
        self.width_fac = width_factor
        self.maxiter = maxiter
        self.quiet = quiet
        self.obnoxious = obnoxious
        Step.__init__(self)

    def __call__(self, struct):
        if self.updater.engine.current_logP is None:
            self.updater.engine.current_logP = self.updater.space.log_posterior(struct)
        z = self.updater.engine.current_logP - np.random.exponential() # log level of slice
        L = -self.width_fac * np.random.random_sample()                # left edge of the slice
        R = L + self.width_fac                                         # right edge
        if self.obnoxious:
            print 'Slice: starting params:', [(p.name, p()) for p in self.updater.space]
            print 'Slice: current level', self.updater.engine.current_logP, '; seeking', z
            print 'L & R: ', L, ' ', R
            print 'Slice: stepping out left'
        for i in range(self.maxiter):
            self.updater.move(L)
            lnew = self.updater.space.log_posterior(struct)
            if self.obnoxious:
                print 'Slice: params:', [(p.name, p()) for p in self.updater.space]
                print 'Slice:', L, lnew
            if lnew <= z:
                break
            L -= self.width_fac;
        else:
            if not self.quiet:
                print "Slice(): warning -- exhausted stepping out (left) loop"
        if self.obnoxious:
            print 'Slice: params:', [p() for p in self.updater.space]
            print 'Slice: stepping out right'
        for i in range(self.maxiter):
            self.updater.move(R)
            lnew = self.updater.space.log_posterior(struct)
            if self.obnoxious:
                print 'Slice: params:', [(p.name, p()) for p in self.updater.space]
                print 'Slice:', R, lnew
            if lnew <= z:
                break
            R += self.width_fac;
        else:
            if not self.quiet:
                print "Slice(): warning -- exhausted stepping out (right) loop"
        if self.obnoxious:
            print 'Slice: stepping in'
        for i in range(self.maxiter):
            x1 = L + (R - L) *  np.random.random_sample()
            self.updater.move(x1)
            self.updater.engine.current_logP = self.updater.space.log_posterior(struct)
            if self.obnoxious:
                print 'Slice: params:', [(p.name, p()) for p in self.updater.space]
                print 'Slice:', x1, self.updater.engine.current_logP
            if self.updater.engine.current_logP < z:
                if x1 < 0:
                    L = x1
                else:
                    R = x1
            else:
                break
        else:
            if not self.quiet:
                print "Slice(): warning -- exhausted stepping in loop"
        if self.obnoxious:
            print 'Slice: completed'


class Metropolis(Step):
    """
    Class implementing the Metropolis (*not* Metropolis-Hastings) proposal algorithm.
    Constructor arguments:
     1. a function or functor (with zero arguments) returning a random proposal distance in units
    of the current estimated posterior width. Positive and negative numbers must be returned with
    equal probability. The default is simply a unit Gaussian random number.
    """
    def __init__(self, proposal_length=np.random.randn, width_factor=2.4):
        self.length = proposal_length
        self.width_fac = width_factor
        self.multiplicity = 0
        Step.__init__(self)

    def __call__(self, struct):
        if self.updater.engine.current_logP is None:
            self.updater.engine.current_logP = self.updater.space.log_posterior(struct)
        self.updater.move(self.width_fac * self.length())
        trial_logP = self.updater.space.log_posterior(struct)
        delta_logP = trial_logP - self.updater.engine.current_logP
        r = np.log(np.random.random_sample())
        if delta_logP > 0.0 or r < delta_logP:
            self.updater.engine.current_logP = trial_logP
            self.multiplicity = 1
        else:
            self.updater.move(0.0)
            self.multiplicity += 1


class randNormalExp(object):
    """
    Functor for providing heavy-tailed proposal lengths: Exponential with probability <ratio>
    and Gaussian with probability 1-<ratio>.
    Constructor arguments: ratio.
    """
    def __init__(self, ratio=0.333333333):
        self.ratio = ratio
    def __call__(self):
        r = np.random.randn()
        if r <= self.ratio:
            if r < 0.5*self.ratio:
                return np.random.exponential()
            else:
                return -np.random.exponential()
        else:
            return np.random.randn()


class randChiExp:
    """
    Functor for providing heavy-tailed proposal lengths: Exponential with probability <ratio>
    and Chi(dof)/sqrt(dof) with probability 1-<ratio>. randChiExp(1/3, 2) is the CosmoMC default.
    Constructor arguments:
     1. ratio
     2. degrees of freedom
    """
    def __init__(self, ratio=0.3333333333, dof=2):
        self.ratio = ratio
        self.dof = 2
    def __call__(self):
        r = np.random.randn()
        if r <= self.ratio:
            if r < 0.5*self.ratio:
                return np.random.exponential()
            else:
                return -np.random.exponential()
        else:
            if r < 0.5*(1.0 + self.ratio):
                return np.sqrt(np.random.chisquare(self.dof) / self.dof)
            else:
                return -np.sqrt(np.random.chisquare(self.dof) / self.dof)


# not really necessary to inherit like this, but what the heck
class Backend(object):
    """
    Abstract base class for chain storage. Do not instantaite directly.
    """
    def __call__(self, space):
        pass


class textBackend(Backend):
    """
    Class to store a chain in a text file.
    Constructor argument: an open Python file object.
    Static function readtoDict( ) loads a chain from such a file into a dictionary.
    """
    def __init__(self, file):
        self.file = file
    def __call__(self, space):
        st = ''
        for p in space:
            st = st + ' ' + str(p())
        st = st + '\n'
        self.file.write(st)
    def close(self):
        self.file.close()
    @classmethod
    def readToDict(cls, filename, quiet=True):
        d = None
        f = open(filename, 'r')
        line = f.readline()
        if line != '':
            d = {}
            keys = line.split()
            try:
                values = [float(key) for key in keys]
                for i, val in enumerate(values):
                    keys[i] = 'V' + str(i + 1)
                    d[keys[i]] = [val]
            except ValueError:
                for key in keys:
                    d[key] = []
            while True:
                line = f.readline()
                if line == '':
                    break
                try:
                    values = [float(word) for word in line.split()]
                    for i, key in enumerate(keys):
                        d[key].append(values[i])
                except ValueError:
                    if not quiet:
                        print "textBackend.readToDict: ignoring line " + line
        f.close()
        return d


class headerTextBackend(Backend):
    """
    Like textBackend, but automatically reads/writes a header line with the parameter names.
    """
    def __init__(self, file, space, writeHeader=True):
        self.fields = [p.name for p in space]
        self.writer = csv.DictWriter(file, self.fields, restval='!', delimiter=' ',
                                     quoting=csv.QUOTE_MINIMAL,)
        if writeHeader is True:
            if sys.version_info < (2, 7):
                self.writer.writer.writerow(self.writer.fieldnames)
            else:
                self.writer.writeheader()

    def __call__(self, space):
        towrite = {}
        for p in space:
            towrite[p.name] = p.value
        self.writer.writerow(towrite)

    @classmethod
    def readToDict(cls, filename, quiet=True):
        db = {}
        reader = csv.DictReader(open(filename), delimiter=' ', quoting=csv.QUOTE_MINIMAL)
        for i, row in enumerate(reader):

            for key in row.keys():
                if key not in db:
                    db[key] = []
                try:
                    db[key].append(float(row[key]))
                except TypeError:
                    print i, row
        for key in db.keys():
            db[key] = np.array(db[key])
        return db


class stdoutBackend(textBackend):
    """
    Class to simply print a chain to the terminal without storing it.
    """
    def __init__(self):
        textBackend.__init__(self, sys.stdout)


class dictBackend(dict, Backend):
    """
    Class to store a chain in a dictionary (inherits dict).

    If a Parameter has a non-empty string-type name attribute, the corresponding key
    is that name, otherise it is a reference to the Parameter object itself.
    """
    def __call__(self, space):
        for p in space:
            key = p
            try:
                if p.name != '':
                    key = p.name
            except:
                pass
            try:
                self[key].append(p())
            except KeyError:
                self[key] = []
                self[key].append(p())


class Engine(list):
    """
    Class to organize Updaters of ParameterSpaces and run the MCMC (inherits list).
    Constructor arguments:
     1. sequence of Updater objects. If Updaters are added any other way, the register_updater( )
    method must be used.
     2. a ParameterSpace of Parameters whose values are to be stored at each step. This need not
    be the same as the ParameterSpace(s) referred to by the Updaters.
     3. a function of one argument to be called after each step (i.e. each time that each Updater
    has been called).
    To run a chain, use the () method. Arguments:
     1. number of iterations (every Updater is called for a single iteration).
     2. an object that is passed to the log_posterior, Updater.on_adapt, and on_step functions.
     3. a sequence of Backend objects where the chain is to be stored.
    """
    # todo: make sure directly assigned Updaters get registered
    def __init__(self, updaterList=(), parameterspace_to_track=None, on_step=None):
        list.__init__(self, updaterList)
        for i, updater in enumerate(self):
            self.register_updater(updater, i)
        self.space = parameterspace_to_track
        self.onStep = on_step
        self.count = 0
        self.current_logP = None

    def __setitem__(self, key, value):
        self[key] = value
        self.register_updater(value, key)

    def __call__(self, number=1, struct=None, backends=(stdoutBackend())):
        try:
            for i in range(number):
                if i % 200 == 0:
                    print 'At Iteration %d' % i
                for updater in self:
                    for j in range(updater.rate):
                        updater(struct)
                self.count += 1
                if not self.onStep is None:
                    self.onStep(struct)
                if not self.space is None:
                    for backend in backends:
                        backend(self.space)
        except KeyboardInterrupt:
            print "Interrupted by keyboard with count = " + str(self.count)

    def register_updater(self, updater, index):
        updater.engine = self
        updater.index = index
        updater.uind = '_' + str(index)


def example(number=None):
    if number == 1:
        print """
# Here is a simple example. As shown it will run in non-parallel mode; comments indicate what 
# to do for parallelization.

from mymc import *
## for MPI
#from mpi4py import MPI
#mpi_rank = MPI.COMM_WORLD.Get_rank()

### Define some parameters.
x = Parameter(name='x')
y = Parameter(name='y')

### This is the object that will be passed to the likelihood function.
### In this simple case, it just holds the parameter objects, but in general it could be anything.
### E.g., usually it would also contain or point to the data being used to constrain the model. 
### A good idea is to write the state of any updaters to a file after each adaptation (using the 
### on_adapt functionality), in which case keeping pointers to the updaters here is convenient. 
### Also commonly useful: a DerivedParameter which holds the value of the posterior log-density 
### for each sample.
class Thing:
    def __init__(self, x, y):
        self.x = x
        self.y = y
thing = Thing(x, y)

### The log-posterior function. Here we just assume a bivariate Gaussian posterior with marginal 
### standard deviations s(x)=2 and s(y)=3, correlation coefficient 0.75, and means <x>=-1, <y>=1.
def post(thing):
    r = 0.75
    sx = 2.0
    sy = 3.0
    mx = -1.0
    my = 1.0
    return -0.5/(1.0-r**2)*((thing.x()-mx)**2/sx**2 + (thing.y()-my)**2/sy**2 - 2.0*r*(thing.x()-mx)/sx*(thing.y()-my)/sy)

### Create a parameter space consisting of x and y, and associate the log-posterior function with it.
space = ParameterSpace([thing.x, thing.y], post)

### If we'd bothered to define a DerivedParameter in Thing which would hold the posterior density, 
### we might want to define a larger ParameterSpace and pass it to the Engine later on to be saved 
### in the Backends (instead of space).
#trace = ParameterSpace([thing.x, thing.y, thing.logP])

### Use slice sampling for robustness. Adapt the proposal distribution every 100 iterations 
### starting with the 100th.
step = Slice()
parallel = None
## for MPI parallelization
# parallel = MPI.COMM_WORLD
## for parallelization via the filesystem, this would have to be set to a different value for 
## each concurrently running instance
#parallel = 1
updater = MultiDimSequentialUpdater(space, step, 100, 100, parallel=parallel)

### Create an Engine and tell it to drive this Updater and to store the values of the free parameters.
engine = Engine([updater], space)

### Store the chain in a text file.
chainfile = open("chain.txt", 'w')
## For filesystem parallelization, each instance should write to a different file.
## For MPI, the same is true, e.g.
#chainfile = open("chain" + str(MPI.COMM_WORLD.Get_rank()) + ".txt", 'w')
backends = [ textBackend(chainfile) ]

### Print the chain to the terminal as well
backends.append(stdoutBackend())

### Run the chain for 10000 iterations
engine(10000, thing, backends)

### Close the text file to clean up.
chainfile.close()

## If this was a parallel run, print the convergence criterion for each parameter.
# print updater.R
"""

    else:
        print """
Usage: example(N), where N is one of:
 1. A very simple example where the posterior is bivariate Gaussian, to illustrate setting 
    up and running the engine.
"""


class ChiSquareLikelihood(object):
    """
    A class to simplify fitting models to Gaussian data.

    Assign a function to the 'priors' attribute to include non-(improper uniform) priors.
    """
    def __init__(self, model, y, err=None, x=None):
        self.model = model
        self.y = y
        if err is None:
            self.err = np.ones(len(y))
        else:
            self.err = err
        if x is None:
            self.x = np.zeros(len(y))
        else:
            self.x = x
        self.chisquare = DerivedParameter(name='ChiSquare')

    def __call__(self, struct):
        try:
            chisq = -2.0 * self.priors(struct)
        except AttributeError:
            chisq = 0.0
        for j, y in enumerate(self.y):
            chisq += ((self.model(self.x[j], struct) - y) / self.err[j])**2
        self.chisquare.value = chisq
        return -0.5 * chisq

# Todo:
# 1. An Updater class that simply goes through an existing sequence, for importance sampling.

