"""Convert a pymc model to something runnable by mymc."""


import csv
import sys
import os
import cPickle
import operator
import pymc
import numpy as np
try:
    from mpi4py import MPI
except ImportError:
    pass
from . import mymc
from . import util


class CompositeParameter(mymc.Parameter):
    def __init__(self, masterobj, index, width=0.1):
        self.masterobj = masterobj
        self.index = index
        self.name = '%s_%d' % (self.masterobj.__name__, self.index)

        self.width = width*np.abs(self())

        print self.name, self.width

    def get_value(self):
        return self.masterobj.value[self.index]

    def set(self, value):
        curval = np.copy(self.masterobj.value)   #to not disturb pymc caching mechanism
        curval[self.index] = value
        self.masterobj.value = curval

    value = property(get_value, set)


class WrapperParameter(mymc.Parameter):
    def __init__(self, masterobj, width=0.1):
        self.masterobj = masterobj
        self.name = self.masterobj.__name__

        self.width = width*np.abs(self())

        print self.name, self.width

    def get_value(self):
        return self.masterobj.value

    def set(self, value):
        self.masterobj.value = value

    value = property(get_value, set)


#################################


class DerivedWrappedParameter(mymc.DerivedParameter):
    def __init__(self, masterobj):
        self.masterobj = masterobj
        self.name = self.masterobj.__name__
    def get_value(self):
        return self.masterobj.value
    value = property(get_value)


##################################


class DerivedAttribute(mymc.DerivedParameter):
    def __init__(self, masterobj, attr):
        self.masterobj = masterobj
        self.attr = attr
        self.name = self.attr

    def get_value(self):
        return getattr(self.masterobj, self.attr)

    value = property(get_value)

#################################


class DerivedFunction(mymc.DerivedParameter):
    def __init__(self, func, name, *args, **kw):
        self.func = func
        self.args = args
        self.kw = kw
        self.name = name

    def get_value(self):
        return self.func(*self.args, **self.kw)

    value = property(get_value)

#################################


def wrapModel(model):

    parameters = []
    deterministics = []

    stochastics = []
    potentials = []
    observed = []

    for s in model.stochastics:
        if s.observed:
            observed.append(s)
            continue

        stochastics.append(s)

        if s.value.shape == ():

            parameters.append(WrapperParameter(s))

        else:
            for i in np.arange(len(s.value)):
                parameters.append(CompositeParameter(s, i))

    for d in model.deterministics:
        if d.keep_trace:
            deterministics.append(DerivedWrappedParameter(d))

    for p in model.potentials:
        potentials.append(p)

    current_names = [x.__name__ for x in stochastics]
    for o in model.observed_stochastics:
        if o.__name__ not in current_names:
            observed.append(o)


    parameters = sorted(parameters, key=operator.attrgetter('name'))
    deterministics = sorted(deterministics, key=operator.attrgetter('name'))

    stochastics = sorted(stochastics, key=operator.attrgetter('__name__'))
    potentials = sorted(potentials, key=operator.attrgetter('__name__'))
    observed = sorted(observed, key=operator.attrgetter('__name__'))

    print [x.__name__ for x in stochastics]
    print [x.__name__ for x in potentials]
    print [x.__name__ for x in observed]

    all_logp = stochastics + potentials + observed

    print [x.__name__ for x in all_logp]

    def posterior(thing):
        try:
            logp = reduce(lambda x, y: x + y.logp, all_logp, 0.)
        except pymc.ZeroProbability as zpexc:
            logp = -np.infty
        return logp

    def likelihood():
        try:
            logp = reduce(lambda x, y: x + y.logp, observed, 0.)
        except pymc.ZeroProbability as zpexc:
            logp = -np.infty
        return logp

    deterministics.append(DerivedFunction(likelihood, 'likelihood'))
    deterministics.append(DerivedFunction(posterior, 'posterior', None))

    space = mymc.ParameterSpace(parameters, posterior)

    trace = mymc.ParameterSpace(deterministics + parameters)

    return space, trace

#################################


class MyMCRunner(object):

    def run(self, manager):

        options = manager.options

        parallel = None
        if not options.singlecore:
            parallel = MPI.COMM_WORLD
            manager.mpi_rank = MPI.COMM_WORLD.Get_rank()
            print 'Rank: ', manager.mpi_rank
            print 'World Size: ', MPI.COMM_WORLD.Get_size()
        else:
            parallel = None
            manager.mpi_rank = 0


        space, trace = wrapModel(manager.model)
        step = mymc.Slice()
        updater = mymc.MultiDimRotationUpdater(space, step, options.adapt_every,
                                               options.adapt_after, parallel=parallel)

        bitsfile = '%s.bits.%d' % (options.outputFile, manager.mpi_rank)
        chainfile = '%s.chain.%d' % (options.outputFile, manager.mpi_rank)

        writeHeader = True
        if options.restore is True:

            #  load previous proposal distribution
            if os.path.exists(bitsfile):
                with open(bitsfile, 'rb') as input:
                    updater.restoreBits(cPickle.load(input))

            ## initialize chain to last sampled value
            if os.path.exists(chainfile):

                writeHeader = False
                chain = util.loadchains([chainfile])
                for param in space:
                    param.set(chain[param.name][0, -1])

        manager.engine = mymc.Engine([updater], trace)

        if os.path.exists(chainfile) and writeHeader is True:
            os.remove(chainfile)

        manager.chainfile = open(chainfile, 'a')
        manager.textout = mymc.headerTextBackend(manager.chainfile, trace, writeHeader=writeHeader)

        backends = [manager.textout]

        manager.engine(options.nsamples, None, backends)

        with open(bitsfile, 'wb') as output:
            cPickle.dump(updater.saveBits(), output)

    def addCLOps(self, parser):
        """Add options."""
        parser.add_option('-s', '--nsamples', dest='nsamples',
                          help='Number of MCMC samples to draw or scan model',
                          default=None, type='int')
        parser.add_option('--adaptevery', dest='adapt_every',
                          help='Adapt MCMC chain every X steps', default=100, type='int')
        parser.add_option('--adaptafter', dest='adapt_after',
                          help='Start adapting MCMC chain aftter X steps', default=100, type='int')
        parser.add_option('--burn', dest='burn',
                          help='Number of MCMC samples to discard before calculated mass statistics',
                          default=10000, type=int)
        parser.add_option('--singlecore', default=False,
                          action='store_true',
                          help='Turn off MPI for test runs on single machines')

    def dump(self, manager):
        """Dump masses."""
        dumpMasses(np.array(manager.chain['mass_15mpc'][manager.options.burn:]),
                   '%s.mass15mpc.%d' % (manager.options.outputFile, manager.mpi_rank))

    def finalize(self, manager):
        """Close all."""
        manager.chainfile.close()


class MyMCMemRunner(object):

    def run(self, manager):

        options = manager.options

        parallel = None
        if not options.singlecore:
            parallel = MPI.COMM_WORLD
            manager.mpi_rank = MPI.COMM_WORLD.Get_rank()
        else:
            parallel = None
            manager.mpi_rank = 0

        space, trace = wrapModel(manager.model)

        step = mymc.Slice()

        if len(space) == 1:
            updater = mymc.CartesianSequentialUpdater(space, step, options.adapt_every,
                                                      options.adapt_after, parallel=parallel)
        else:
            updater = mymc.MultiDimRotationUpdater(space, step, options.adapt_every,
                                                   options.adapt_after, parallel=parallel)

        manager.engine = mymc.Engine([updater], trace)
        manager.chain = mymc.dictBackend()
        backends = [manager.chain]

        manager.engine(options.nsamples, None, backends)

    def addCLOps(self, parser):
        pass

    def dump(self, manager):

        options = manager.options
        if not 'outputFile' in options:
            return

        with open(options.outputFile, 'w') as output:

            chain = manager.chain

            fields = chain.keys()
            nrows = len(chain[fields[0]])

            writer = csv.DictWriter(output, fields, restval='!', delimiter=' ',
                                    quoting=csv.QUOTE_MINIMAL, )

            if sys.version_info < (2, 7):
                writer.writer.writerow(self.writer.fieldnames)
            else:
                writer.writeheader()

            for i in range(nrows):
                towrite = {}
                for field in fields:
                    towrite[field] = chain[field][i]
                writer.writerow(towrite)

    def finalize(self, manager):
        pass


def dumpMasses(masses, outputFile):

    with open('%s.mass.pkl' % outputFile, 'wb') as output:
        cPickle.dump(masses, output)

    mean = np.mean(masses)
    stddev = np.std(masses)
    quantiles = pymc.utils.quantiles(masses, qlist=[2.5, 15.8, 25, 50, 75, 84.1, 97.5])
    hpd68 = pymc.utils.hpd(masses, 0.32)
    hpd95 = pymc.utils.hpd(masses, 0.05)
    ml, (m, p) = ci.maxDensityConfidenceRegion(masses)
    lml, (lm, lp) = ci.maxDensityConfidenceRegion(np.log10(masses))

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
