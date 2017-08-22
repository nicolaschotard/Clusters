"""Drives Maxlike Mass calculations."""

from __future__ import with_statement
import optparse
import inspect
import os
import time
from . import datamanager
from . import util


usage = "To be defined later"


class Controller(datamanager.DataManager):

    def __init__(self, modelbuilder, filehandler, runmethod):

        super(Controller, self).__init__()

        self.modelbuilder = modelbuilder
        self.filehandler = filehandler
        self.runmethod = runmethod

        self.dump_strategies()

        self.options = None
        self.args = None

    def dump_strategies(self):

        def comment_string(strategy):

            if inspect.isfunction(strategy) or inspect.ismethod(strategy):
                srcfile = inspect.getfile(strategy)
                moddate = time.asctime(
                    time.localtime(os.stat(srcfile).st_mtime))
                return '%s defined in %s (Mod: %s)' % (strategy.__name__, srcfile, moddate)
            elif inspect.ismodule(strategy):
                srcfile = inspect.getfile(strategy)
                moddate = time.asctime(
                    time.localtime(os.stat(srcfile).st_mtime))
                return 'Defined in %s (Mod: %s)' % (srcfile, moddate)
            else:
                srcfile = inspect.getfile(strategy.__class__)
                moddate = time.asctime(
                    time.localtime(os.stat(srcfile).st_mtime))
                return '%s defined in %s (Mod: %s)' % \
                    (strategy.__class__.__name__, srcfile, moddate)

        comment = '''
Configured With:
    ModelBuilder : %(modelbuilder)s
    FileHandler  : %(filehandler)s
    RunMethod    : %(runmethod)s
''' % {'modelbuilder': comment_string(self.modelbuilder),
            'filehandler': comment_string(self.filehandler),
            'runmethod': comment_string(self.runmethod)}
        self.comment(comment)

    def run_all(self):

        self.parse_cl()
        self.load()
        self.run()
        self.dump()
        self.finalize()

    def addCLOps(self, parser):

        parser.add_option('-o', '--outfile', dest='outputFile',
                          help='Basename for output results', metavar='FILE')

    def parse_cl(self):

        parser = optparse.OptionParser(usage=usage)

        def add_group(strategy, printedname):

            if hasattr(strategy, 'addCLOps'):
                optgroup = optparse.OptionGroup(
                    parser, '%s Options' % printedname, '')
                strategy.addCLOps(optgroup)
                parser.add_option_group(optgroup)

        add_group(self, 'Controller')
        add_group(self.modelbuilder, 'Model')
        add_group(self.filehandler, 'File Handler')
        add_group(self.runmethod, 'Run Method')

        varoptions = util.VarContainer()
        options, args = parser.parse_args()
        for key, val in vars(options).iteritems():
            varoptions[key] = val

        self.replace('options', varoptions)
        self.replace('args', args)

        return options, args

    def load(self, options=None, args=None):

        if options is not None:
            self.replace('options', options)
        if args is not None:
            self.replace('args', args)

        self.comment('load called with options = %s\n\tand args = %s' %
                     (str(self.options), str(self.args)))

        if hasattr(self.modelbuilder, 'load'):
            self.modelbuilder.load(self)

        self.filehandler.readData(self)

        if hasattr(self.filehandler, 'cuts'):
            for cut in self.filehandler.cuts:
                self.update(cut, {'inputcat': 'filter', 'pz': '__getitem__'})

        if hasattr(self.modelbuilder, 'cuts'):
            for cut in self.modelbuilder.cuts:
                self.update(cut, {'inputcat': 'filter', 'pz': '__getitem__'})

        self.model = self.modelbuilder.createModel(self)

    def run(self):

        self.runmethod.run(self)

    def dump(self):

        self.runmethod.dump(self)

        outputFile = self.options.outputFile

        self.history.dump('%s.log' % outputFile)

    def finalize(self):

        self._do_if_there(self.filehandler, 'finalize')
        self._do_if_there(self.modelbuilder, 'finalize')
        self._do_if_there(self.runmethod, 'finalize')

    def _do_if_there(self, strategy, funcname):

        try:
            getattr(strategy, funcname)(self)
        except AttributeError:
            pass
