#############################
# Handles loading files for a simulation run
#############################

from __future__ import with_statement
import ldac, cPickle, os, subprocess, copy, sys
import numpy as np
import astropy.io.fits as pyfits
import shearprofile as sp, bashreader
import maxlike_subaru_filehandler as msf, maxlike_general_filehandler, varcontainer, utilities

#############################

class SubaruSecureFilehandler(msf.SubaruFilehandler):

    def addCLOps(self, parser):

        super(SubaruSecureFilehandler, self).addCLOps(parser)

        parser.add_option('--workdir', dest='workdir',
                          help='location to read/write all files')
        parser.add_option('-i', '--incatalog', dest='incatalog',
                          help='Input catalog')



    #############################

    def createOptions(self, workdir, incatalog, *args, **keywords):

        options, args = super(SubaruSecureFilehandler, self).createOptions(*args, **keywords)

        options.workdir = workdir
        options.incatalog = incatalog


        return options, args

    #############################

    def buildFileNames(self, manager, newoptions):

        options = manager.options


        newoptions.lensingcat = options.incatalog

        newoptions.bpzfile = '%s/%s.%s.bpz.tab' % (options.workdir, options.cluster,
                                                   options.filter)

        newoptions.inputPDZ = '%s/%s.%s.pdz.cat' % (options.workdir, options.cluster,
                                                    options.filter)

        newoptions.neighborcat = '%s/%s.%s.%s.neighbors.cat' % (options.workdir, options.cluster, 
                                                                options.filter, options.image)

        

    #############################
