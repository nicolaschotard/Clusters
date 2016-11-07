#############################
# Handles loading files for a simulation run
#############################

import ldac, cPickle, numpy as np
import astropy.io.fits as pyfits
import shearprofile as sp
import maxlike_subaru_filehandler as msf
import readtxtfile

#############################

__cvs_id__ = "$Id$"

#############################

def loadBetaSamples(options):

    inputfile = '%s/%s.%s.%s.betas.cat' % (options.workdir, 
                                             options.cluster, 
                                             options.filter, 
                                             options.image)
    betasamples = readtxtfile.readtxtfile(inputfile)

    return betasamples

#############################

def loadContamFrac(options):

    inputfile = '%s/contamfrac.dat' % options.workdir

    contamfracmean = None
    contamfracsigma = None

    rawdat = readtxtfile.readtxtfile(inputfile)
    for line in rawdat:
        if line[0] == 'average:':
            contamfracmean = float(line[1])
        if line[1] == 'sigma:':
            contamfracsigma = float(line[1])

    if contamfracmean is None or contamfracsigma is None:
        raise msf.WeirdnessException

    return contamfracmean, contamfracsig

#############################

def loadNShells(options):

    inputfile = '%s/%s.%s.%s.nshells' % (options.workdir, 
                                             options.cluster, 
                                             options.filter, 
                                             options.image)

    rawdat = readtxtfile(inputfile)

    return rawdat[0][0]

#############################

class ColorcutFilehandler(object):

    def addCLOps(self, parser):

        parser.add_option('--workdir', dest='workdir',
                          help='location to read/write all files')
        parser.add_option('-i', '--incatalog', dest='incatalog',
                          help='Input catalog')

        parser.add_option('-c', '--cluster', dest='cluster',
                          help='Cluster name')
        parser.add_option('-f', '--filter', dest='filter',
                          help='Detection and lensing filter')
        parser.add_option('-v', '--image', dest='image',
                          help='Lensing image, (eg. good, gabodsid3088)')
        parser.add_option('--ldaclensing', dest='ldaclensing',
                          help='location of the ldaclensing directory', default='../ldaclensing')


#############################
        
    def createOptions(self, workdir, incatalog, cluster, filter, image, ldaclensing='../ldaclensing', 
                      *args, **keywords):

        options.workdir = workdir
        options.incatalog = incatalog
        options.cluster = cluster
        options.filter = filter
        options.image = image
        options.ldaclensing = ldaclensing



        return options, args

#########################

    def buildFileNames(self, manager):

        options = manager.options


        options.lensingcat = options.incatalog

        options.neighborcat = '%s/%s.%s.%s.neighbors.cat' % (options.workdir, options.cluster, 
                                                                options.filter, options.image)

        

    #############################


    def readData(self, manager):

        options = manager.options


        self.buildFileNames(manager, options)


        manager.open('lensingcat', options.lensingcat, ldac.openObjectFile)

        if options.nearestneighbor is True:
            nncat = ldac.openObjectFile(options.neighborcat)
            manager.nearestneighbors = nncat.matchById(manager.lensingcat)

            


        options.psfsize = msf.readPSFSize(options.workdir, options.cluster, options.filter, options.image)

        options.zcluster = msf.parseZCluster(options.cluster)
        options.r500 = msf.readR500(options.cluster)
        options.centerx, newoptions.centery = msf.readClusterCenters(options.cluster)







        options.pixscale = 0.2

        options.xcol = 'Xpos'
        options.ycol = 'Ypos'

        options.g1col = 'gs1'
        options.g2col = 'gs2'

        options.sizecol = 'rh'

        options.centerx = 5000
        options.centery = 5000


        options.snratio = 'snratio_scaled1'



        manager.zcluster = options.zcluster
        manager.r500 = options.r500
        manager.psfsize = options.psfsize
        manager.pixscale = options.pixscale


        if 'lensingcat' not in manager:
            manager.open('lensingcat', options.lensingcat, ldac.openObjectFile)

        r_arc, E, B = sp.calcTangentialShear(cat = manager.lensingcat, 
                                             center = (options.centerx, options.centery),
                                             pixscale = options.pixscale,
                                             xcol = options.xcol,
                                             ycol = options.ycol,
                                             g1col = options.g1col,
                                             g2col = options.g2col)

        r_pix = r_arc / options.pixscale

        r_mpc = r_arc * (1./3600.) * (np.pi / 180. ) * sp.angulardist(options.zcluster)

        size = manager.lensingcat[options.sizecol] / options.psfsize
        snratio = manager.lensingcat[options.snratio]


        cols = [pyfits.Column(name='SeqNr', format = 'J', array = manager.lensingcat['SeqNr']),
                pyfits.Column(name = 'r_mpc', format = 'E', array = r_mpc),
                pyfits.Column(name = 'size', format = 'E', array = size),
                pyfits.Column(name = 'snratio', format = 'E', array = snratio),
                pyfits.Column(name = 'ghats', format = 'E', array = E),
                pyfits.Column(name = 'B', format = 'E',  array = B)]

        manager.store('inputcat', ldac.LDACCat(pyfits.BinTableHDU.from_columns(pyfits.ColDefs(cols))))


        # Need to load beta, beta2 distributions

        beta_samples = loadBetaSamples(options)
        manager.store('beta_samples', beta_samples)        

        # Need to load contam frac info

        manager.store('contamfrac', loadContamFrac(options))

        # need to load nshells for each cluster

        manager.store('nshells', loadNShells(options))

        ####


                 
