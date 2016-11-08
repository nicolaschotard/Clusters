###################################
# File handler for an astropy table file format
# Assumes that all relevant cuts have been applied to the file already
####################################

import numpy as np
import astropy.table as table
import astropy.io.fits as pyfits
import varcontainer
import nfwutils, clusterTools
import ldac

####################################



####################################

class AstropyTableFilehandler(object):

    ######

    def __init__(self):

        self.cuts =  []

    ######

    def addCLOps(self, parser):
        #not designed to be run from the command line

        raise NotImplementedError
    
    ######

    def createOptions(self, cluster, zcluster,
                      lensingcat,
                      pdzfile,
                      cluster_ra,
                      cluster_dec,
                      raCol = 'coord_ra_deg',
                      decCol = 'coord_dec_deg',
                      g1Col = 'ext_shapeHSM_HsmShapeRegauss_e1',
                      g2Col = 'ext_shapeHSM_HsmShapeRegauss_e2',
                      options = None, args = None):

        if options is None:
            options = varcontainer.VarContainer()

        options.cluster = cluster
        options.zcluster = zcluster
        options.lensingcat = lensingcat
        options.cluster_ra = cluster_ra
        options.cluster_dec = cluster_dec
        options.pdzfile = pdzfile


        options.raCol      = raCol
        options.decCol     = decCol
        options.g1Col      = g1Col
        options.g2Col      = g2Col


        return options, args

    ######

    def readData(self, manager):

        options = manager.options

        manager.lensingcat = options.lensingcat
        
        manager.clustername = options.cluster
        manager.zcluster    = options.zcluster



        r_arcmin, E, B = calcTangentialShear(cat = manager.lensingcat, 
                                             center = (options.cluster_ra, options.cluster_dec),
                                             raCol = options.raCol,
                                             decCol = options.decCol,
                                             g1Col = options.g1Col,
                                             g2Col = options.g2Col)


        r_mpc = r_arcmin * (1./60.) * (np.pi / 180. ) * nfwutils.global_cosmology.angulardist(options.zcluster)

#        size = manager.lensingcat[options.sizeCol] / options.psfsize
#        snratio = manager.lensingcat[options.snratioCol]

        manager.open('pdzcat', options.pdzfile, table.Table.read, path='pdz_values')
        manager.open('pdzrange', options.pdzfile, table.Table.read, path='pdz_bins')

        manager.matched_pdzcat = matchById(manager.pdzcat, manager.lensingcat, 'id', 'objectId')
        print manager.matched_pdzcat.keys()
        manager.pz = manager.matched_pdzcat['pdz']  #area normalized, ie density function

        z_b = manager.matched_pdzcat['Z_BEST']

        cols = [pyfits.Column(name='SeqNr', format = 'J', array = manager.lensingcat['id']),
                pyfits.Column(name = 'r_mpc', format = 'E', array = r_mpc),
#                pyfits.Column(name = 'size', format = 'E', array = size),
#                pyfits.Column(name = 'snratio', format = 'E', array = snratio),
                pyfits.Column(name = 'z_b', format = 'E', array = z_b),
                pyfits.Column(name = 'ghats', format = 'E', array = E),
                pyfits.Column(name = 'B', format = 'E',  array = B)]

        manager.store('inputcat', ldac.LDACCat(pyfits.BinTableHDU.from_columns(pyfits.ColDefs(cols))))





#############



def calcTangentialShear(cat, center, raCol, decCol, g1Col, g2Col):

    cluster_ra, cluster_dec = center
    ra = cat[raCol]
    dec = cat[decCol]
    e1 = cat[g1Col]
    e2 = cat[g2Col]

    posangle = clusterTools.positionAngle(ra, dec, cluster_ra, cluster_dec) #radians
    r_arcmin = clusterTools.greatCircleDistance(ra, dec, cluster_ra, cluster_dec)*60

    cos2phi = np.cos(2*posangle)
    sin2phi = np.sin(2*posangle)

    E = -(e1*cos2phi+e2*sin2phi)

    b1 =  e2
    b2 = -e1
    B = -(b1*cos2phi+b2*sin2phi)

    return r_arcmin, E, B


##########


def matchById(firstcat, othercat, otherid='SeqNr', selfid='SeqNr'):
    '''Returns a subset of this catalog, that matches the order of the provided catalog'''



    order = {}
    for i, x in enumerate(firstcat[selfid]):
        order[x] = i

    keepOrder = []
    for x in othercat[otherid]:
        if x in order:
            keepOrder.append(order[x])

    keep = np.array(keepOrder)
    matched = firstcat[keep]
    return matched
