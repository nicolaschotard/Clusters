###################################
# File handler for an astropy table file format
# Assumes that all relevant cuts have been applied to the file already
####################################

import numpy as np
import astropy.table as table
import astropy.io.fits as pyfits
from . import varcontainer
from . import nfwutils
from . import sphereGeometry
from . import ldac

####################################


class AstropyTableFilehandler(object):

    ######

    def __init__(self):

        self.cuts = []

    ######

    def addCLOps(self, parser):
        #not designed to be run from the command line

        raise NotImplementedError

    ######

    def createOptions(self, cluster, zcluster,
                      cat, mconfig,
                      cluster_ra,
                      cluster_dec,
                      raCol='coord_ra_deg',
                      decCol='coord_dec_deg',
                      g1Col='ext_shapeHSM_HsmShapeRegauss_e1',
                      g2Col='ext_shapeHSM_HsmShapeRegauss_e2',
                      sizeCol='rh',
                      snratioCol='snratio_scaled1',
                      wtg_shearcal=False,
                      psfsize=None,
                      logprior=False,
                      options=None, args=None):

        
        if options is None:
            options = varcontainer.VarContainer()

        options.cluster = cluster
        options.zcluster = zcluster
        options.cluster_ra = cluster_ra
        options.cluster_dec = cluster_dec
        options.wtg_shearcal = wtg_shearcal
        options.logprior = logprior
        
        if wtg_shearcal:
            options.psfsize = psfsize
            options.sizeCol = sizeCol
            options.snratioCol = snratioCol
        
        options.cat = cat
        options.mconfig=mconfig

        options.raCol = raCol
        options.decCol = decCol
        options.g1Col = g1Col
        options.g2Col = g2Col

        return options, args

    ######

    def readData(self, manager):

        options = manager.options

        manager.lensingcat = options.cat['deepCoadd_meas']
        manager.zcat = options.cat[options.mconfig['zconfig']]

        manager.clustername = options.cluster
        manager.zcluster = options.zcluster

        manager.logprior = options.logprior
        manager.wtg_shearcal = options.wtg_shearcal
        
        r_arcmin, E, B, phi = calcTangentialShear(cat=manager.lensingcat,
                                             center=(options.cluster_ra, options.cluster_dec),
                                             raCol=options.raCol,
                                             decCol=options.decCol,
                                             g1Col=options.g1Col,
                                             g2Col=options.g2Col)

        
        r_mpc = r_arcmin * (1. / 60.) * (np.pi / 180.) * \
                nfwutils.global_cosmology.angulardist(options.zcluster)

        if options.wtg_shearcal:
            manager.psfsize = options.psfsize
            size = manager.lensingcat[options.sizeCol] / options.psfsize
            snratio = manager.lensingcat[options.snratioCol]

##       old version                
#        manager.open('pdzcat', options.pdzfile, table.Table.read, path=options.prefix + 'pdz_values')
#        manager.open('pdzrange', options.pdzfile, table.Table.read, path=options.prefix + 'pdz_bins')
#        manager.replace('pdzrange', lambda: manager.pdzrange['zbins'])

        manager.pdzrange = manager.zcat['zbins'][0]  # all objects have same zbins, take the first one
        manager.replace('pdzrange', lambda: manager.pdzrange) 
        
        # only keep 'i' filter
        if 'filter' in manager.lensingcat.keys():
            manager.replace('lensingcat', manager.lensingcat[manager.lensingcat["filter"] == 'i'])

        # redshift cut
#       if 'z_flag_pdz_' + options.prefix[:-1] in manager.lensingcat.keys():
        if 'flag_' + options.mconfig['zconfig'] in options.cat.keys():
            if 'zflagconfig' in options.mconfig and options.mconfig['zflagconfig'] == 'pdz':
                print "Using pdz flag", len(manager.lensingcat[options.cat["flag_" + options.mconfig['zconfig']]['flag_z_pdz'] == True])
                manager.replace('lensingcat',
                                manager.lensingcat[options.cat["flag_" + options.mconfig['zconfig']]['flag_z_pdz'] == True])
            elif 'zflagconfig' in options.mconfig and options.mconfig['zflagconfig'] == 'hard':
                print "Using hard flag",len(manager.lensingcat[options.cat["flag_" + options.mconfig['zconfig']]['flag_z_hard'] == True])
                manager.replace('lensingcat',
                                 manager.lensingcat[options.cat["flag_" + options.mconfig['zconfig']]['flag_z_hard'] == True])

        manager.matched_zcat = matchById(manager.zcat, manager.lensingcat, 'id', 'objectId')
        manager.pz = manager.matched_zcat['pdz']  # area normalized, ie density function

        z_b = manager.matched_zcat['Z_BEST']

        if manager.wtg_shearcal:
            cols = [pyfits.Column(name='SeqNr', format='J', array=manager.lensingcat['id']),
                    pyfits.Column(name='r_mpc', format='E', array=r_mpc),
                    pyfits.Column(name='size', format='E', array=size),
                    pyfits.Column(name='snratio', format='E', array=snratio),
                    pyfits.Column(name='z_b', format='E', array=z_b),
                    pyfits.Column(name='ghats', format='E', array=E),
                    pyfits.Column(name='B', format='E', array=B),
                    pyfits.Column(name='phi', format='E', array=phi)]
        else:
            cols = [pyfits.Column(name='SeqNr', format='J', array=manager.lensingcat['id']),
                    pyfits.Column(name='r_mpc', format='E', array=r_mpc),
                    pyfits.Column(name='z_b', format='E', array=z_b),
                    pyfits.Column(name='ghats', format='E', array=E),
                    pyfits.Column(name='B', format='E', array=B)]
 
        manager.store('inputcat',
                      ldac.LDACCat(pyfits.BinTableHDU.from_columns(pyfits.ColDefs(cols))))

        
        
#############


def calcTangentialShear(cat, center, raCol, decCol, g1Col, g2Col):

    cluster_ra, cluster_dec = center
    ra = cat[raCol]
    dec = cat[decCol]
    e1 = cat[g1Col]
    e2 = cat[g2Col]
    
#    posangle = ((np.pi / 2.) - sphereGeometry.positionAngle(ra, dec, cluster_ra, cluster_dec)) #radians
    posangle = -((np.pi / 2.) - sphereGeometry.positionAngle(ra, dec, cluster_ra, cluster_dec)) #radians... need minus sign to get WTG result !
    r_arcmin = sphereGeometry.greatCircleDistance(ra, dec, cluster_ra, cluster_dec) * 60

    cos2phi = np.cos(2 * posangle)
    sin2phi = np.sin(2 * posangle)

    E = -(e1 * cos2phi + e2 * sin2phi)
    B = e1 * sin2phi - e2 * cos2phi

    return r_arcmin, E, B, posangle


##########


def matchById(firstcat, othercat, otherid='SeqNr', selfid='SeqNr'):
    """Returns a subset of this catalog, that matches the order of the provided catalog."""
    order = {}
    for i, x in enumerate(firstcat[selfid]):
        order[x] = i

    keeporder = []
    for x in othercat[otherid]:
        if x in order:
            keeporder.append(order[x])

    keep = np.array(keeporder)
    matched = firstcat[keep]
    print "INFO: %i matched galaxies kept" % len(matched)
    return matched
