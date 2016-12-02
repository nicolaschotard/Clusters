#############################
# Handles loading files for a simulation run
#############################

import numpy as np
import astropy.io.fits as pyfits
from . import ldac
from . import pdzfile_utils
import shearprofile as sp


#############################

class GeneralFilehandler(object):

    def addCLOps(parser):

        parser.add_option('-l', '--lensingcat', dest='lensingcat',
                          help='Lensing catalog with shear information')
        parser.add_option('-b', '--bpzfile', dest='bpzfile',
                          help='BPZ file containing summary photo-z info')
        parser.add_option('-p', '--pdzfile', dest='inputPDZ',
                          help='Parsed BPZ pdz file')
        parser.add_option('--psf-size', dest='psfsize',
                          type='float',
                          help='Size of the PSF in pixels')
        parser.add_option('-z', '--zcluster', dest='zcluster',
                          type='float',
                          help='Cluster redshift')
        parser.add_option('--pixscale', dest='pixscale',
                          type='float',
                          help='Pixel scale in arcseconds / pix')
        parser.add_option('--r500', dest='r500',
                          help='R500 value; for translating results into a mass',
                          type='float')
        parser.add_option('--xcol', dest='xcol', default='Xpos',
                          help='X pixel location column in the lensing catalog')
        parser.add_option('--ycol', dest='ycol', default='Ypos',
                          help='Y pixel location column in the lensing catalog')
        parser.add_option('--g1col', dest='g1col', default='gs1',
                          help='shear1 column in the lensing catalog')
        parser.add_option('--g2col', dest='g2col', default='gs2',
                          help='shear2 column in the lensing catalog')
        parser.add_option('--sizecol', dest='sizecol', default='rh',
                          help='Column of size of objects in pixels')
        parser.add_option('--centerx', dest='centerx', default=5000,
                          type='float',
                          help='Cluster x center in pixels')
        parser.add_option('--centery', dest='centery', default=5000,
                          type='float',
                          help='Cluster y center in pixels')
        parser.add_option('--snratiocol', dest='snratio', default='snratio',
                          help='Snratio column name')
        parser.add_option('--zbcol', dest='zbcol', default='BPZ_Z_B',
                          help='Name of z best column in the BPZ catalog')


    #############################

    def readData(manager):

        options = manager.options
        args = manager.args

        manager.zcluster = options.zcluster
        manager.r500 = options.r500
        manager.psfsize = options.psfsize
        manager.pixscale = options.pixscale


        if 'lensingcat' not in manager:
            manager.open('lensingcat', options.lensingcat, ldac.openObjectFile)

        r_arc, E, B = sp.calcTangentialShear(cat=manager.lensingcat,
                                             center=(options.centerx, options.centery),
                                             pixscale=options.pixscale,
                                             xcol=options.xcol,
                                             ycol=options.ycol,
                                             g1col=options.g1col,
                                             g2col=options.g2col)

        r_pix = r_arc / options.pixscale

        r_mpc = r_arc * (1. / 3600.) * (np.pi / 180. ) * sp.angulardist(options.zcluster)

        size = manager.lensingcat[options.sizecol] / options.psfsize
        snratio = manager.lensingcat[options.snratio]

        manager.open('bpzcat', options.bpzfile, ldac.openObjectFile, 'STDTAB')

        manager.store('matched_bpzcat', manager.bpzcat.matchById, manager.lensingcat)
        z_b = manager.matched_bpzcat['BPZ_Z_B']
        odds = manager.matched_bpzcat['BPZ_ODDS']
        nfilt = manager.matched_bpzcat['NFILT']
        galtype = manager.matched_bpzcat['BPZ_T_B']

        cols = [pyfits.Column(name='SeqNr', format='J', array=manager.lensingcat['SeqNr']),
                pyfits.Column(name='r_mpc', format='E', array=r_mpc),
                pyfits.Column(name='size', format='E', array=size),
                pyfits.Column(name='snratio', format='E', array=snratio),
                pyfits.Column(name='z_b', format='E', array=z_b),
                pyfits.Column(name='z_t', format='E', array=galtype),
                pyfits.Column(name='odds', format='E', array=odds),
                pyfits.Column(name='nfilt', format='I', array=nfilt),
                pyfits.Column(name='ghats', format='E', array=E),
                pyfits.Column(name='B', format='E', array=B)]

        manager.store('inputcat', ldac.LDACCat(pyfits.BinTableHDU.from_columns(pyfits.ColDefs(cols))))

        manager.open('pdzmanager', options.inputPDZ, pdzfile_utils.PDZManager.open)
        manager.store('pdzrange pdz'.split(),
                      manager.pdzmanager.associatePDZ(manager.lensingcat['SeqNr']))

        manager.replace('pdzmanager', None)
        manager.replace('matched_bpz', None)
        manager.replace('bpzcat', None)
