###################################
# File handler for the compact, single fits-file storage format
# Assumes that all relevant cuts have been applied to the file already
####################################

import astropy.io.fits as pyfits
import varcontainer, ldac, pdzfile_utils

####################################

def createCompactFile(clustername, manager, outfile, clobber = False):

    inputcat = manager.inputcat
    header = inputcat.hdu.header
    header.update('cluster', clustername)
    header.update('zcluster', manager.zcluster)
    header.update('psfsize', manager.psfsize)
    header.update('pixscale', manager.pixscale)

    pdzcat = pdzfile_utils.createPDZcat(inputcat['SeqNr'], manager.pdzrange, manager.pdz)
    pdzcat.hdu.header.update('EXTNAME', 'PDZ')

    hdulist = pyfits.HDUList([pyfits.PrimaryHDU(), inputcat.hdu, pdzcat.hdu])

    hdulist.writeto(outfile, clobber = clobber)

    


####################################

class CompactFilehandler(object):

    ######

    def __init__(self):

        self.cuts =  []

    ######

    def addCLOps(self, parser):

        parser.add_option('-i', '--inputcat',
                          help = 'path and name of input catalog to use')

    ######

    def createOptions(self, inputcat, options = None, args = None):

        if options is None:
            options = varcontainer.VarContainer()

        options.inputcat = inputcat

        return options, args

    ######

    def readData(self, manager):

        options = manager.options

        manager.open('inputcat', options.inputcat, ldac.openObjectFile)
        inputcat = manager.inputcat

        manager.clustername = inputcat.hdu.header['cluster']
        manager.zcluster = inputcat.hdu.header['zcluster']
        manager.psfsize = inputcat.hdu.header['psfsize']
        manager.pixscale = inputcat.hdu.header['pixscale']

        manager.open('pdzmanager', options.inputcat, pdzfile_utils.PDZManager.open, 'PDZ')
        manager.store('pdzrange pdz'.split(), manager.pdzmanager.associatePDZ(manager.inputcat['SeqNr']))

        manager.replace('pdzmanager', None)
