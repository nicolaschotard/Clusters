###################################
# File handler for an astropy table file format
# Assumes that all relevant cuts have been applied to the file already
####################################

import astropy.table as table
import varcontainer

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
                      sheartable,
                      pdzfile,
                      psfsize = None,
                      pixscale = None,
                      options = None, args = None):

        if options is None:
            options = varcontainer.VarContainer()

        options.cluster = cluster
        options.zcluster = zcluster
        options.sheartable = sheartable
        options.pdzfile = pdzfile
        options.psfsize = psfsize
        options.pixscale = pixscale

        return options, args

    ######

    def readData(self, manager):

        options = manager.options

        manager.sheartable = options.sheartable
        
        manager.clustername = options.cluster
        manager.zcluster    = options.zcluster
        manager.psfsize     = options.psfsize
        manager.pixscale    = options.pixscale


        manager.open('pdzmanager', options.pdzfile, pdzfile_utils.ClustersPDZManager.open)
        manager.store('pdzrange pdz'.split(), manager.pdzmanager.associatePDZ(manager.sheartable['objectId']))

        manager.replace('pdzmanager', None)


