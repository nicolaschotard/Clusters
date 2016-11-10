#########################
# Utilities to deal with Clusters pdz file ops
#########################

import sys, re, cPickle
import numpy as np
import astropy.table



##########################

class PDZManager(object):

    def __init__(self, pdzcat):

        self.idkey = 'objectId'

        self.pdzcat = pdzcat
        self.pdzrange = None
        self.index = None

        self._buildIndex()

    #################


    def __getitem__(self, key):
        #key is the object ID
        #returns pdz array
        #index maps ID to element number

        if self.index is None or key not in self.index:
            raise IndexError

        return self.pdzcat[self.index[key]][1]
    
    ####################

    def _buildIndex(self):

        self.index = {}
        for i, id in enumerate(self.pdzcat[self.idkey]):
            self.index[id] = i

    #####################


    def associatePDZ(self, z_ids):

        arrangedPDZ = []
        for id in z_ids:
            arrangedPDZ.append(self[id])

        finalPDZs = np.column_stack(arrangedPDZ).transpose()

        return self.pdzrange, finalPDZs


    ######################

    @classmethod
    def open(cls, pdzfile):

        pdzrange = astropy.table.Table.read(pdzfile, path='pdz_bins')
        pdzvals = astropy.table.Table.read(pdzfile, path='pdz_values')
        
        return cls(pdzrange, rawpdzvals)    

###############################


        


