#!/usr/bin/env python
###############
# @file ldac.py
# @author Douglas Applegate
# @date 9/2/2008
#
# @brief Utilities to make accessing LDAC cats easier
###############

from __future__ import with_statement
import unittest
import astropy.io.fits as pyfits, numpy

#######################

class MismatchedKeysException(Exception): pass

########################

class LDACCat(object):

    def __init__(self, hdu):
        self.hdu = hdu
        self.sourcefile = None
        self.hdu.data is not None   #make sure data is read into memory
        if not 'EXTNAME' in hdu.header:
            self.hdu.header['EXTNAME'] = 'OBJECTS'
        
    def __len__(self):
        return len(self.hdu.data)

    def __getitem__(self, key):

        if type(key) == type(5) or \
                type(key) == type(slice(5)):
            return self.hdu.data[key]

        if type(key) == type("a"):
            try:
                return self.hdu.data.field(key)
            except AttributeError:
                raise KeyError(key)

        raise TypeError

    def __setitem__(self, key, val):
        raise NotImplementedError

    def __delitem__(self, key):
        raise NotImplementedError

    def keys(self):
        return self.hdu.columns.names

    def __iter__(self):
        return self.hdu.data.__iter__()

    def __contains__(self, item):
        return item in self.keys()

    def has_key(self, key):
        return self.__contains__(key)

    def filter(self, mask):
        newcat = LDACCat(pyfits.BinTableHDU(data=self.hdu.data[mask]))
        newcat.hdu.header['EXTNAME'] = self.hdu.header['EXTNAME']
        newcat.sourcefile = self.sourcefile
        return newcat

    def extractColumn(self, key):

        for col in self.hdu.columns.data:
            if col.name == key:
                return col

        return None

    def saveas(self, file, clobber=False):
        self.hdu.writeto(file, clobber = clobber)
        self.sourcefile = file

    def append(self, othercat):

        if len(self.keys()) != len(othercat.keys()):
            raise MismatchedKeysException((self.keys(), othercat.keys()))
        for key in self.keys():
            if key not in othercat.keys():
                raise MismatchedKeysException((self.keys(), othercat.keys()))

        cols = [ pyfits.Column(name = key, format = format, array = numpy.hstack([self[key], othercat[key]])) \
                     for key,format in zip(self.hdu.columns.names, self.hdu.columns.formats) ]

        newcat = LDACCat(pyfits.BinTableHDU.from_columns(pyfits.ColDefs(cols)))

        for key in self.hdu.header.keys():
            if key not in newcat.hdu.header:
                newcat.hdu.header[key] = self.hdu.header[key]

        newcat.hdu.header['EXTNAME'] = self.hdu.header['EXTNAME']




        return newcat


    def matchById(self, othercat, otherid='SeqNr', selfid='SeqNr'):
        '''Returns a subset of this catalog, that matches the order of the provided catalog'''

        order = {}
        for i, x in enumerate(self[selfid]):
            order[x] = i

        keepOrder = []
        for x in othercat[otherid]:
            if x in order:
                keepOrder.append(order[x])

        keep = numpy.array(keepOrder)
        matched = self.filter(keep)
        return matched



def openObjects(hdulist, table='OBJECTS'):

    for hdu in hdulist:
        try:
            if table == hdu.header['EXTNAME']:
                return LDACCat(hdu)
        except KeyError:
            pass

    return None
    

def openObjectFile(filename, table='OBJECTS'):

    hdulist = pyfits.open(filename)
    if hdulist is None:
        return None
    cat = openObjects(hdulist, table)
    #hdulist.close()
    if cat is None:
        return None
    cat.sourcefile = filename
    return cat

def matchCommonSubsets(cat1, cat2, cat1id='SeqNr', cat2id = 'SeqNr'):

    cat1order = {}
    for i, x in enumerate(cat1[cat1id]):
        cat1order[x] = i
        
    
    cat1KeepOrder = []
    cat2Keep = []
    for x in cat2[cat2id]:
        if x in cat1order:
            cat1KeepOrder.append(cat1order[x])
            cat2Keep.append(True)
        else:
            cat2Keep.append(False)

    cat1keep = numpy.array(cat1KeepOrder)
    cat2keep = numpy.array(cat2Keep)
    cat1matched = cat1.filter(cat1keep)
    cat2matched = cat2.filter(cat2keep)
    return cat1matched, cat2matched






################################################
# TESTING
################################################

class TestComponents(unittest.TestCase):        
        
    def testExtractColumn(self):

        keys = 'FLUX_APER1 FLUXERR_APER1 BLANK1 BLANK2 MAG_APER1 MAGERR_APER1'.split()
        cols = [pyfits.Column(name = k, 
                                    format = 'E', 
                                    array = numpy.ones(30)) \
                          for k in keys]
        cat = LDACCat(pyfits.BinTableHDU.from_columns(pyfits.ColDefs(cols)))

        extractedCol = cat.extractColumn('BLANK1')
        self.assertTrue((extractedCol.array == cols[2].array).all())
        self.assertEqual(extractedCol.name, cols[2].name)
        self.assertEqual(extractedCol.format, cols[2].format)

    #####

    def testAppend_Mismatch(self):

        keys = 'a b c'.split()
        cols = [pyfits.Column(name = k, 
                              format = 'E', 
                              array = numpy.ones(30)) \
                    for k in keys]
        cat = LDACCat(pyfits.BinTableHDU.from_columns(pyfits.ColDefs(cols)))

        keys = 'b e'.split()
        cols = [pyfits.Column(name = k, 
                              format = 'E', 
                              array = numpy.ones(30)) \
                    for k in keys]
        cat2 = LDACCat(pyfits.BinTableHDU.from_columns(pyfits.ColDefs(cols)))

        self.assertRaises(MismatchedKeysException, lambda : cat.append(cat2))

    ######

    def testAppend(self):

        keys = 'a b c'.split()
        cols = [pyfits.Column(name = k, 
                              format = 'E', 
                              array = numpy.ones(30)) \
                    for k in keys]
        cat = LDACCat(pyfits.BinTableHDU.from_columns(pyfits.ColDefs(cols)))
        cat.hdu.header['A'] = 25
        cat.hdu.header['B'] = 'sss'
        cat.hdu.header['EXTNAME'] = 'STUFF'

        keys = 'a b c'.split()
        cols = [pyfits.Column(name = k, 
                              format = 'E', 
                              array = numpy.zeros(30)) \
                    for k in keys]
        cat2 = LDACCat(pyfits.BinTableHDU.from_columns(pyfits.ColDefs(cols)))
        cat2.hdu.header['EXTNAME'] = 'OBJECTS'
        cat2.hdu.header['A'] = 27
        
        cat3 = cat.append(cat2)

        self.assertEquals(len(cat3), len(cat2) + len(cat))

        self.assertEquals(cat3.keys(), cat.keys())

        testcolumn = numpy.hstack([numpy.ones(30), numpy.zeros(30)])

        for key in cat3.keys():
            self.assertTrue((testcolumn == cat3[key]).all())


        self.assertEquals(cat3.hdu.header['EXTNAME'], 'STUFF')
        self.assertEquals(cat3.hdu.header['A'], 25)
        self.assertEquals(cat3.hdu.header['B'], 'sss')


    ########

    def testMatchById(self):

        cat1 = LDACCat(pyfits.BinTableHDU.from_columns(pyfits.ColDefs([pyfits.Column(name = 'SeqNr',
                                                                      format = 'K',
                                                                      array = numpy.arange(30))])))
        cat2 = LDACCat(pyfits.BinTableHDU.from_columns(pyfits.ColDefs([pyfits.Column(name = 'SeqNr',
                                                                      format = 'K', \
                               array = numpy.random.permutation(numpy.arange(20, 40)))])))

        matchedCat = cat1.matchById(cat2)

        self.assertTrue((cat2['SeqNr'][cat2['SeqNr'] < 30] == matchedCat['SeqNr']).all())

        
        matchedCat = cat2.matchById(cat1)

        self.assertTrue((cat1['SeqNr'][cat1['SeqNr'] >= 20] == matchedCat['SeqNr']).all())

    ############

    def testMatchCommonSubsets(self):

        cat1 = LDACCat(pyfits.BinTableHDU.from_columns(pyfits.ColDefs([pyfits.Column(name = 'SeqNr',
                                                                      format = 'K',
                                                                      array = numpy.arange(30))])))
        cat2 = LDACCat(pyfits.BinTableHDU.from_columns(pyfits.ColDefs([pyfits.Column(name = 'SeqNr',
                                                                      format = 'K', \
                               array = numpy.random.permutation(numpy.arange(20, 40)))])))


        
        matched1, matched2 = matchCommonSubsets(cat1, cat2)

        self.assertTrue((matched1['SeqNr'] == matched2['SeqNr']).all())

        


                     

        

#######################################
        
def test():
    
    testcases = [TestComponents]
    
    suite = unittest.TestSuite(map(unittest.TestLoader().loadTestsFromTestCase,
                                   testcases))
    unittest.TextTestRunner(verbosity=2).run(suite)



#####################################################
# COMMAND LINE RUNNABLE
#####################################################

if __name__ == '__main__':

    test()
