#!/usr/bin/env python
#################################################

from __future__ import with_statement
import tempfile, os, re, unittest, datetime
import ephem, numpy as np, scipy
import astropy.io.fits as pyfits
import ldac

#################################################

################################################################
# Exceptions
################################################################

class FluxKeyException(Exception): pass

class NoPhotometricCalibrationException(Exception): pass

class UnrecognizedFilterException(Exception): pass

##################################################################


def getTempFile(dir='/tmp', prefix='',suffix=''):
    thefile, thefilename = tempfile.mkstemp(dir=dir,
                                          prefix=prefix,
                                          suffix=suffix)
    os.close(thefile)
    return thefilename


#################################

def compare(x, y):
    if x['SEEING'] < y['SEEING']:
        return 1    
    
    elif x['SEEING'] == y['SEEING']:
        return 0
    else:
        return -1

##################################

_stdfiltermatch = re.compile('(\w+)-(.+?)-(\d)-(.+)')
def parseFilter(filter):

    match = _stdfiltermatch.match(filter)
    if match is None:
        raise UnrecognizedFilterException(filter)
    instrum = match.group(1)
    config = match.group(2)
    chipid = int(match.group(3))
    stdfilter = match.group(4)

    return instrum, config, chipid, stdfilter

####################################
# EXTINCTIONS
####################################

subaru_extinctions = {
    'W-J-U' : 5.434,   #guess
    'W-J-B' : 4.031,
    'W-J-V' : 3.128,
    'W-C-RC' : 2.548,
    'W-C-IC' : 1.867,
    'W-S-I+' : 1.867,
    'W-S-Z+' : 1.481
}

#from http://www.astro.princeton.edu/~schlegel/dust/dustpub/CodeIDL/README.IDL
megaprime_extinctions = {
    'u' : 5.155,
    'g' : 3.793,
    'r' : 2.751,
    'i' : 2.086,
    'z' : 1.479
}

special_extinctions = {
    'B' : 4.315,
    'I' : 1.867,
    'K' : 0.367
}

wht_extinctions = {
    'B' : 4.315,
    'U' : 5.434

}

extinctions = {
    'SUBARU' : subaru_extinctions,
    'MEGAPRIME' : megaprime_extinctions,
    'SPECIAL' : special_extinctions,
    'WHT' : wht_extinctions
}


def getExtinction(filter):
    #expect INSTRUM-CONFIG-CHIPID-FILTER filters

    instrum, config, chipid, stdfilter = parseFilter(filter)

    return extinctions[instrum][stdfilter]

##################################

flux_search = re.compile('^FLUX_')
fluxerr_search = re.compile('^FLUXERR_')
mag_search = re.compile('^MAG_')
magerr_search = re.compile('MAGERR_')
rejectKeys = "FLUX_RADIUS".split()
def sortFluxKeys(keylist):

    fluxkeys = []
    fluxerrkeys = []
    temp_magonlykeys = []
    otherkeys = []
    for key in keylist:
        if flux_search.match(key) and key not in rejectKeys:
            fluxkeys.append(key)
        elif fluxerr_search.match(key):
            fluxerrkeys.append(key)
        elif mag_search.match(key):
            temp_magonlykeys.append(key)
        elif magerr_search.match(key):
            continue
        else:
            otherkeys.append(key)

    magonlykeys = []
    for key in temp_magonlykeys:
        match = re.match('^MAG_(.+)', key)
        magtype = match.group(1)
        fluxkey = 'FLUX_%s' % magtype
        if fluxkey not in keylist:
            print fluxkey
            magonlykeys.append(key)
            magerr_key = 'MAGERR_%s' % magtype
            otherkeys.append(magerr_key)

    return fluxkeys, fluxerrkeys, magonlykeys, otherkeys

#########################################

_fluxsplit = re.compile('^FLUX_(\w+)(?:-(.+))?')

def extractFluxType(fluxkey):

    match = _fluxsplit.match(fluxkey)
    if match is None:
        raise FluxKeyException('Cannot parse fluxkey: %s' % fluxkey)

    suffix = match.group(1)
    return suffix

############################################

_magsplit = re.compile('^MAG_(\w+)(?:-(.+))?')

def extractMagType(magkey):

    match = _magsplit.match(magkey)
    if match is None:
        raise FluxKeyException('Cannot parse magkey: %s' % magkey)

    suffix = match.group(1)
    return suffix

############################################

def calc_seeing(file,PIXSCALE):
    #set up bins
    binsize = 0.03 
    nbins = int((3.0-0.3)/binsize+0.5)
    import scipy
    bin = scipy.zeros(nbins)


    # for each line get fwhm
    for line in open(file,'r').readlines():
        tokens = line.split()
        fwhm_obj = float(tokens[2])
        flag = float(tokens[3])

        
        # make sure flag is zero and the seeing is reasonable
        if 3.0 > fwhm_obj*PIXSCALE > 0.3 and flag == 0:
            
            actubin = int((fwhm_obj * PIXSCALE - 0.3)/binsize)
            bin[actubin] += 1
    # find max
    max = 0
    k = 0
    nobjs = 0
    for i in range(nbins):
        nobjs += bin[i]
        if bin[i]>max:
            k=i
            max = bin[i]

    # set the fwhm
    fwhm = 0.3 + k*binsize

    # check that its ok
    if nobjs < 100:
        fwhm = -999
        
    # print to screen
      
    return fwhm

###############################

def run(command,to_delete=[]):
    import os
    for file in to_delete: 
        os.system('rm ' + file)
    print command
    #raw_input()
    os.system(command)



def make_filters_info(filters):
    filters_dat = [['B','b',4.031,0],['W-J-B','b',4.031,0],['W-J-V','v',3.128,1],['W-C-RC','r',2.548,2],['W-C-IC','i',1.867,3],['W-S-I+','i',1.867,3],['I','i',1.867,3],['W-S-Z+','z',1.481,4]]

    filters_info = []
    for filter in filters:
    	for filter_dat in filters_dat:
    		if filter_dat[0] == filter:
    			filters_info.append(filter_dat)
   
    return filters_info

############################

def get_header_info(file): 
    # Gets the appropriate header info for image file 

    header = pyfits.getheader(file)

    GAIN = header['GAIN']

    if header.has_key('CDELT1'):
        CDELT1 = header['CDELT1']
        CDELT2 = header['CDELT2']
    else:
        CDELT1 = header['CD1_1']
        CDELT2 = header['CD2_2']

    EXPTIME = header['EXPTIME']

    PIXSCALE = (abs(float(CDELT1)) + abs(float(CDELT2))) / (2.0) * 3600.0
    
    return GAIN, PIXSCALE, EXPTIME 

################################

def get_header_kw(file,kws=[]):
    import commands, string
    print kws, file
    print 'kws'
    return_mat = {} 
    for kw in kws + ['GAIN','CD1_1','CD2_2','CD2_1','CD1_2','CDELT1','CDELT2']:
        command = "dfits " + file + " | fitsort -d " + kw + " | awk '{print $2}' "
        #print command
        value = commands.getoutput(command)     
        return_mat[kw] = value

    if string.find(return_mat['CDELT1'],'/') == -1 and len(return_mat['CDELT1'].replace(' ','')) > 1:
        PIXSCALE = (abs(float(return_mat['CDELT1'])) + abs(float(return_mat['CDELT2']))) / (2.0) * 3600.0

        return_mat['PIXSCALE'] = PIXSCALE
    elif string.find(return_mat['CD1_1'],'/') == -1:
        PIXSCALE = (abs(float(return_mat['CD1_1'])) + abs(float(return_mat['CD2_2']))) / (2.0) * 3600.0     
        if PIXSCALE == 0: 
            PIXSCALE = (abs(float(return_mat['CD2_1'])) + abs(float(return_mat['CD1_2']))) / (2.0) * 3600.0

        return_mat['PIXSCALE'] = PIXSCALE

    return return_mat 


########################################

def convert_to_galactic(alpha,delta):

    gallongs = []
    gallats = []
    for coord_in_ra, coord_in_dec in zip(alpha, delta):

        coord = ephem.Equatorial( str(coord_in_ra*(24./360.)), str(coord_in_dec), epoch='2000') # input needs to be in HOURS as a STRING
        g = ephem.Galactic(coord, epoch='2000') # output is in degrees not hours--it's latitude/longitude
 
        spt = re.split('\:',str(g.lat))

        gallat = float(spt[0]) / abs(float(spt[0])) * (abs(float(spt[0])) + float(spt[1])/60. + float(spt[2])/3600. )
        spt = re.split('\:',str(g.long))

        gallong = float(spt[0]) / abs(float(spt[0])) *  (abs(float(spt[0])) + float(spt[1])/60. + float(spt[2])/3600. )
        
        gallongs.append(gallong)
        gallats.append(gallat)

    return scipy.array(gallongs), scipy.array(gallats)

#######################################

def getDust(alpha, delta):

    dustpos_filename = getTempFile(dir='/tmp', prefix='dustpos')
    ebvgal_filename = getTempFile(dir='/tmp', prefix='ebvgal')
    
    def cleanUp():

        if os.path.exists(dustpos_filename):
            os.remove(dustpos_filename)
        if os.path.exists(ebvgal_filename):
            os.remove(ebvgal_filename)

    #############
        
    try:

        gallong, gallat = convert_to_galactic(alpha, delta)

        with open(dustpos_filename, 'w') as dustpos_file:
            for pos in zip(gallong, gallat):
                dustpos_file.write('%f %f\n' % pos)
            
        os.system('dust_getval interp=y ipath=/nfs/slac/g/ki/ki04/pkelly/DUST/maps infile=%s outfile=%s noloop=y' % (dustpos_filename, ebvgal_filename)) 

        _isFileCorrupted = not (os.path.exists(ebvgal_filename) and os.path.getsize(ebvgal_filename) > 0)
        if _isFileCorrupted:
            raise RuntimeError('Dust file corrupted or not produced!')

        with open(ebvgal_filename) as ebvgal_file:
            ebv = scipy.array([float(lines.split()[2]) for lines in ebvgal_file.readlines()])

        return ebv


    finally:
        cleanUp()
    
#############################

def dust_correct(table,tempfile, tempdir='/tmp', clean=False):
    import os,re
    ebvfile = getTempFile(dir=tempdir,prefix='outebv.')
    ebvgal = getTempFile(dir=tempdir,prefix='outgal.')
    tempmap = getTempFile(dir=tempdir,prefix='tempmap.')
    templines = getTempFile(dir=tempdir,prefix='templines.')

    command = "ldactoasc -b  -i " + table + " -t OBJECTS -k ALPHA_J2000 DELTA_J2000 > " + ebvfile 
    print command
    os.system(command)
    print 'done'

    # convert coordinates to galactic coordinates
    file = open(ebvgal,'w')
    print ebvfile
    for line in open(ebvfile,'r').readlines():
        res = re.split('\s+',line)
        if res[0] == '': res = res[1:]
        #print res
        gallong, gallat = convert_to_galactic(float(res[0]),float(res[1]))
        file.write(str(gallong) + ' ' + str(gallat) + '\n') 
    file.close()

    print 'running dust'
    os.system('dust_getval interp=y ipath=/nfs/slac/g/ki/ki04/pkelly/DUST/maps infile=' + ebvgal + ' outfile=' + tempmap + ' noloop=y') 
    
    if os.path.exists(tempmap) and os.path.getsize(tempmap) > 0:
        print 'dust done'
    else:
        raise RuntimeError('Dust file corrupted or not produced!')

    command = "gawk '{print $3}' " + tempmap + " > " + templines 
    print command
    os.system(command)
    command = "asctoldac -i " + templines + " -o " + tempfile + " -c " + "./photconf/EBV.conf -t OBJECTS "
    print command
    os.system(command)
    command = "ldacjoinkey -o " + tempfile + " -i " + table + " -p " + tempfile + " -t OBJECTS -k EBV" 
    print command
    os.system(command)

    if clean:
        try:
            os.remove(ebvfile)
            os.remove(ebvgal)
            os.remove(tempmap)
            os.remove(templines)
        except OSError:
            pass

def color_index(color, colors):
    i = 0
    for c in colors:
        if c == color: 
            break
        i += 1
    return i


def convert_modelname_to_array(model_name):
    ''' convert the model name into a list of dictionaries with the terms and names '''
    import re
    print model_name
    red = re.split('P',model_name)
    model = []  
    for term in red:
        redterm = re.split('T',term)
        model.append({'name':term,'term':redterm})
    print model

    return model

def apply_color_term(model,calib,table,filter,magnitude):
    print calib
    dict = {'zp':1.,'bandcomp':table.field(calib['bandcomp%' +filter] + 'mag'),'color1':table.field(calib['color1which%'+filter]),'color2':table.field(calib['color2which%'+filter])}

    print magnitude + '_' + filter 
    from copy import copy 
    data = copy(dict['bandcomp']) # table.field(magnitude + '_' + filter)

    print data[0:10]

    
    for i in range(len(model)):    
        if model[i]['name'] != 'zp':
            term = model[i]['term']                                                                             
            term_name = reduce(lambda x,y: x  + 'T' + y,term)
            print term_name, term, calib[term_name + '%' + filter]
            data -= float(calib[term_name + '%' + filter]) * reduce(lambda x,y: x * y,[dict[z] for z in term]) 

    print data[0:10]

    return data

def zp_correct(model,calib,table,filter,magnitude):
    print calib
    dict = {'zp':1.}

    print magnitude + '_' + filter 
    data = table.field( magnitude + '_' + filter)

    print data[0:10]
    
    for i in range(len(model)):    
        if model[i]['name'] == 'zp':
            term = model[i]['term']                                                                             
            term_name = reduce(lambda x,y: x  + 'T' + y,term)
            print term_name, term, calib[term_name + '%' + filter]
            data += float(calib[term_name + '%' + filter]) * reduce(lambda x,y: x * y,[dict[z] for z in term]) 
    print data[0:10]
    return data

''' remove the color term from standard sdss magnitudes '''
def color_std_correct(model,calib,table,filter,magnitude,color1which):

    dict = {'zp':1.,'color1':table.field(color1which)}

    print magnitude + '_' + filter , '$$$$'
    from copy import copy 
    data = copy(table.field(magnitude))

    return data











    # note this!!!!!!!!!!!!!!! what the heck is going on???
    print data[0:10]
   
    print model, '$$$$$$$'

    for i in range(len(model)):    
        print model[i]['name']
        if model[i]['name'] != 'zp':
            term = model[i]['term']  
            term_name = reduce(lambda x,y: x  + 'T' + y,term)
            print term_name, term, calib[term_name+'_star_'], '!!!!'
            data += float(calib[term_name + '_star_']) * reduce(lambda x,y: x * y,[dict[z] for z in term])
    print data[0:10]

    print data[0:10], table.field(magnitude)[0:10], table.field(color1which)[0:10]
    return data

def zp_dust_correct(model,calib,input_data,filter,magnitude,color='yes',doDust=True):


    filters = {'B':4.031,'W-J-B':4.031,'W-J-V':3.128,'W-C-RC':2.548,'I':1.867,'W-C-IC':1.867,'W-S-Z+':1.481}
    #ab_filters = {'B':4.031,'W-J-B':4.031,'W-J-V':3.128,'W-C-RC':2.548,'I':1.867,'W-C-IC':1.867,'W-S-Z+':1.481}


    print calib
    if doDust:
        dict = {'zp':1.,'EBV':input_data.field('EBV')}
    else:
        dict = {'zp':1}

    print magnitude + '_' + filter 

    ''' allow input to be a table or a single value '''
    if type(input_data) != type(3) and color=='yes':
        data = input_data.field( magnitude + '_' + filter)
        thresh = input_data.field('MU_THRESHOLD_' + filter)
        print data[0:10]
    elif type(input_data) != type(3) and color!='yes':
        data = input_data.field( magnitude)
        thresh = input_data.field('MU_THRESHOLD')
        print data[0:10]
    else:
        data = input_data
        print data

    for i in range(len(model)):    
        if model[i]['name'] == 'zp':
            term = model[i]['term']                                                                             
            term_name = reduce(lambda x,y: x  + 'T' + y,term)
            print 'a:',term_name, term, calib[term_name + '%' + filter]

#            print data[0:50]
#            print thresh[0:50]
            data += float(calib[term_name + '%' + filter]) * reduce(lambda x,y: x * y,[dict[z] for z in term]) 
            thresh += float(calib[term_name + '%' + filter]) * reduce(lambda x,y: x * y,[dict[z] for z in term]) 
#            print data[0:50]
#            print thresh[0:50]
#            print filter

        if doDust:
            data -= dict['EBV'] * filters[filter] # dust 
            thresh -= dict['EBV'] * filters[filter] # dust 

    if type(input_data) != type(3):
        print data[0:10]
    else:
        print data

    return data, thresh

def get_calibrations(cluster,filters,aper_index = 0, vars_fit=None):
    #aper_index aperature size number from original SE extraction
    #vars_fit? whats it used for?
    #chanings filters to mean the INSTRUM-CONFIG-chipid-filter type filters
    import MySQLdb
    db2 = MySQLdb.connect(db='subaru', user='weaklensing', passwd='darkmatter', host='ki-sr01')
    c = db2.cursor()
    c.execute("describe photometry_db") 
    results = c.fetchall()
    columns = []
    for column in results:
       columns.append(column[0])
    cal_dict = {} 
    for filter in filters:
        
        command = "select * from photometry_db where BONN_TARGET='%(cluster)s' AND fitID='%(filter)_%(aper)d'" % {'cluster' : cluster,'filter' : filter, 'aper' : aper_index}

        print command
        c.execute(command)
        results = c.fetchall()
        if len(results) > 0:

            use_phot = results[-1]
            result_lines = []
            for result in results:
                line = {}
                for i in range(len(columns)):
                    line[columns[i]] = result[i]
                    result_lines.append(line)
                    use = result_lines[-1]
            
            for key in use.keys(): 
                cal_dict[key.replace('_star_','') + '%' + filter] = use[key]
                print key, use[key]
            cal_dict['model_name%' + filter] = 'zpPcolor1'
    return cal_dict

def get_calibrations_threesecond(cluster,filters,vars_fit=[2,2,3,2,2]):
    import MySQLdb
    db2 = MySQLdb.connect(db='subaru', user='weaklensing', passwd='darkmatter', host='ki-sr01')
    c = db2.cursor()
    c.execute("describe photometry_db") 
    results = c.fetchall()
    columns = []
    for column in results:
       columns.append(column[0])
        
    #filters = [['B','b',4.031],['W-J-B','b',4.031],['W-J-V','v',3.128],['W-C-RC','r',2.548],['W-C-IC','i',1.867],['W-S-Z+','z',1.481]]
       
    cal_dict = {} 
    ip = 0
    print filters
    for filter in filters:
        print ip
        var=vars_fit[ip]
        ip += 1
        long_filter = filter[0]
        short_filter = filter[1]
        extcoeff = filter[2]
        
        command = "select * from photometry_db where BONN_FILTER='" + str(long_filter) + "' AND BONN_TARGET='" + str(cluster) + "'" # AND cluster!='illumination'" 

        #command = "select * from photometry_db where filter like '" + str(long_filter) + "%' and cluster='MACS0417-11' and  OBJECT like '%0417c%'" 
        print command
        c.execute(command)
        results = c.fetchall()
        if len(results) > 0:
            #print results                                                                                                                    
            use_phot = results[-1]
            result_lines = []
            for result in results:
                line = {}
                for i in range(len(columns)):
                    line[columns[i]] = result[i]
                    result_lines.append(line)
                    use = result_lines[-1]
            
                    #print use
                       
            for key in use.keys(): 
                cal_dict[key.replace('_star_','') + '%' + long_filter] = use[key]
                print key, use[key]
            cal_dict['model_name%' + long_filter] = 'zpPcolor1'
        #print cal_dict['zp%W']
    return cal_dict



##########################################
# TESTING
##########################################

class TestComponents(unittest.TestCase):

    def testSortFluxKeys(self):

        keys = 'FLUX_APER1 FLUXERR_APER1 NPIXELS MAG_FAKE MAGERR_FAKE MAG_APER1 MAGERR_APER1 FLUX_AUTO2 FLUXERR_AUTO2'.split()
        fluxkeys, fluxerrkeys, magonlykeys, otherkeys = sortFluxKeys(keys)
        self.assertEquals(fluxkeys, 'FLUX_APER1 FLUX_AUTO2'.split())
        self.assertEquals(fluxerrkeys, 'FLUXERR_APER1 FLUXERR_AUTO2'.split())
        self.assertEquals(magonlykeys, 'MAG_FAKE'.split())
        self.assertEquals(otherkeys, 'NPIXELS MAGERR_FAKE'.split())

    ################

    def testSortFluxKeys_badones(self):

        keys = "FLUX_APER FLUX_APER1 MAG_APER FLUX_KRON FLUX_AUTO FLUX_RADIUS".split()
        fluxkeys, fluxerrkeys, magonlykeys, otherkeys = sortFluxKeys(keys)
        self.assertEquals(fluxkeys, "FLUX_APER FLUX_APER1 FLUX_KRON FLUX_AUTO".split())
        self.assertEquals(magonlykeys, [])
        self.assertEquals(otherkeys, ["FLUX_RADIUS"])

    #################

    def testExtractFluxType(self):

        errkey = extractFluxType('FLUX_APER1')
        self.assertEquals(errkey, 'APER1')

    def testExtractFluxType_Long(self):
        errkey = extractFluxType('FLUX_APER1-SUBARU-COADD-1-W-J-V')
        self.assertEquals(errkey, 'APER1')


    def testConstructFluxErrKey_Badkey(self):

        self.assertRaises(FluxKeyException, 
                          extractFluxType, 'FAKE1_FAKE')



    #################

    def testConvertToGalactic(self):

        gallong, gallat = convert_to_galactic(scipy.array([0]), scipy.array([0]))

        self.assertEqual(gallong.shape, (1,))
        self.assertEqual(gallat.shape, (1,))

        print gallong
        print gallat
        self.assertTrue( (np.abs(gallong - 96.3373) < 1e-4).all() )
        self.assertTrue( (np.abs(gallat + 60.1886) < 1e-4).all() )

    ######################

    def testGetDust(self):

        alpha = scipy.random.standard_normal(30)
        delta = scipy.random.standard_normal(30)

        ebv  = getDust(alpha, delta)

        self.assertTrue(ebv.shape, (30,))

    ########################

    def testParseFilter(self):

        instrum, config, chipid, stdfilter = parseFilter('SUBARU-10_2-1-W-C-RC')
        self.assertEqual(instrum, 'SUBARU')
        self.assertEqual(config, '10_2')
        self.assertEqual(chipid, 1)
        self.assertEqual(stdfilter, 'W-C-RC')

        instrum, config, chipid, stdfilter = parseFilter('MEGAPRIME-0-1-u')
        self.assertEqual(instrum, 'MEGAPRIME')
        self.assertEqual(config, '0')
        self.assertEqual(chipid, 1)
        self.assertEqual(stdfilter, 'u')

    def testParseFilter_UnrecognizedFilterException(self):

        self.assertRaises(UnrecognizedFilterException, parseFilter, 'RANDOM-V')

    #########################


    ############################

    def testGetExtinction(self):

        extinction = getExtinction('SUBARU-10_2-1-W-C-RC')
        self.assertEqual(type(extinction), type(5.))

    ############################

    def testGetExtinction_Coadd(self):

        extinction = getExtinction('SUBARU-COADD-1-W-J-V')
        self.assertEquals(extinction, 3.128)

        extinction = getExtinction('MEGAPRIME-COADD-1-u')
        self.assertEquals(extinction, 5.155)

    ############################




##############################################
        


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
