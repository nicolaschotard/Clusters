#!/usr/bin/env python

import os, yaml, cPickle
import numpy as N
from argparse import ArgumentParser

from Clusters import data

import lsst.daf.persistence as dafPersist
import lsst.afw.geom as afwGeom

if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument('config', help='Configuration (yaml) file')
    parser.add_argument("--output",
                        help="Name of the output file (pkl file)")
    args = parser.parse_args()

    config = yaml.load(open(args.config))
    if args.output is None:
        args.output = os.path.basename(args.config).replace('.yaml', '_output.pkl')

    filters = config['filters']
    
    raClust = afwGeom.degToRad(config['ra'])
    deClust = afwGeom.degToRad(config['dec'])

    print "INFO: Working on cluster %s (z=%.4f)" % (config['cluster'],
                                                    config['redshift'])
    print "INFO: Working on filters", filters
    print "INFO: Butler located under %s" % config['butler']
    
    # Initialize butler
    print "INFO: Initializing the butler..."
    butler = dafPersist.Butler(config['butler'])
    print "INFO: done."
    
    # Initialize some lists
    print "INFO: Initializing variables"
    mags = {f: [] for f in filters}
    mags_sigma = {f: [] for f in filters}
    ell = {f: {'e1': [], 'e2': []} for f in filters}
    coords = {'ra': [], 'dec': [], 'id': []}
    resolution = {f: [] for f in filters}
    xSrc, ySrc = [], []
    
    print "INFO: Looping over the input patches"
    for i, patch in enumerate(config['patches']):
    
        print "INFO:   Working on patch", patch
        meas = {f: data.get_from_butler(butler, 'deepCoadd_meas', f, patch) for f in filters} # meas_
        forced = {f: data.get_from_butler(butler, 'deepCoadd_forced_src', f, patch) for f in filters}
        calexp = {f: data.get_from_butler(butler, 'deepCoadd_calexp', f, patch) for f in filters} # md_
        calib = {f: calexp[f].getCalib() for f in filters}
    
        if i == 0:
            schema_m = meas['r'].getSchema()
            schema = forced['r'].getSchema()
    
            # Get keys from the measurement catalog
            # The following is not strictly necessary as one could use the get("key_name")
            # method to access values in the
            # catalogs, but it is much more efficient to use get(key)
            fluxKey = schema["modelfit_CModel_flux"].asKey()
            fluxSigmaKey = schema["modelfit_CModel_fluxSigma"].asKey()
            fluxFlagKey = schema["modelfit_CModel_flag"].asKey()
            extKey = schema_m["base_ClassificationExtendedness_value"].asKey()
            extFlagKey = schema_m["base_ClassificationExtendedness_flag"].asKey()
            
            e1Key = schema_m["ext_shapeHSM_HsmShapeRegauss_e1"].asKey()
            e2Key = schema_m["ext_shapeHSM_HsmShapeRegauss_e2"].asKey()
            regaussFlagKey = schema_m["ext_shapeHSM_HsmShapeRegauss_flag"].asKey()
            resKey = schema_m["ext_shapeHSM_HsmShapeRegauss_resolution"].asKey()
    
            raKey = schema["coord_ra"].asKey()
            decKey = schema["coord_dec"].asKey()
    
            # Get cluster center position in pixel coordinates
            wcs = calexp['i'].getWcs()
            xClust, yClust = wcs.skyToPixel(afwGeom.Angle(raClust), afwGeom.Angle(deClust))
    
        # Loop over deblended sources in the filters forcedPhotCoadd catalogs
        rejected = {'star': 0, 'flag_flux':0, 'pos_flux': 0, 'stn': 0, 'gauss': 0}
        for i in range(len(forced[filters[0]])):
    
            # Select galaxies (and reject stars)
            if meas['r'][i].get(extFlagKey) or meas['r'][i].get(extKey) < 0.5:
                rejected['star'] += 1
                continue
    
            # Select sources which have a proper flux value in r, g and i bands
            # Notice that it would not be strictly necessary with forced photometry
            if N.any([forced[f][i].get(fluxFlagKey) for f in filters]):
                rejected['flag_flux'] += 1
                continue
    
            # Check the flux value, which must be > 0
            fluxes = {f: forced[f][i].get(fluxKey) for f in filters}
            fluxes_sigma = {f: forced[f][i].get(fluxSigmaKey) for f in filters}
            if any([fluxes[f] <= 0. for f in fluxes]):
                rejected['pos_flux'] += 1
                continue
    
            # Check the signal to noise (stn) value, which must be > 10
            stns = [forced[f][i].get(fluxKey)/forced[f][i].get(fluxSigmaKey) for f in filters
                    if forced[f][i].get(fluxSigmaKey) != 0]
            if any([stn < 10. for stn in stns]):
                rejected['stn'] += 1
                continue
    
            # Gauss regulerarization flag?
            if meas['r'][i].get(regaussFlagKey) or meas['r'][i].get(regaussFlagKey):
                rejected['gauss'] += 1
                continue
    
            # Get filter dependent values
            for f in filters:
    
                # Need to use a calibobject in order to convert flux to magnitude
                m, sm = calib[f].getMagnitude(fluxes[f], fluxes_sigma[f])
                mags[f].append(m)
                mags_sigma[f].append(sm)
                
                # Get ellipticities
                ell[f]['e1'].append(meas[f][i].get(e1Key))
                ell[f]['e2'].append(meas[f][i].get(e2Key))
    
                # Get resolution
                resolution[f].append(meas[f][i].get(resKey))
    
            ra = forced[filters[0]][i].get(raKey)
            dec = forced[filters[0]][i].get(decKey)
            x, y = wcs.skyToPixel(ra, dec)
            coords['ra'].append(float(ra))
            coords['dec'].append(float(dec))
            coords['id'].append(int(meas[filters[0]][i].get("id")))
            xSrc.append(x)
            ySrc.append(y)
    
        print "          Total # of sources:", len(forced[filters[0]])
        print "          Selected sources:", len(mags['r'])
        print "          Rejected sources:", len(forced[filters[0]])-len(mags['r'])
        for r in rejected:
            print "            %s: %i" % (r, rejected[r])

    # Make them all numpy array
    d = map(data.from_list_to_array, [mags, mags_sigma, ell, resolution, xSrc, ySrc, coords])

    # And dump them into a file
    cPickle.dump(d, open(args.output, 'w'))
    print "INFO: Data saved in", args.output 
