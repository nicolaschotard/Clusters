#!/usr/bin/env python

import os, yaml, cPickle
import numpy as N
from argparse import ArgumentParser

from Clusters import extinction

import lsst.afw.geom as afwGeom

if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument('config', help='Configuration (yaml) file')
    parser.add_argument('input', help='Input data file (output of clusters_data.py)')
    parser.add_argument("--output",
                        help="Name of the output file (pkl file)")
    parser.add_argument("--plot", action='store_true', default=False,
                        help="Make some plots")
    args = parser.parse_args()

    config = yaml.load(open(args.config))
    if args.output is None:
        args.output = os.path.basename(args.config).replace('.yaml', '_extcorr.pkl')

    filters = config['filters']
    
    raClust = afwGeom.degToRad(config['ra'])
    deClust = afwGeom.degToRad(config['dec'])

    print "INFO: Working on cluster %s (z=%.4f)" % (config['cluster'],
                                                    config['redshift'])
    print "INFO: Working on filters", filters
    
    # And dump them into a file
    data = cPickle.load(open(args.input, 'r'))
    mags, mags_sigma, ell, resolution, xSrc, ySrc, coords = data
    
    allra = [afwGeom.radToDeg(ra) for ra in coords['ra']]
    alldec = [afwGeom.radToDeg(dec) for dec in coords['dec']]
    ebv_sfd = N.array(extinction.query(allra, alldec, coordsys='equ', mode='sfd')['EBV_SFD'])
    albd = extinction.from_ebv_sfd_TO_megacam_albd(ebv_sfd)
    for f in mags:
        if f == 'i':
            mags[f] -= albd['i_new']
        else:
            mags[f] -= albd[f]

    data = [mags, mags_sigma, ell, resolution, xSrc, ySrc, coords]
    cPickle.dump(data, open(args.output, 'w'))
    print "INFO: Milky Way dust extinction correctino applied"
    print "INFO: Data saved in", args.output

    if args.plot:
        print "INFO: Making some plots"
        extinction.plots(allra, alldec, ebv_sfd, albd, filters=['u', 'g', 'r', 'i_old', 'i_new', 'z'],
                         title='Dust extinction map, %s, %i sources' % (config['cluster'], len(allra)),
                         figname=config['cluster'])
