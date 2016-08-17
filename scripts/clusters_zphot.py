#!/usr/bin/env python

import os
import sys
import yaml
import cPickle
from argparse import ArgumentParser

from Clusters import zphot

def doplot(data):
    """Make a few plots."""
    print "INFO: Making some plots"
    data.hist('Z_BEST', min=0, nbins=100, xlabel='Photometric redshift',
              figname=config['cluster'],
              title="LEPHARE photo-z for %s (%i sources)" % \
              (config['cluster'], data.nsources), zclust=config['redshift'])
    data.hist('CHI_BEST', nbins=100, max=100, figname=config['cluster'],
              title="LEPHARE photo-z for %s (%i sources)" % \
              (config['cluster'], data.nsources))
    data.plot('CHI_BEST', 'Z_BEST', miny=0, figname=config['cluster'])
    data.plot_map(title="LEPHARE photometric redshift map for %s (%i sources)"%\
                  (config['cluster'], data.nsources), figname=config['cluster'],
                  zmin=args.zmin, zmax=args.zmax)
    zphot.P.show()

if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument('config', help='Configuration (yaml) file')
    parser.add_argument('input', help='Input data file')
    parser.add_argument("--output",
                        help="Name of the output file (pkl file)")
    parser.add_argument("--data",
                        help="LEPHARE output file, used for the plots")
    parser.add_argument("--zpara",
                        help="LEPHARE configuration file (zphot.para)")
    parser.add_argument("--plot", action='store_true', default=False,
                        help="Make some plots")
    parser.add_argument("--zmin", type=float, default=0,
                        help="Redshift cut used to plot the map (min value)")
    parser.add_argument("--zmax", type=float, default=999,
                        help="Redshift cut used to plot the map (max value)")
    args = parser.parse_args()

    config = yaml.load(open(args.config))
    if args.output is None:
        args.output = os.path.basename(args.config).replace('.yaml', '_lephare_output.pkl')

    filters = config['filters']

    print "INFO: Working on cluster %s (z=%.4f)" % (config['cluster'],
                                                    config['redshift'])
    print "INFO: Working on filters", filters

    if args.data is not None:
        doplot(zphot.LEPHARO(args.data, args.data.replace('out', 'all')))
        sys.exit()

    # And dump them into a file
    mags, mags_sigma, ell, resolution, xSrc, ySrc, coords = cPickle.load(open(args.input, 'r'))

    zp = zphot.LEPHARE(zphot.dict_to_array(mags, filters=filters),
                       zphot.dict_to_array(mags_sigma, filters=filters),
                       config['cluster'], filters=filters, zpara=args.zpara,
                       RA=coords['ra'], DEC=coords['dec'], ID=coords['id'])
    zp.check_config()
    zp.run()

    cPickle.dump(zp.data_out.data_dict, open(args.output, 'w'))
    print "INFO: LEPHARE data saved in", args.output

    if args.plot:
        doplot(zp.data_out)
