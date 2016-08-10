"""Main entry points for scripts."""

import os
import cPickle
import numpy as N
from argparse import ArgumentParser

import lsst.afw.geom as afwGeom

from . import data as D
from . import extinction as E

def load_data(argv=None):
    """
    Load data from the DM stack butler
    """

    description = """Load data from the DM stack butler."""
    prog = "clusters_data"
    usage = """usage: %s [options] datadir""" % prog

    parser = ArgumentParser(prog=prog, usage=usage, description=description)
    parser.add_argument('config', help='Configuration (yaml) file')
    parser.add_argument("--output",
                        help="Name of the output file (pkl file)")
    parser.add_argument("--overwrite", action="store_true", default=False,
                        help="Overwrite the output files if they exist already")
    parser.add_argument("--overwrite", action="store_true", default=False,
                        help="Overwrite the output files if they exist already")
    args = parser.parse_args(argv)
    
    config = D.load_config(args.config)

    if args.output is None:
        output = os.path.basename(args.config).replace('.yaml', '_data.hdf5')
        output_filtered = os.path.basename(args.config).replace('.yaml', '_filtered_data.hdf5')
    else:
        output = args.output
        output_filtered = "filtered_" + args.output
    
    print "INFO: Working on cluster %s (z=%.4f)" % (config['cluster'],
                                                    config['redshift'])
    print "INFO: Working on filters", config['filters']
    print "INFO: Butler located under %s" % config['butler']

    d = D.get_all_data(config['butler'], config['patches'],
                          config['filters'], add_extra=True)
    df = D.filter_table(d)
    D.write_data(d, output, overwrite=args.overwrite)
    D.write_data(df, output_filtered, overwrite=args.overwrite)

def extinction(argv=None):

    parser = ArgumentParser()
    parser.add_argument('config', help='Configuration (yaml) file')
    parser.add_argument('input', help='Input data file (output of clusters_data.py, i.e, a hdf5 file)')
    parser.add_argument("--output",
                        help="Name of the output file (pkl file)")
    parser.add_argument("--plot", action='store_true', default=False,
                        help="Make some plots")
    args = parser.parse_args(argv)

    config = D.load_config(args.config)
    if args.output is None:
        args.output = os.path.basename(args.config).replace('.yaml', '_extcorr.pkl')
    
    print "INFO: Working on cluster %s (z=%.4f)" % (config['cluster'], config['redshift'])
    print "INFO: Working on filters", config['filters']
    
    # Load the data
    d = D.read_data(args.input)
    df = d['forced']
    # Get the data we need
    mags = {f: df['modelfit_CModel_mag'][df['filter'] == f] for f in config['filters']}
    mags_sigma = {f: df['modelfit_CModel_magSigma'][df['filter'] == f] for f in config['filters']}
    coords = {'ra': df['coord_ra'][df['filter'] == config['filters'][0]],
              'dec': df['coord_dec'][df['filter'] == config['filters'][0]],
              'id': df['objectId'][df['filter'] == config['filters'][0]]}
    
    allra = [afwGeom.radToDeg(ra) for ra in coords['ra']]
    alldec = [afwGeom.radToDeg(dec) for dec in coords['dec']]
    ebv_sfd = N.array(E.query(allra, alldec, coordsys='equ', mode='sfd')['EBV_SFD'])
    albd = E.from_ebv_sfd_TO_megacam_albd(ebv_sfd)
    for f in mags:
        if f == 'i':
            mags[f] -= albd['i_new']
        else:
            mags[f] -= albd[f]


    data = [mags, mags_sigma, coords]
    cPickle.dump(data, open(args.output, 'w'))
    print "INFO: Milky Way dust extinction correctino applied"
    print "INFO: Data saved in", args.output

    if args.plot:
        print "INFO: Making some plots"
        E.plots(allra, alldec, ebv_sfd, albd, filters=['u', 'g', 'r', 'i_old', 'i_new', 'z'],
                title='Dust extinction map, %s, %i sources' % (config['cluster'], len(allra)),
                figname=config['cluster'])
