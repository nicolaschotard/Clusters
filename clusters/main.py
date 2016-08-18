"""Main entry points for scripts."""

import os
from argparse import ArgumentParser

from astropy.table import Table, hstack

import lsst.afw.geom as afwGeom

from . import data as D
from . import extinction as E


def load_data(argv=None):
    """Load data from the DM stack butler."""
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
    """Get color excess E(B-V) and store it in the data table for further use."""
    parser = ArgumentParser()
    parser.add_argument('config', help='Configuration (yaml) file')
    parser.add_argument('input', help='Input data file: output of clusters_data.py, i.e, hdf5 file')
    parser.add_argument("--output",
                        help="Name of the output file (pkl file)")
    parser.add_argument("--overwrite", action="store_true", default=False,
                        help="Overwrite the output files if they exist already")
    parser.add_argument("--plot", action='store_true', default=False,
                        help="Make some plots")
    args = parser.parse_args(argv)

    config = D.load_config(args.config)
    if args.output is None:
        args.output = os.path.basename(args.input).replace('.hdf5', '_extinction.hdf5')

    print "INFO: Working on cluster %s (z=%.4f)" % (config['cluster'], config['redshift'])
    print "INFO: Working on filters", config['filters']

    # Load the data
    d = D.read_data(args.input)['forced']

    # Get the coordinates
    ras = [afwGeom.radToDeg(ra) for ra in d['coord_ra']]
    decs = [afwGeom.radToDeg(dec) for dec in d['coord_dec']]

    # Query for E(b-v) and compute the extinction
    ebmv = {'ebv_sfd': E.query(ras, decs, coordsys='equ', mode='sfd')['EBV_SFD']}
    albds = {}
    for k in ebmv:
        albd = E.from_ebv_sfd_TO_megacam_albd(ebmv[k])
        albds.update({k.replace('ebv_', 'albd_%s_' % f): albd[f] for f in albd})

    # Create a new table and save it
    new_tab = hstack([d['objectId', 'coord_ra', 'coord_dec', 'filter'],
                      Table(ebmv), Table(albds)], join_type='inner')
    new_tab.write(args.output, path='extinction', compression=True,
                  serialize_meta=True, overwrite=args.overwrite)
    print "INFO: Milky Way dust extinction correctino applied"
    print "INFO: Data saved in", args.output

    # Make some plots if asked
    if args.plot:
        print "INFO: Making some plots"
        filt = new_tab['filter'] == config['filters'][0]
        E.plots(new_tab['coord_ra'][filt],
                new_tab['coord_dec'][filt],
                new_tab['ebv_sfd'], albds['albd_sfd'][filt],
                filters=['u', 'g', 'r', 'i_old', 'i_new', 'z'],
                title='Dust extinction map, %s, %i sources' % (config['cluster'],
                                                               len(new_tab['coord_ra'][filt])),
                figname=config['cluster'])
