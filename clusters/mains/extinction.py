"""Main entry points for scripts."""

from astropy.table import Table, hstack
from extinctions import reddening
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from .. import data as cdata
from .. import extinction as cextinction


def extinction(argv=None):
    """Get color excess E(B-V) and store it in the data table for further use."""
    description = """Get color excess E(B-V) and store it in the data table for further use."""
    prog = "clusters_extinction.py"
    usage = """%s [options] config input""" % prog

    parser = ArgumentParser(prog=prog, usage=usage, description=description,
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('config', help='Configuration (yaml) file')
    parser.add_argument('input', help='Input data file: output of clusters_data.py, i.e, hdf5 file')
    parser.add_argument("--output",
                        help="Name of the output file (hdf5 file)")
    parser.add_argument("--overwrite", action="store_true", default=False,
                        help="Overwrite the paths in the output file if they exist already")
    parser.add_argument("--plot", action='store_true', default=False,
                        help="Make some plots")
    args = parser.parse_args(argv)

    config = cdata.load_config(args.config)
    if args.output is None:
        args.output = args.input

    print "INFO: Working on cluster %s (z=%.4f)" % (config['cluster'], config['redshift'])
    print "INFO: Working on filters", config['filter']

    # Load the data
    data = cdata.read_hdf5(args.input, path='deepCoadd_meas', dic=False)

    # Query for E(b-v) and compute the extinction
    red = reddening.Reddening(data['coord_ra_deg'].tolist(), data['coord_dec_deg'].tolist())
    ebmv = {'ebv_sfd': red.from_argonaut()}

    albds = {}
    for k in ebmv:
        albd = cextinction.from_ebv_sfd_to_megacam_albd(ebmv[k])
        albds.update({k.replace('ebv_', 'albd_%s_' % f): albd[f] for f in albd})

    # Create a new table and save it
    new_tab = hstack([data['id', 'coord_ra', 'coord_dec', 'filter'],
                      Table(ebmv), Table(albds)], join_type='inner')
#    new_tab.write(args.output, path='extinction', compression=True,
#                  serialize_meta=True, overwrite=args.overwrite)

    cdata.overwrite_or_append(args.output, 'extinction', new_tab, overwrite=args.overwrite)

    print "INFO: Milky Way dust extinction correction applied"
    print "INFO: Data saved in", args.output

    # Make some plots if asked
    if args.plot:
        print "INFO: Making some plots"
        filt = new_tab['filter'] == config['filter'][0]
        cextinction.plots(new_tab['coord_ra'][filt],
                          new_tab['coord_dec'][filt],
                          new_tab['ebv_sfd'], albds['albd_g_sfd'][filt],
                          title='Dust extinction map, %s, %i sources, g_sfd' %
                          (config['cluster'], len(new_tab['coord_ra'][filt])),
                          figname=config['cluster'])
