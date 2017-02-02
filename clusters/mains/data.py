"""Main entry points for scripts."""

import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from .. import data as cdata


def load_data(argv=None):
    """Load data from the DM stack butler."""
    description = """Load data from the DM stack butler."""
    prog = "clusters_data.py"
    usage = """%s [options] config""" % prog

    parser = ArgumentParser(prog=prog, usage=usage, description=description,
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('config', help='Configuration (yaml) file')
    parser.add_argument("--output",
                        help="Name of the output file (hdf5 file)")
    parser.add_argument("--catalogs",
                        default='forced_src,deepCoadd_meas,deepCoadd_forced_src',
                        help="List of catalogs to load (coma separated)")
    parser.add_argument("--overwrite", action="store_true", default=False,
                        help="Overwrite the output files if they exist already")
    parser.add_argument("--show", action="store_true", default=False,
                        help="Show and save the list of available keys in the catalogs, and exit.")
    args = parser.parse_args(argv)

    config = cdata.load_config(args.config)

    if args.output is None:
        output = os.path.basename(args.config).replace('.yaml', '_data.hdf5')
        output_filtered = os.path.basename(args.config).replace('.yaml', '_filtered_data.hdf5')
    else:
        output = args.output if args.output.endswith('.hdf5') else args.output + ".hdf5"
        output_filtered = output.replace('.hdf5', '_filtered_data.hdf5')

    if not args.overwrite and (os.path.exists(output) or os.path.exists(output_filtered)):
        raise IOError("Output(s) already exist(s). Remove them or use overwrite=True.")

    print "\nINFO: Working on:"
    print "  - cluster %s (z=%.4f)" % (config['cluster'], config['redshift'])
    print "  - filters", config['filter']
    print "  - patches", config['patch']
    print "INFO: Butler located under %s" % config['butler']

    data = cdata.Catalogs(config['butler'])

    if args.show:
        data.show_keys(args.catalogs.split(','))
        return
    config['output_name'] = output
    config['overwrite'] = args.overwrite
    data.load_catalogs(args.catalogs.split(','), matchid=True, **config)
    print "\nINFO: Applying filters on the data to keep a clean sample of galaxies"
    catalogs = cdata.read_hdf5(output)
    data = cdata.Catalogs(config['butler'])
    data.catalogs = cdata.filter_table(catalogs)
    data.save_catalogs(output_filtered, overwrite=args.overwrite, delete_catalog=True)
