"""Main entry points for scripts."""

from astropy.table import Table, hstack
import yaml
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from .. import background
from .. import data as cdata

def getbackground(argv=None):
    """Get a cluster background galaxies."""
    description = """Get a cluster background galaxies."""
    prog = "clusters_getbackground.py"
    usage = """%s [options] config input""" % prog

    parser = ArgumentParser(prog=prog, usage=usage, description=description,
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('config', help="Configuration (yaml) file")
    parser.add_argument('data', help="Catalogue file data.hdf5"
                        " (*_hdf5 output of clusters_data)")
    parser.add_argument("--zdata",
                        help="Photometric redshift data file (including pdz information) if not "
                             "in data.hdf5")
    parser.add_argument("--output",
                        help="Filename for the shear catalogue with bkg flags to also be stored "
                        "in a seperate file. If not specified, flags are added to data.hdf5")
    parser.add_argument("--zmin", type=float,
                        help="Minimum redshift for photoz hard cut")
    parser.add_argument("--zmax", type=float,
                        help="Maximum redshift for photoz hard cut")
    parser.add_argument("--thresh_prob", type=float,
                        help="Threshod redshift probability to select galaxy (in percent)")
    parser.add_argument("--plot", default=False, action='store_true', help="Make some plots")
    parser.add_argument("--overwrite", default=False, action='store_true',
                        help="Overwrite the paths in the output file if they exist already")
    parser.add_argument("--rs", default=False, action='store_true',
                        help="Also slect galaxy based on red sequence")

    args = parser.parse_args(argv)

    config = yaml.load(open(args.config))

    if args.zmin is None:
        args.zmin = config['redshift'] + 0.1
    if args.zmax is None:
        args.zmax = 1.25
    if args.thresh_prob is None:
        args.thresh_prob = 1.
    if args.output is None:
        args.output = args.data
    if args.zdata is None:
        args.zdata = args.data

    data = cdata.read_hdf5(args.data)
    zdata = cdata.read_hdf5(args.zdata)

    # If the user did not define a configuration to run the photoz,
    # add default one to the config dictionary
    if not 'zphot' in config:
        config['zphot']={'zphot_ref':{}}

    # Loop over all zphot configurations found in config.yaml file
    for k in config['zphot'].keys():
        z_config = config['zphot'][k]
        z_flag1, z_flag2 = background.get_zphot_background(config, zdata[k],
                                                           zmin=args.zmin,
                                                           zmax=args.zmax,
                                                           z_config=z_config,
                                                           thresh=args.thresh_prob,
                                                           plot=args.plot)                                                
        new_tab = hstack([Table([zdata[k]['id' if 'id' in data[k].keys() else 'objectId']], names=['id' if 'id' in data[k].keys() else 'objectId']),
                          Table([z_flag1], names=['flag_z_hard']),
                          Table([z_flag2], names=['flag_z_pdz'])],
                         join_type='inner')

        cdata.overwrite_or_append(args.output, 'flag_' + k, new_tab, overwrite=args.overwrite)

    if args.rs:
        rs_flag = background.get_rs_background(config, data['deepCoadd_forced_src'])
        new_tab = hstack([Table([data['deepCoadd_forced_src']['id' if 'id' in data.keys() else 'objectId']], names=['id' if 'id' in data.keys() else 'objectId']),
                          Table([rs_flag], names=['flag_rs'])],
                         join_type='inner')
        cdata.overwrite_or_append(args.output, 'flag_rs', Table([rs_flag]))
