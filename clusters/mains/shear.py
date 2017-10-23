"""Main entry points for scripts."""


from __future__ import print_function
import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from .. import utils as cutils
from .. import shear as cshear


def shear(argv=None):
    """Compute the shear."""
    description = """Compute the shear."""
    prog = "clusters_shear.py"
    usage = """%s [options] config input""" % prog

    parser = ArgumentParser(prog=prog, usage=usage, description=description,
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('config', help='Configuration (yaml) file')
    parser.add_argument(
        'input', help='Input data file: output of clusters_data.py, i.e, hdf5 file')
    parser.add_argument("--output",
                        help="Name of the output file (hdf5 file)")
    parser.add_argument("--overwrite", action="store_true", default=False,
                        help="Overwrite the output files if they exist already")
    parser.add_argument("--plot", action='store_true', default=False,
                        help="Make some plots")
    parser.add_argument("--step", type=int,
                        help="Step to use while making the kappa maps",
                        default=200)
    args = parser.parse_args(argv)

    config = cutils.load_config(args.config)
    if args.output is None:
        args.output = os.path.basename(
            args.input).replace('.hdf5', '_shear.hdf5')
        if not args.overwrite and os.path.exists(args.output):
            raise IOError(
                "Output already exists. Remove them or use --overwrite.")

    print("INFO: Working on cluster %s (z=%.4f)" %
          (config['cluster'], config['redshift']))
    print("INFO: Working on filters", config['filter'])

    # Load the data
    data = cutils.read_hdf5(args.input)
    meas = data['deepCoadd_meas']
    # converts astropy Table to the right wcs format
    wcs = cutils.load_wcs(data['wcs'])
    xclust, yclust = cutils.skycoord_to_pixel(
        [config['ra'], config['dec']], wcs)
    cshear.analysis(meas, float(xclust), float(yclust),
                    config=config, datafile=args.input, step=args.step)
