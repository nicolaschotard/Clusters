"""Main entry points for scripts."""


from __future__ import print_function
import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from .. import data as cdata
from pzmassfitter import dmstackdriver


def mass(argv=None):
    """Compute cluster mass"""
    description = """Compute the mass."""
    prog = "clusters_mass.py"
    usage = """%s [options] config input""" % prog

    parser = ArgumentParser(prog=prog, usage=usage, description=description,
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('config', help='Configuration (yaml) file')
    parser.add_argument('input', help='Input data file: output of clusters_data.py, i.e, hdf5 file')
    parser.add_argument("--output",
                        help="Name of the output file (hdf5 file)")
    parser.add_argument("--overwrite", action="store_true", default=False,
                        help="Overwrite the output files if they exist already")
    parser.add_argument("--plot", action='store_true', default=False,
                        help="Make some plots")
    parser.add_argument("--nsamples", default=10000, type=int,
                        help="Number of sample to run")
    parser.add_argument("--testing", action="store_true", default=False,
                        help="Simplify model for testing purposes")
    args = parser.parse_args(argv)

    config = cdata.load_config(args.config)

    # Select the zphot configuration to use for mass estimation
    # If not specified in yaml file, order the zphot configuration names alphabetically
    # and take the first one.
    
    mconfig = config['mass'] if 'mass' in config else {'zconfig':'zphot_ref'}
    tag ='' if 'zflagconfig' not in mconfig else '_'+mconfig['zflagconfig']
    
    print("Cluster mass computed using ", mconfig,
          " configuration for photoz estimation and backgroud selection")


    print("INFO: Working on cluster %s (z=%.4f)" % (config['cluster'], config['redshift']))
    print("INFO: Working on filters", config['filter'])

    # Load the data
    data = cdata.read_hdf5(args.input)

    cluster = config['cluster']
    zcluster = config['redshift']
    cluster_ra = config['ra']
    cluster_dec = config['dec']

    # WTG quantities for shear calibration
    wtg_shearcal = False if 'wtg_shearcal' not in mconfig else config['mass']['wtg_shearcal']
    psfsize = None if 'psfsize' not in mconfig else config['mass']['psfsize']

    # Choose lin or log sampling for the mass
    mprior = 'lin' if 'mprior' not in mconfig else config['mass']['mprior']
    if mprior == 'lin':
        logprior = False
    else:
        logprior = True

    if args.output is None:
        args.output = args.input.replace('.hdf5', '_mass'+mprior+'_cal'+str(wtg_shearcal)+'_'+mconfig['zconfig']+tag+'.hdf5')
        if not args.overwrite and os.path.exists(args.output):
            raise IOError("Output already exists. Remove them or use --overwrite.")

    ###let's assume that all quality cuts were made previously

    if args.testing:
        print('TESTING!!!!')
        masscontroller = dmstackdriver.makeTestingController()
        options, cmdargs = masscontroller.modelbuilder.createOptions(concentration=4.)
        options, cmdargs = masscontroller.runmethod.createOptions(outputFile=args.output,
                                                                  options=options,
                                                                  args=cmdargs)

    else:
        masscontroller = dmstackdriver.controller
        options, cmdargs = masscontroller.modelbuilder.createOptions()
        options, cmdargs = masscontroller.runmethod.createOptions(outputFile=args.output,
                                                                  nsamples=args.nsamples,
                                                                  burn=2000,
                                                                  options=options,
                                                                  args=cmdargs)

    options, cmdargs = masscontroller.filehandler.createOptions(cluster=cluster,
                                                                zcluster=zcluster,
                                                                cat=data,
                                                                mconfig=mconfig, 
                                                                cluster_ra=cluster_ra,
                                                                cluster_dec=cluster_dec,
                                                                wtg_shearcal=wtg_shearcal,
                                                                psfsize=psfsize,
                                                                logprior=logprior,
                                                                options=options,
                                                                args=cmdargs)

    masscontroller.load(options, args)
    masscontroller.run()
    masscontroller.dump()
    masscontroller.finalize()
