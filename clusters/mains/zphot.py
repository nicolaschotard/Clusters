"""Main entry points for scripts."""

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from .. import data as cdata
from .. import zphot as czphot


def photometric_redshift(argv=None):
    """Comput photometric redshift using LEPHARE."""
    parser = ArgumentParser(prog="clusters_zphot.py",
                            usage="clusters_zphot.py [options] config input",
                            description="Comput photometric redshift using LEPHARE.",
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('config', help='Configuration (yaml) file')
    parser.add_argument('input', help='Input data file')
    parser.add_argument("--output",
                        help="Name of the output file (hdf5 file)")
    parser.add_argument("--extinction",
                        help="Output of clusters_extinction (hdf5 file)."
                        "Use to compute the extinction-corrected magnitudes.")
    parser.add_argument("--plot", action='store_true', default=False,
                        help="Make some plots")
    parser.add_argument("--zrange", default="0,999",
                        help="Redshift range used to plot the map (min,max)")
    parser.add_argument("--mag", type=str, default='modelfit_CModel_mag',
                        help="Magnitude name [default]")
    parser.add_argument("--data",
                        help="Photoz output file, used for the analysis only (plots)")
    parser.add_argument("--overwrite", action="store_true", default=False,
                        help="Overwrite the paths in the output file if they exist already")
    args = parser.parse_args(argv)

    config = cdata.load_config(args.config)
    if args.data is not None:
        doplot(czphot.LEPHARO(args.data, args.data.replace('_zphot', '')),
               config, float(args.zrange.split(',')[0]), float(args.zrange.split(',')[1]))
        sys.exit()

    if args.output is None:  # if no output name specified, append table to input file
        args.output = args.input

    print "INFO: Working on cluster %s (z=%.4f)" % (config['cluster'], config['redshift'])

    # Load the data
    print "INFO: Loading the data from", args.input
    data = cdata.read_hdf5(args.input)['deepCoadd_forced_src']

    # Compute extinction-corrected magitudes
    if args.extinction is not None:
        print "INFO: Computing extinction-corrected magnitude for", args.mag
        edata = cdata.read_hdf5(args.extinction)['extinction']
        cdata.correct_for_extinction(data, edata, mag=args.mag)
        args.mag += "_extcorr"

    # Make sure the selected magnitude does exist in the data table
    if args.mag not in data.keys():
        raise IOError("%s is not a column of the input table" % args.mag)

    # If the user did not define a configuration to run the photoz,
    # add default one to the config dictionary
    if not 'zphot' in config:
        config['zphot']={'zphot_ref':{}} 
    
    # Loop over all zphot configurations present in the config.yaml file
    for zconfig in config['zphot'].keys():
        zcode = config['zphot'][zconfig]['code'] if 'code' in config['zphot'][zconfig].keys() \
                                                 else 'lephare'

        # If a spectroscopic sample is provided, LEPHARE/BPZ will run using the adaptative method
        # (zero points determination); still to be implemented for BPZ...
   
        zpara = config['zphot'][zconfig]['zpara'] if 'zpara' in config['zphot'][zconfig] else None
        spectro_file = config['zphot'][zconfig]['zspectro_file'] if 'zspectro_file' \
            in config['zphot'][zconfig] else None
        kwargs = {'basename': config['cluster'],
                  'filters': [f for f in config['filter'] if f in set(data['filter'].tolist())],
                  'ra': data['coord_ra_deg'][data['filter'] == config['filter'][0]],
                  'dec': data['coord_dec_deg'][data['filter'] == config['filter'][0]],
                  'id': data['objectId'][data['filter'] == config['filter'][0]]}
        path = zconfig
        print "INFO: Running", zcode, "using configuration from", zpara, spectro_file

        if zcode == 'bpz':  # Run BPZ
            zphot = czphot.BPZ([data[args.mag][data['filter'] == f] for f in kwargs['filters']],
                               [data[args.mag.replace("_extcorr", "") + "Sigma"][data['filter'] == f]
                                for f in kwargs['filters']],
                               zpara=zpara, spectro_file=spectro_file, **kwargs)

        if zcode == 'lephare': # Run LEPHARE
            zphot = czphot.LEPHARE([data[args.mag][data['filter'] == f] for f in kwargs['filters']],
                                   [data[args.mag.replace("_extcorr", "") + 
                                         "Sigma"][data['filter'] == f] for f in kwargs['filters']],
                                   zpara=zpara, spectro_file=spectro_file, **kwargs)
            zphot.check_config()

        zphot.run()
        zphot.data_out.save_zphot(args.output, path, overwrite=args.overwrite)

    # Plot
    if args.plot:
        doplot(zphot.data_out, config,
               float(args.zrange.split(',')[0]),
               float(args.zrange.split(',')[1]))

    if args.plot:
        czphot.P.show()
