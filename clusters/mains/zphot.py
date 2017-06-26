"""Main entry points for scripts."""

import sys, pdb
import numpy as N
from astropy.table import Table
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
    parser.add_argument("--extinction", default=False, action='store_true',
                        help="Compute the extinction-corrected magnitude usings "
                        "clusters_extinction outputs.")
    parser.add_argument("--dustmap", default='sfd',
                        help="Dustmap name used to compute the extinction-corrected magnitudes")
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
    parser.add_argument("--zeropoints", type=str,
                        help="Input file with zero points for some or all filters"
                        ". A three columns file: # filt zp zp_err")
    args = parser.parse_args(argv)

    config = cdata.load_config(args.config)
    if args.data is not None:
        doplot(czphot.LEPHARO(args.data, args.data.replace('_zphot', '')),
               config, float(args.zrange.split(',')[0]), float(args.zrange.split(',')[1]))
        sys.exit()

    if args.output is None:  # if no output name specified, append table to input file
        args.output = args.input

    print "INFO: Working on cluster %s (z=%.4f)" % (config['cluster'], config['redshift'])

    if not 'sim' in config:
        config['sim'] = {'flag':False}

        # Load the data
        print "INFO: Loading the data from", args.input
        tables = cdata.read_hdf5(args.input)
        data = tables['deepCoadd_forced_src']

        # Compute extinction-corrected magitudes
        if args.extinction and 'extinction' in tables.keys():
            print "INFO: Computing extinction-corrected magnitude for", args.mag, \
              "using the '%s' dust map" % args.dustmap
            cdata.correct_for_extinction(data, tables['extinction'], mag=args.mag, ext=args.dustmap)
        args.mag += "_extcorr"

        # Make sure the selected magnitude does exist in the data table
        if args.mag not in data.keys():
            raise IOError("%s is not a column of the input table" % args.mag)
        
        # Apply zeropoints?
        if args.zeropoints is not None:
            print "INFO: Applying zeropoints for the follwoing filter:"
            fzpoint, zpoint, dzpoint = N.loadtxt(open(args.zeropoints), unpack=True, dtype='string')
            for fz, zp, dzp in zip(fzpoint, zpoint, dzpoint):
                if fz not in data['filter']:
                    print " - WARNING: %s not in the filter list" % fz
                else:
                    print " - correcting %s mags with zp = %.5f +/- %.5f" % (fz, float(zp), float(dzp))
                    data[args.mag][data['filter'] == fz] += float(zp)
                    merr = data[args.mag.replace("_extcorr", "") + "Sigma"][data['filter'] == fz]
                    new_err = N.sqrt(merr ** 2 + float(dzp) ** 2)
                    data[args.mag.replace("_extcorr", "") + "Sigma"][data['filter'] == fz] = new_err

    # If the user did not define a configuration to run the photoz,
    # add default one to the config dictionary
    if not 'zphot' in config:
        config['zphot'] = {'zphot_ref':{}}
    
    if not config['sim']['flag']: # we're not dealing with simulation data

        # Loop over all zphot configurations present in the config.yaml file
        for zconfig in config['zphot'].keys():
            zcode = config['zphot'][zconfig]['code'] if 'code' in config['zphot'][zconfig].keys() \
                                                     else 'lephare'

            # If a spectroscopic sample is provided, LEPHARE/BPZ will run using the adaptative
            # method (zero points determination); still to be implemented for BPZ...

            zpara = config['zphot'][zconfig]['zpara'] \
                    if 'zpara' in config['zphot'][zconfig] else None
            spectro_file = config['zphot'][zconfig]['zspectro_file'] if 'zspectro_file' \
                           in config['zphot'][zconfig] else None
            kwargs = {'basename': config['cluster'],
                      'filters': [f for f in config['filter'] if f in set(data['filter'].tolist())],
                      'ra': data['coord_ra_deg'][data['filter'] == config['filter'][0]],
                      'dec': data['coord_dec_deg'][data['filter'] == config['filter'][0]],
                      'id': data['id' if 'id' in data.keys() else 'objectId'][data['filter'] == \
                                                                              config['filter'][0]]}
            path = zconfig
            print "INFO: Running", zcode, "using configuration from", zpara, spectro_file

            if zcode == 'bpz':  # Run BPZ
                zphot = czphot.BPZ([data[args.mag][data['filter'] == f] for f in kwargs['filters']],
                                   [data[args.mag.replace("_extcorr", "") + \
                                         "Sigma"][data['filter'] == f]
                                    for f in kwargs['filters']],
                                   zpara=zpara, spectro_file=spectro_file, **kwargs)

            if zcode == 'lephare': # Run LEPHARE
                zphot = czphot.LEPHARE([data[args.mag][data['filter'] == f]
                                        for f in kwargs['filters']],
                                       [data[args.mag.replace("_extcorr", "") +
                                             "Sigma"][data['filter'] == f]
                                        for f in kwargs['filters']],
                                       zpara=zpara, spectro_file=spectro_file, **kwargs)
                zphot.check_config()

        zphot.run()
        zphot.data_out.save_zphot(args.output, path, overwrite=args.overwrite)

    else:
        # we are dealing with simulated data --> make fake p(z) from real z.
        # assumes the z information is in a 2 column (id, z) txt file,
        # where id corresponds to DM stack ObjectId
        path = 'zphot_ref'
        data_z_sim = N.loadtxt(config['sim']['zfile'],
                               dtype={'names':('id', 'z'), 'formats':('i8', 'f8')},
                               unpack=True, comments="#")
        min_pdz = 0
        max_pdz = 4
        pdz_step = 0.01
        zrange = N.arange(min_pdz, max_pdz, pdz_step)
        id_sim = data_z_sim[0]
        z_sim = data_z_sim[1]
        zbinsgrid = N.vstack(len(id_sim)*[zrange])
        pdz = [N.zeros(len(zrange))]
        for i, z in enumerate(z_sim):
            pdz_tmp = N.zeros(len(zrange))
            pdz_tmp[N.digitize(z, zrange)] = 1/pdz_step # normalised so that int_zmin^zmax pdz = 1
            pdz = pdz_tmp if i == 0 else N.vstack((pdz, pdz_tmp))
        #pdb.set_trace()
        pdz_values = Table([id_sim, z_sim, pdz, zbinsgrid],
                           names=('objectId', 'Z_BEST', 'pdz', 'zbins'))
        cdata.overwrite_or_append(args.output, path, pdz_values, overwrite=True)

    # Plot
    if args.plot:
        doplot(zphot.data_out, config,
               float(args.zrange.split(',')[0]),
               float(args.zrange.split(',')[1]))

    if args.plot:
        czphot.P.show()
