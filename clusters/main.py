"""Main entry points for scripts."""

import numpy as N
import os
import sys
from astropy.table import Table, Column, hstack
from extinctions import reddening
import yaml
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from . import data as cdata
from . import extinction as cextinction
from . import zphot as czphot
from . import shear as cshear
from . import background
from pzmassfitter import dmstackdriver

#import pdb


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
#        args.output = os.path.basename(args.input).replace('.hdf5', '_extinction.hdf5')
        args.output = args.input
#        if not args.overwrite and os.path.exists(args.output):
#            raise IOError("Output already exists. Remove them or use --overwrite.")

    print "INFO: Working on cluster %s (z=%.4f)" % (config['cluster'], config['redshift'])
    print "INFO: Working on filters", config['filter']

    # Load the data
    data = cdata.read_hdf5(args.input, path='deepCoadd_forced_src', dic=False)

    # Query for E(b-v) and compute the extinction
    red = reddening.Reddening(data['coord_ra_deg'].tolist(), data['coord_dec_deg'].tolist())
    ebmv = {'ebv_sfd': red.from_argonaut()}

    albds = {}
    for k in ebmv:
        albd = cextinction.from_ebv_sfd_to_megacam_albd(ebmv[k])
        albds.update({k.replace('ebv_', 'albd_%s_' % f): albd[f] for f in albd})

    # Create a new table and save it
    new_tab = hstack([data['objectId', 'coord_ra', 'coord_dec', 'filter'],
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


def doplot(data, config, zmin=0, zmax=999):
    """Make a few plots."""
    print "INFO: Making some plots"
    data.hist('Z_BEST', minv=0, nbins=100, xlabel='Photometric_redshift',
              title="LEPHARE photo-z for %s (%i sources)" %
              (config['cluster'], data.nsources), zclust=config['redshift'])
    data.hist('CHI_BEST', nbins=100, maxv=100,
              title="LEPHARE photo-z for %s (%i sources)" % (config['cluster'], data.nsources))
    data.plot('CHI_BEST', 'Z_BEST', miny=0)
    data.plot_map(title="LEPHARE photometric redshift map for %s (%i sources)" %
                  (config['cluster'], data.nsources), zmin=zmin, zmax=zmax)


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

    if args.output is None: # if no output name specified, append table to input file
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


    # Loop over all zphot codes present in the config.yaml file
    for zconfig in config['zphot'].keys():
        zcode = config['zphot'][zconfig]['code'] if 'code' in config['zphot'][zconfig].keys() else 'lephare'
        # If a spectroscopic sample is provided, LEPHARE/BPZ will run using the adaptative method
        # (zero points determination); still to be implemented for BPZ...
        # spectro_file = config['zphot'][zcode]['zspectro_file'] if 'zspectro_file' in config['zphot'][zcode] else None

        zpara = config['zphot'][zconfig]['zpara'] if 'zpara' in config['zphot'][zconfig] else None
        spectro_file = config['zphot'][zconfig]['zspectro_file'] if 'zspectro_file' in config['zphot'][zconfig] else None
        kwargs = {'basename': config['cluster'],
                    'filters': [f for f in config['filter'] if f in set(data['filter'].tolist())],
                    'ra': data['coord_ra_deg'][data['filter'] == config['filter'][0]],
                    'dec': data['coord_dec_deg'][data['filter'] == config['filter'][0]],
                    'id': data['objectId'][data['filter'] == config['filter'][0]]}
        path = zconfig
        print "INFO: Running", zcode, "using configuration from", zpara, spectro_file

        if zcode == 'bpz': # Run BPZ
            zphot = czphot.BPZ([data[args.mag][data['filter'] == f] for f in kwargs['filters']],
                               [data[args.mag.replace("_extcorr", "") + "Sigma"][data['filter'] == f]
                                for f in kwargs['filters']],
                                zpara=zpara, spectro_file=spectro_file, **kwargs)
               
        if zcode == 'lephare': # Run LEPHARE
            zphot = czphot.LEPHARE([data[args.mag][data['filter'] == f] for f in kwargs['filters']],
                                       [data[args.mag.replace("_extcorr", "") + \
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
        args.zmin = 0.
    if args.zmax is None:
        args.zmax = config['redshift'] + 0.1
    if args.thresh_prob is None:
        args.thresh_prob = 1.
    if args.output is None:
        args.output = args.data
    if args.zdata is None:
        args.zdata = args.data
        
    data = cdata.read_hdf5(args.data)
    zdata = cdata.read_hdf5(args.zdata)

    # Loop over all zphot configurations found in config.yaml file
    for k in config['zphot'].keys():
        z_config = config['zphot'][k]
        z_flag1, z_flag2 = background.get_zphot_background(config, zdata[k],
                                                         zmin=args.zmin,
                                                         zmax=args.zmax,
                                                         z_config=z_config,
                                                         thresh=args.thresh_prob,
                                                         plot=args.plot)
        new_tab=hstack([Table([zdata[k]['objectId']], names=['objectId']), Table([z_flag1], names=['flag_z_hard']),
                           Table([z_flag2], names=['flag_z_pdz'])],
                           join_type='inner')
            
        cdata.overwrite_or_append(args.output, 'flag_'+k, new_tab, overwrite=args.overwrite)
            
    if args.rs:
        rs_flag = background.get_rs_background(config, data['deepCoadd_forced_src'])
        new_tab=hstack([Table([zdata[k]['objectId']], names=['objectId']), Table([rs_flag], names=['flag_rs'])],
                        join_type='inner')
        cdata.overwrite_or_append(args.output, 'flag_rs', Table([rs_flag]))

def shear(argv=None):
    """Compute the shear."""
    description = """Compute the shear."""
    prog = "clusters_shear.py"
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
    args = parser.parse_args(argv)

    config = cdata.load_config(args.config)
    if args.output is None:
        args.output = os.path.basename(args.input).replace('.hdf5', '_shear.hdf5')
        if not args.overwrite and os.path.exists(args.output):
            raise IOError("Output already exists. Remove them or use --overwrite.")

    print "INFO: Working on cluster %s (z=%.4f)" % (config['cluster'], config['redshift'])
    print "INFO: Working on filters", config['filter']

    # Load the data
    data = cdata.read_hdf5(args.input)
    meas = data['deepCoadd_meas']
    wcs = data['wcs']
    xclust, yclust = cdata.skycoord_to_pixel([config['ra'], config['dec']], wcs)
    cshear.analysis(meas, float(xclust), float(yclust))



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
    if args.output is None:
        args.output = args.input.replace('.hdf5', '_mass.hdf5')
        if not args.overwrite and os.path.exists(args.output):
            raise IOError("Output already exists. Remove them or use --overwrite.")


    print "INFO: Working on cluster %s (z=%.4f)" % (config['cluster'], config['redshift'])
    print "INFO: Working on filters", config['filter']

    # Load the data
    data = cdata.read_hdf5(args.input)
#    meas = data['deepCoadd_meas']

    cluster = config['cluster']
    zcluster = config['redshift']
    cluster_ra = config['ra']
    cluster_dec = config['dec']

    ###let's assume that all quality cuts were made previously

    if args.testing:
        print 'TESTING!!!!'
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
                                                                zconfig=config['zphot'].keys()[0], # Run only on the first zphot configuration
                                                                cluster_ra=cluster_ra,
                                                                cluster_dec=cluster_dec,
                                                                options=options,
                                                                args=cmdargs)

    masscontroller.load(options, args)
    masscontroller.run()
    masscontroller.dump()
    masscontroller.finalize()



def pipeline(argv=None):
    """Pipeline for a standard cluster analysis."""
    description = """Run standard comand lines for full cluster analysis."""
    prog = "clusters_pipeline.py"
    usage = """%s [options] config""" % prog

    parser = ArgumentParser(prog=prog, usage=usage, description=description,
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('config', help='Configuration (yaml) file')
    parser.add_argument("--output",
                        help="Name of the output file (hdf5 file)")
    parser.add_argument("--overwrite", action="store_true", default=False,
                        help="Overwrite the output files if they exist already")
    args = parser.parse_args(argv)

    load_data()
    extinction(argv=None)
    print "TBD", args.config
