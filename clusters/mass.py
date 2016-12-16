"""Mass analysis."""

import seaborn
import pylab
import pandas
import cPickle
from astropy.table import Column
from . import data as data
from . import shear


def build_table(config='MACSJ2243.3-0935.yaml', datafile='MACSJ2243.3-0935_filtered_data.hdf5'):
    """Load data and write an ascii file containing enough info to run Yves mass estimator."""
    table = data.read_hdf5(datafile, path='deepCoadd_meas', dic=False)
    xclust, yclust = shear.xy_clust(data.load_config(config),
                                    data.load_wcs(data.read_hdf5(datafile, path='wcs', dic=False)))
    e1i = table["ext_shapeHSM_HsmShapeRegauss_e1"][table['filter'] == 'i']
    e2i = table["ext_shapeHSM_HsmShapeRegauss_e2"][table['filter'] == 'i']
    sigmai = table["ext_shapeHSM_HsmShapeRegauss_sigma"][table['filter'] == 'i']
    distx = table["x_Src"][table['filter'] == 'r'] - xclust
    disty = table["y_Src"][table['filter'] == 'r'] - yclust
    objectid = table['id'][table['filter'] == 'i']

    # Apply cuts
    e1i, e2i, distx, disty, objectid, sigmai = [x[quality_cuts(table)] for x in [e1i, e2i, distx,
                                                                                 disty, objectid,
                                                                                 sigmai]]

    params = data.read_hdf5("MACSJ2243.3-0935_filtered_data_zphot.hdf5",
                            path='bpz', dic=False)['coord_ra_deg', 'coord_dec_deg', 'Z_B',
                                                   'objectId', 'Z_B_MIN', 'Z_B_MAX']

    params = params[[i for i, oid in enumerate(params['objectId']) if oid in objectid]]
    params.add_columns([Column(data=shear.compute_shear(e1i, e2i, distx, disty)[0], name='Tshear'),
                        Column(data=shear.compute_shear(e1i, e2i, distx, disty)[1], name='Cshear'),
                        Column(data=e1i, name='Ell-1-ir'),
                        Column(data=e2i, name='Ell-2-i'),
                        Column(data=sigmai, name='Sigma-ell')])
    params.write("data_yves.txt", format='ascii')


def quality_cuts(table):
    """Apply some quality cuts."""
    rfilter = table['filter'] == 'r'
    ifilter = table['filter'] == 'r'
    filt = table['modelfit_CModel_mag'][rfilter] < 23.5
    # resolution cut
    filt &= table['ext_shapeHSM_HsmShapeRegauss_resolution'][ifilter] > 0.3
    # ellipticity cut
    filt &= (abs(table["ext_shapeHSM_HsmShapeRegauss_e1"][ifilter]) < 1) & \
            (abs(table["ext_shapeHSM_HsmShapeRegauss_e2"][ifilter]) < 1)
    # er ~= ei
    filt &= (abs(table["ext_shapeHSM_HsmShapeRegauss_e1"][rfilter] - \
                 table["ext_shapeHSM_HsmShapeRegauss_e1"][ifilter]) < 0.5) & \
            (abs(table["ext_shapeHSM_HsmShapeRegauss_e2"][rfilter] - \
                 table["ext_shapeHSM_HsmShapeRegauss_e2"][ifilter]) < 0.5)
    return filt


def plot_pzmassfitter_output(datafile):
    """Te datafile is the ***.chain.pkl file."""
    d = cPickle.load(open(datafile))
    df = pandas.DataFrame(d)
    g = seaborn.PairGrid(df, diag_sharey=False)
    g.map_lower(pylab.scatter)
    g.map_diag(pylab.hist, bins=50)
    
    #x = pandas.DataFrame(d['mdelta'])
    #y = pandas.DataFrame(d['log10concentration'])
    #seaborn.jointplot(x, y, kind="kde", size=7, space=0)
    pylab.show()
    