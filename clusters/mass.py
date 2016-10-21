"""Mass analysis."""

import numpy
from . import data as data
from . import shear
from astropy.table import Column


def build_table():
    config = data.load_config('MACSJ2243.3-0935.yaml')
    catalogs = data.read_hdf5('MACSJ2243.3-0935_filtered_data.hdf5')
    table = catalogs['deepCoadd_meas']
    xclust, yclust = shear.xy_clust(config, data.load_wcs(catalogs['wcs']))
    e1r = table["ext_shapeHSM_HsmShapeRegauss_e1"][table['filter'] == 'r']
    e2r = table["ext_shapeHSM_HsmShapeRegauss_e2"][table['filter'] == 'r']
    e1i = table["ext_shapeHSM_HsmShapeRegauss_e1"][table['filter'] == 'i']
    e2i = table["ext_shapeHSM_HsmShapeRegauss_e2"][table['filter'] == 'i']
    distx = table["x_Src"][table['filter'] == 'r'] - xclust
    disty = table["y_Src"][table['filter'] == 'r'] - yclust
    objectid = table['id'][table['filter'] == 'i']

    # Quality cuts
    # magnitude cut
    filt = table['modelfit_CModel_mag'][table['filter'] == 'r'] < 23.5
    # resolution cut
    filt &= table['ext_shapeHSM_HsmShapeRegauss_resolution'][table['filter'] == 'i'] > 0.3
    # ellipticity cut
    filt &= (abs(e1i) < 1) & (abs(e2i) < 1)
    # er ~= ei
    filt &= (abs(e1r - e1i) < 0.5) & (abs(e2r - e2i) < 0.5)

    # Apply cuts
    e1i, e2i, distx, disty, objectid = [x[filt] for x in [e1i, e2i, distx, disty, objectid]]

    # Comput the shear
    gamt, gamc, dist = shear.compute_shear(e1i, e2i, distx, disty)

    zphot = data.read_hdf5("MACSJ2243.3-0935_filtered_data_zphot.hdf5")['zphot_zphot_megacam']
    params = zphot['coord_ra_deg', 'coord_dec_deg', 'Z_BEST', 'objectId']

    filt = [i for i, oid in enumerate(params['objectId']) if oid in objectid]
    params = params[filt]
    params.add_columns([Column(data=gamt, name='Tshear'), Column(data=gamc, name='Cshear')])
    params.write("data_yves.txt", format='ascii')
