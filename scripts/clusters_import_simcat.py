import pyfits
import astropy.table as table
import numpy as np
import pdb
# Create hdf5 Clusters pipeline-compatible file/astropy table from simulated catalogues
# Paths are hard-coded in this short script. Edit as necessary.

# WTG shear catalog where quality cuts have already been done
#file_shear = '/Users/combet/Desktop/src.fits'

# Name of pipeline-ready output file
srcfiles=[]
catfiles=[]
for i in np.arange(8):
    srcfile = '/Users/combet/RECHERCHE/LSST/RTF/analysis/SIM/DM_cat/src0'+str(i)+'.fits'
    srcfiles.append(srcfile)
    catfile = '/Users/combet/RECHERCHE/LSST/RTF/analysis/SIM/correct_ra/cat0'+str(i)+'_sim_leftra.hdf5'  # to be used with "clusters_mass.py config.yaml cat_wtg.hdf5"
    catfiles.append(catfile)


    
for i,fileshear in enumerate(srcfiles):
# Extract Shear
    cat_sim = pyfits.open(fileshear)
    ra_pix = cat_sim[1].data['base_GaussianCentroid_X']
    dec_pix = cat_sim[1].data['base_GaussianCentroid_y']
    ra = np.abs(ra_pix*0.2/3600.-0.2) # assumes a 0.2 arcmin pixel size
    dec = dec_pix*0.2/3600.
#    ra = cat_sim[1].data['coord_ra']
#    dec = cat_sim[1].data['coord_dec']
    e1 = cat_sim[1].data['ext_shapeHSM_HsmShapeRegauss_e1']
    e2 = cat_sim[1].data['ext_shapeHSM_HsmShapeRegauss_e2']
    obj_id = cat_sim[1].data['id']

# Write pipeline-compatible hdf5 file
    deepCoadd_meas = table.Table([obj_id, ra, dec, e1, e2], names=('id', 'coord_ra_deg', 'coord_dec_deg', 'ext_shapeHSM_HsmShapeRegauss_e1', 'ext_shapeHSM_HsmShapeRegauss_e2'))
    deepCoadd_meas.write(catfiles[i], path='deepCoadd_meas', overwrite=True)

    zsim = np.zeros(len(ra))+1.5
    data = np.array([obj_id, zsim])
    pdb.set_trace()
    np.savetxt('/Users/combet/RECHERCHE/LSST/RTF/analysis/SIM/sim0'+str(i)+'_z.txt', data.T, fmt=['%i','%f'])

#catfile = '/Users/combet/RECHERCHE/LSST/RTF/analysis/SIM/cat_sim_reshuffle.hdf5'  # to be used with "clusters_mass.py config.yaml cat_wtg.hdf5"
#np.random.shuffle(e1)
#np.random.shuffle(e2)

#deepCoadd_meas_fake = table.Table([id, ra, dec, e1, e2], names=('id', 'coord_ra_deg', 'coord_dec_deg', 'ext_shapeHSM_HsmShapeRegauss_e1', 'ext_shapeHSM_HsmShapeRegauss_e2'))

#deepCoadd_meas_fake.write(catfile, path='deepCoadd_meas', overwrite=True)

