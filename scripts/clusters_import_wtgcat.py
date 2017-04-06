import numpy as np
import yaml
import pyfits
import astropy.table as table

# Create hdf5 Clusters pipeline-compatible file/astropy table from WtG catalogues
# Paths are hard-coded in this short script. Edit as necessary.

# WTG shear catalog where quality cuts have already been done
file_shear = '/Users/combet/RECHERCHE/LSST/RTF/WtG_catalogues/MACS2243-09.W-J-V.gab1728-rot1.cut_lensing.cat'
#file_shear = '/Users/combet/RECHERCHE/LSST/RTF/WtG_catalogues/MACS0018+16.W-J-V.gab1729-rot2.cut_lensing.cat'

# WTG redshift point estimates and pdz files
file_bpz = '/Users/combet/RECHERCHE/LSST/RTF/WtG_catalogues/MACS2243-09.W-J-V.bpz.tab'
file_bpzpdz = '/Users/combet/RECHERCHE/LSST/RTF/WtG_catalogues/MACS2243-09.W-J-V.pdz.cat'
#file_bpz = '/Users/combet/RECHERCHE/LSST/RTF/WtG_catalogues/MACS0018+16.W-J-V.bpz.tab'
#file_bpzpdz = '/Users/combet/RECHERCHE/LSST/RTF/WtG_catalogues/MACS0018+16.W-J-V.pdz.cat'

# Name of pipeline-ready output file
catfile = '/Users/combet/RECHERCHE/LSST/RTF/analysis/MACS2243/cat_wtg.hdf5'  # to be used with "clusters_mass.py config.yaml cat_wtg.hdf5"

# Extract Shear
cat_wtg = pyfits.open(file_shear)
ra = cat_wtg[4].data['ALPHA_J2000']
dec = cat_wtg[4].data['DELTA_J2000']
e1 = cat_wtg[4].data['gs1']
e2 = cat_wtg[4].data['gs2']
id1 = cat_wtg[4].data['SeqNr']

dict_cat_wtg={objidi:(rai,deci,e1i,e2i) for (objidi,rai,deci,e1i,e2i) in zip(id1,ra,dec,e1,e2)}

# Extract redshift information
cat_bpz_wtg = pyfits.open(file_bpz)
cat_pdz_wtg = pyfits.open(file_bpzpdz)


id2 = cat_bpz_wtg[1].data['SeqNr']
zb = cat_bpz_wtg[1].data['BPZ_Z_B']
dict_cat_zb_wtg={objidi:(zbi) for (objidi,zbi) in zip(id2,zb)}

id3 = cat_pdz_wtg[1].data['SeqNr']
pdz = cat_pdz_wtg[1].data['pdz']/cat_pdz_wtg[1].header['PDZSTEP']
dict_cat_pdz_wtg={objidi:(pdzi) for (objidi,pdzi) in zip(id3,pdz)}

zrange = np.arange(cat_pdz_wtg[1].header['MINPDZ'],cat_pdz_wtg[1].header['MAXPDZ'],cat_pdz_wtg[1].header['PDZSTEP'])
zbinsgrid = np.vstack(len(id1)*[zrange])

# Find which indices in cat_bpz_wtg and cat_pdz_wtg correspond to the objects in cat_wtg
good_SeqNr2 = []
[good_SeqNr2.append(np.where(id2 == id1[i])[0][0]) for i in np.arange(len(id1))]

good_SeqNr3 = []
[good_SeqNr3.append(np.where(id3 == id1[i])[0][0]) for i in np.arange(len(id1))]

# redshift flags
z_cl = 0.447
flag_hard = (zb[good_SeqNr2] > z_cl + 0.1) & (zb[good_SeqNr2] < 1.25)
cut = zrange < z_cl + 0.1
thresh = 1  # threshold probability for the galaxy to be located below z_cl + 0.1
flag_pdz = np.array([np.trapz(pdzi[cut], zrange[cut]) * 100. < thresh for pdzi in pdz[good_SeqNr3]])

# size and snratio, required when running pzmassfitter with STEP2 shear calibration
size = cat_wtg[4].data['rh']
snratio = cat_wtg[4].data['snratio_scaled1']

# Write pipeline-compatible hdf5 file
deepCoadd_meas = table.Table([id1, ra, dec, e1, e2, size, snratio], names=('id', 'coord_ra_deg', 'coord_dec_deg', 'ext_shapeHSM_HsmShapeRegauss_e1', 'ext_shapeHSM_HsmShapeRegauss_e2', 'rh', 'snratio_scaled1'))
pdz_values = table.Table([id1, zb[good_SeqNr2], pdz[good_SeqNr3], zbinsgrid], names=('objectId', 'Z_BEST', 'pdz', 'zbins'))
flag=table.Table([id1, flag_hard,flag_pdz], names=('objectId','flag_z_hard','flag_z_pdz'))

deepCoadd_meas.write(catfile, path='deepCoadd_meas', overwrite=True)
pdz_values.write(catfile, path='zphot_ref',append=True)
flag.write(catfile, path='flag_zphot_ref',append=True)
