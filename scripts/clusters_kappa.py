#!/usr/bin/env python
# Python code to take multiple shear mapping codes and a catalog from the CFHT
# clusters pipeline and produce multiple versions of 2-d maps

# read in the catalog: TODO--turn this code to accept a catalog from command
# line or from a config file
from clusters import data
import numpy as np
#f = "/home/chotard/Work/scripts/analysis/test_Cluster/MACSJ2243.3-0935_filtered_data.hdf5"
f = "/Volumes/clustersdata/MACSJ2243.3-0935_filtered_data.hdf5"
d = data.read_hdf5(f)
fc = d['deepCoadd_meas']

# read in the variables.  TODO add a switch that determines which shears to
# read in (or to do multiple ones.  For now we just use the regauss ones)
x = fc['x_Src']
y = fc['y_Src']
e1 = fc['ext_shapeHSM_HsmShapeRegauss_e1']
e2 = fc['ext_shapeHSM_HsmShapeRegauss_e2']
#
# do I need to convert these explicitly to numpy arrays?
#
# old code, ignore
#np.savetxt('filename', np.c_[x,y,e1,e2])
# also read in raw moments as a diagnostic
#ixx = fc['ext_shapeHSM_HsmSourceMoments_xx']
#ixy = fc['ext_shapeHSM_HsmSourceMoments_xy']
#iyy = fc['ext_shapeHSM_HsmSourceMoments_yy']

#old code, ignore
#np.savetxt('test2.txt', np.c_[x,y,ixx,iyy,ixy])

# define inner and outer cutoff radii for invlens algorithm
# TODO these should not be hardcoded, but should be input from the config
# file.  These are deweights to remove the quadratic noisedivergence for galaxies 
# with dist=0 and the logarithmic divergence for the noise from gals as r->infty.
# values are stored in pixels!
rinner = 500.0 
router = 8000.0
#
# set the subsampling step size.  This is the shrink factor by which the output
# image is reduced relative to the input image.  So a step of 50 will take a
# 10,000x10,000 image and turn it into a 200x200 map.
# TODO make this a parameter file input.
# TODO remember when we write out the WCS, that cd1_1,cd1_2,cd2_1,and cd2_2
# are multiplied by step, and crpix1 and crpix2 are divided by step!
#
step = 50
#
# now calculate the minimum and maximum x and y of catalog
# and use this to determine the size of image.
xmin = np.amin(x)
xmax = np.amax(x)
ymin = np.amin(y)
ymax = np.amax(y)
sizex = xmax - xmin
sizey = ymax - ymin 
nxpoints = int(sizex/step)
nypoints = int(sizey/step)
npoints = nxpoints * nypoints

# this sets up the weights for the invlens algorithm
# store the squares of rinner, router so we don't have to keep computing!
ri2 = rinner*rinner
ro2 = router*router
rmax = np.sqrt(sizex*sizex+sizey*sizey)
# no objects are ever separated by more than rmax, so we never need to
# store cutoff weights for r> rmax as an array so we don't have to calculate it on the fly.
irmax = int(rmax)
# create an empty weight array for each method.  We want emty weights so that odd values do not
# get imported by accident
wt = np.zeros(irmax)
wtmat = np.zeros(irmax)
wtpot = np.zeros(irmax)
wtint = np.zeros(irmax)
wtapmass = np.zeros(irmax)
#
# this part of the code sets up the Maturi et al filter for aperture mass
tanha=6.0
tanhb=150.0
tanhc=50.0
tanhd=47.0
tanhxc=0.1
theta0 = 3000.0
# TODO: the theta0 parameter should be input from the configuration file
#
#
# The Aperture mass radius from Schneider et al. 1998 is set to 3000 TODO: put it in the config file
#
aprad = 3000
apr2 = aprad*aprad

# now populate all the weight arrays

for i in range(irmax):
 r2 = i
 icut = 1 - np.exp(-r2/(2.0*ri2))
 ocut = np.exp(-r2/(2.0*ro2))
 wt[i] = icut * ocut / r2;
 tanhx = np.sqrt(r2/(rinner*rinner))
 icut2 = np.tanh(tanhx/tanhxc)/((tanhx/tanhxc)*(1+np.exp(tanha-tanhb*tanhx)+np.exp(tanhc*tanhx-tanhd)))
 wtmat[i] = icut2
 wtpot[i] = icut * ocut / np.sqrt(r2)
 wtint[i] = icut * ocut
 wtapmass[i] = 6.0/np.pi * (1.0 - r2/apr2)*(r2/apr2) 
 if (r2 > apr2):
  wtapmass[i] = 0
 if (r2 > (theta0*theta0)):
  wtmat[i] = 0
#
# now let's actually calculate the shear.  The algorithm for all of these images is pretty simple.  
# 1) make a grid of (x,y) locations
# 2) at each location, loop over all galaxies and calculate the tangential shear (and 45-degree shear)
#    calculated from a rotation of e_1 and e_2.  First, the direction connecting the location and the 
#    position is calculated, and the sin(2*phi) and cos(2*phi) terms are  calculated.
#    then    etan = -1.0*( e1 * cos2phi + e2 * sin2phi)
#    and   ecross = ( e2 * cos2phi - e1 * sin2phi)
# 3) The rotated ellipticities are multiplied by a Weight, which depends on the type of map.  
# 4) The sum of the weighted ellipticities is divided by the sum of the weights to provide an output
# first, loop over  x and y points in output map
  #
#nxpoints, nypoints = 50, 50
invlensmap=np.zeros((nxpoints,nypoints))
inv45map=np.zeros((nxpoints,nypoints))
maturi=np.zeros((nxpoints,nypoints))
maturi45=np.zeros((nxpoints,nypoints))
apmassmap=np.zeros((nxpoints,nypoints))
apmass45map=np.zeros((nxpoints,nypoints))
potmap=np.zeros((nxpoints,nypoints))
pot45map=np.zeros((nxpoints,nypoints))
intmap=np.zeros((nxpoints,nypoints))
int45map=np.zeros((nxpoints,nypoints))
for nxp in range(nxpoints):
 yp = ymin + (nxp+0.5)*step
 for nyp in range(nypoints):
  xp = xmin + (nyp+0.5)
# now loop over all the objects in your catalog
# the new version of the code treats all objects as numpy arrays implicitly.  Will it work?
  dx = x - xp
  dy = y - yp
  r2 = dx * dx + dy * dy
  r = np.sqrt(r2)
  nr = np.array(r.tolist(),dtype=int)
  cos2phi = (dx*dx - dy*dy) / r2
  sin2phi = 2.0*dx*dy / r2
  etan = -1.0*( e1 * cos2phi + e2 * sin2phi)
  ecross = ( e2 * cos2phi - e1 * sin2phi)
  invlens = wt[nr] * etan
  inv45 = wt[nr] * ecross
  mat = wtmat[nr] * etan
  mat45 = wtmat[nr] * ecross
  apmass = wtapmass[nr] * etan
  apmass45 = wtapmass[nr] * ecross
  pot = wtpot[nr] * etan
  pot45 = wtpot[nr] * ecross
  integ = wtint[nr] * etan
  int45 = wt[nr] * ecross
  #end loop over objects
  invlensmap[nxp][nyp] = np.sum(invlens) /  np.sum(wt[nr])
  inv45map[nxp][nyp] = np.sum(inv45) / np.sum(wt[nr])
  maturi[nxp][nyp] = np.sum(mat) / np.sum(wtmat[nr])
  maturi45[nxp][nyp] = np.sum(mat45) / np.sum(wtmat[nr])
  apmassmap[nxp][nyp] = np.sum(apmass) / np.sum(wtapmass[nr])
  apmass45map[nxp][nyp] = np.sum(apmass45) / np.sum(wtapmass[nr])
  potmap[nxp][nyp] = np.sum(pot) / np.sum(wtpot[nr])
  pot45map[nxp][nyp] = np.sum(pot45) / np.sum(wtpot[nr])
  intmap[nxp][nyp] = np.sum(integ) / np.sum(wtint[nr])
  int45map[nxp][nyp] = np.sum(int45) / np.sum(wtint[nr])
  #end x loop
 #end y loop
#
# Now we have to write the files out as fits files.
import astropy.io.fits as pyfits
hdu = pyfits.PrimaryHDU(invlensmap)
# to modify header, use this syntax
#hdu.header["cd1_1"]=XXXX
hdu.writeto("invlens.fits",clobber=True)

hdu = pyfits.PrimaryHDU(inv45map)
hdu.writeto("inv45.fits", clobber=True)

hdu = pyfits.PrimaryHDU(maturi)
hdu.writeto("maturi.fits", clobber=True)

hdu = pyfits.PrimaryHDU(maturi45)
hdu.writeto("maturi45.fits", clobber=True)

hdu = pyfits.PrimaryHDU(apmassmap)
hdu.writeto("apmass.fits", clobber=True)

hdu = pyfits.PrimaryHDU(apmass45map)
hdu.writeto("apmass45.fits", clobber=True)

hdu = pyfits.PrimaryHDU(potmap)
hdu.writeto("potential.fits", clobber=True)

hdu = pyfits.PrimaryHDU(pot45map)
hdu.writeto("potential45.fits", clobber=True)

hdu = pyfits.PrimaryHDU(intmap)
hdu.writeto("integralshear.fits", clobber=True)

hdu = pyfits.PrimaryHDU(int45map)
hdu.writeto("integralshear45.fits", clobber=True)
