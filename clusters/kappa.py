"""Shear analysis."""

import numpy
import pylab
from astropy.table import Table, Column
from . import data as cdata


def load_data(datafile):
    """Load the needed data."""
    return cdata.read_hdf5(datafile, path='deepCoadd_meas', dic=False)


def get_data(cat, e1='ext_shapeHSM_HsmShapeRegauss_e1', e2='ext_shapeHSM_HsmShapeRegauss_e2'):
    """Read in the variables.

    Input: astropy table for the deepcoad_meas catalog
    Return: x, y, e1, e2

    .. TODO::  add a switch that determines which shears to read in
    (or to do multiple ones.  For now we just use the regauss ones.
    """
    filt = (abs(cat[e1]) < 1.2) & (abs(cat[e2] < 1.2))
    if 'z_flag_pdz' in cat:
        filt &= cat['z_flag_pdz']
    # Comput ethe size of the image in pixel
    sx = max(cat['x_Src']) - min(cat['x_Src'])
    sy = max(cat['y_Src']) - min(cat['y_Src'])
    return [cat[k][filt] for k in ['x_Src', 'y_Src', e1, e2]] + [sx[filt], sy[filt]]


def aperture_mass_maturi_filter(tanhx, **kwargs):
    """Maturi et al filter for aperture mass."""
    tanha = kwargs.get("tanha", 6.0)
    tanhb = kwargs.get('tanhb', 150.0)
    tanhc = kwargs.get('tanhc', 50.0)
    tanhd = kwargs.get('tanhd', 47.0)
    tanhxc = kwargs.get('tanhxc', 0.1)
    return numpy.tanh(tanhx / tanhxc) / ((tanhx / tanhxc) * \
                                         (1 + numpy.exp(tanha - tanhb * tanhx) + \
                                          numpy.exp(tanhc * tanhx-tanhd)))


def analysis(datafile, rinner=500.0, router=8000.0, step=200, theta0=6000.0,
             aprad=6000.0):
    """Compute the kappa maps.

    Define inner and outer cutoff radii for invlens algorithm file. These are
    deweights to remove the quadratic noise divergence for galaxies right at the
    points being considered and the logarithmic divergence for the noise from
    gals as r->infty. Values are stored in pixels!
    :param int step: Subsampling step size.  This is the shrink factor by which the
                     output image is reduced relative to the input image.  So a step
                     of 50 will take a 10,000x10,000 image and turn it into a 200x200
                     map.
                     .. TODO:: remember when we write out the WCS, that
                     cd1_1,cd1_2,cd2_1,and cd2_2 are multiplied by step, and
                     crpix1 and crpix2 are divided by step!
                     Step set to 200 to make testing faster
    :param frloat aprad: Aperture mass radius from Schneider et al. 98 is set to 3000

    """
    # Load the catalog
    cat = load_data(datafile)

    # Get the data from the catalog
    x, y, e1, e2, sizex, sizey = get_data(cat)

    # Determine the size and number of points in the image
    nxpoints = int(sizex / step)
    nypoints = int(sizey / step)

    # this sets up the weights for the invlens algorithm
    # store the squares of rinner, router so we don't have to keep computing!
    ri2 = rinner*rinner
    ro2 = router*router
    rmax = numpy.sqrt(sizex*sizex+sizey*sizey)
    # no objects are ever separated by more than rmax, so we never need to
    # store cutoff weights for r> rmax as an array so we don't have to calculate it on the fly.
    irmax = int(rmax+0.5)
    # create an empty weight array for each method.  We want emty weights so that
    # odd values do not get imported by accident
    wt = numpy.zeros(irmax)
    wtmat = numpy.zeros(irmax)
    wtpot = numpy.zeros(irmax)
    wtint = numpy.zeros(irmax)
    wtapmass = numpy.zeros(irmax)

    # The Aperture mass radius from Schneider et al. 1998 is set to 3000
    apr2 = aprad * aprad

    # now populate all the weight arrays
    for i in range(irmax):
        r2 = i * i
        icut = 1 - numpy.exp( -r2 / (2.0 * ri2))
        ocut = numpy.exp( -r2 / (2.0 * ro2))
        icut2 = aperture_mass_maturi_filter(numpy.sqrt(r2 / (rinner * rinner)))
        wt[i] = icut * ocut / r2
        wtmat[i] = icut2
        wtpot[i] = icut * ocut / numpy.sqrt(r2)
        wtint[i] = icut * ocut
        wtapmass[i] = 6.0/numpy.pi * (1.0 - r2/apr2)*(r2/apr2) 
        if (r2 > apr2):
            wtapmass[i] = 0
        if (r2 > (theta0*theta0)):
            wtmat[i] = 0

    wt[0]=0
    wtpot[0]=0
    wtint[0]=0
    wtmat[0]=0
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
    invlensmap=numpy.zeros((nypoints,nxpoints))
    inv45map=numpy.zeros((nypoints,nxpoints))
    maturi=numpy.zeros((nypoints,nxpoints))
    maturi45=numpy.zeros((nypoints,nxpoints))
    apmassmap=numpy.zeros((nypoints,nxpoints))
    apmass45map=numpy.zeros((nypoints,nxpoints))
    potmap=numpy.zeros((nypoints,nxpoints))
    pot45map=numpy.zeros((nypoints,nxpoints))
    intmap=numpy.zeros((nypoints,nxpoints))
    int45map=numpy.zeros((nypoints,nxpoints))
    for nxp in range(nxpoints):
        xp = min(x) + (nxp+0.5)*step
        for nyp in range(nypoints):
            yp = min(y) + (nyp+0.5)*step
            # now loop over all the objects in your catalog
            # the new version of the code treats all objects as numpy arrays implicitly.  Will it work?
            dx = x - xp
            dy = y - yp
            r2 = dx * dx + dy * dy
            r = numpy.sqrt(r2)
            nr = numpy.array(r.tolist(),dtype=int)
            cos2phi = (dx*dx - dy*dy) / r2
            sin2phi = 2.0*dx*dy / r2
            etan = -1.0*( e1 * cos2phi + e2 * sin2phi)
            ecross = -1.0 * ( e2 * cos2phi - e1 * sin2phi)
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
            invlensmap[nyp][nxp] = numpy.sum(invlens) /  numpy.sum(wt[nr])
            inv45map[nyp][nxp] = numpy.sum(inv45) / numpy.sum(wt[nr])
            maturi[nyp][nxp] = numpy.sum(mat) / numpy.sum(wtmat[nr])
            maturi45[nyp][nxp] = numpy.sum(mat45) / numpy.sum(wtmat[nr])
            apmassmap[nyp][nxp] = numpy.sum(apmass) / numpy.sum(wtapmass[nr])
            apmass45map[nyp][nxp] = numpy.sum(apmass45) / numpy.sum(wtapmass[nr])
            potmap[nyp][nxp] = numpy.sum(pot) / numpy.sum(wtpot[nr])
            pot45map[nyp][nxp] = numpy.sum(pot45) / numpy.sum(wtpot[nr])
            intmap[nyp][nxp] = numpy.sum(integ) / numpy.sum(wtint[nr])
            int45map[nyp][nxp] = numpy.sum(int45) / numpy.sum(wtint[nr])
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