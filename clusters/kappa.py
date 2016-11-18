"""Kappa analysis."""

import numpy as np
import astropy.io.fits as pyfits
from . import data as cdata


class Kappa(object):

    """Kappa analysis of a shear map."""
    def __init__(self, xsrc, ysrc, sch1, sch2, **kwargs):
        """Kappa analysis of a shear map.

        :param list xsrc: x coordinates
        :param list xsrc: y coordinates
        :param list sch1: schear 1
        :param list sch2: schear 2
        :param list filt: filter to apply
        """
        assert len(xsrc) == len(ysrc) == len(sch1) == len(sch2)

        # Make sure all list are actually numpy arrays
        xsrc, ysrc, sch1, sch2 = [np.array(x) for x in [xsrc, ysrc, sch1, sch2]]

        # Do we have to filter?
        if kwargs.get('filt', None) is not None:
            assert len(kwargs['filt']) == len(xsrc)
            self._idata = {'xsrc': xsrc, 'ysrc': ysrc,
                           'sch1': sch1, 'sch2': sch2,
                           'flag': np.array(kwargs.get('filt'))}
            print "Filtering data"
            print " - Before cut: %i sources" % len(xsrc)
            xsrc, ysrc, sch1, sch2 = [arr[np.array(kwargs.get('filt'))]
                                      for arr in [xsrc, ysrc, sch1, sch2]]
            print " - After cut: %i sources" % len(xsrc)
            self.filtered = True
        else:
            self._idata = None
            self.filtered = False

        # Save the data
        self.data = {'xsrc': xsrc, 'ysrc': ysrc,
                     'sch1': sch1, 'sch2': sch2}

        # Get or compute some useful parameters
        self.parameters = {
            # size of the image
            'sizex': max(xsrc) - min(xsrc),
            'sizey': max(ysrc) - min(ysrc),
            # step, useful for test purpose
            'step': kwargs.get('step', 1),
            # number of points for each axis
            'nxpoints': int((max(xsrc) - min(xsrc)) / kwargs.get('step', 1)),
            'nypoints': int((max(ysrc) - min(ysrc)) / kwargs.get('step', 1)),
            # inner and outer radius
            'rinner': kwargs.get('rinner', 500.0),
            'router': kwargs.get('router', 8000.0),
            # maximum radius
            'rmax': np.sqrt((max(xsrc) - min(xsrc))**2 + (max(ysrc) - min(ysrc))**2),
            'theta0': kwargs.get('theta0', 6000.0),
            'aprad': kwargs.get('aprad', 6000.0)
        }

        # Define needed dictionnary
        self.maps = self.weights = {}

        # Get the weights and the maps
        self.get_weights()
        self._get_kappa()

    def get_weights(self):
        """Set up the weights for the invlens algorithm.

        No objects are ever separated by more than rmax, so we never need to store cutoff weights
        for r > rmax as an array so we don't have to calculate it on the fly.
        """
        r2 = np.arange(int(self.parameters['rmax'] + 0.5))**2
        icut = 1 - np.exp(-r2 / (2.0 * self.parameters['rinner']**2))
        ocut = np.exp(-r2 / (2.0 * self.parameters['router']**2))
        icut2 = aperture_mass_maturi_filter(np.sqrt(r2 / self.parameters['rinner']**2))
        # Is the following formula the right one? 6.0 / pi * ... or (6 / pi) * ...?
        wtapmass = 6.0 / np.pi * (1.0 - r2/self.parameters['aprad']**2) * \
                   (r2/self.parameters['aprad']**2)

        # Store them
        self.weights = {"invlens": icut * ocut / r2,
                        "invlens45": icut * ocut / r2,
                        "maturi": icut2,
                        "maturi45": icut2,
                        "apmass": wtapmass,
                        "apmass45": wtapmass,
                        "potential": icut * ocut / np.sqrt(r2),
                        "potential45": icut * ocut / np.sqrt(r2),
                        "integralshear": icut * ocut,
                        "integralshear45": icut * ocut}

        # Other settings
        for weight in self.weights:
            if weight.startswith('maturi'):
                self.weights[weight][np.argwhere(r2 > self.parameters['theta0']**2)] = 0
            if weight.startswith('apmass'):
                self.weights[weight][np.argwhere(r2 > self.parameters['aprad']**2)] = 0
            else:
                self.weights[weight][0] = 0

    def _get_axis_3dgrid(self, axis='x'):
        """Construct a grid around the data map.

        Move the grid so its botom left corner match the galaxy which is at the bottom left corner
        Here, we build a grid either for the x or the y position.
        """
        if axis == 'x':
            cmin = min(self.data['xsrc'])
            carange = np.arange(self.parameters['nxpoints']) + 0.5
            ncpoints = self.parameters['nxpoints']
            nopoints = self.parameters['nypoints']
        else:
            cmin = min(self.data['ysrc'])
            carange = np.arange(self.parameters['nypoints']) + 0.5
            ncpoints = self.parameters['nypoints']
            nopoints = self.parameters['nxpoints']
        cgrid = (cmin +  carange * self.parameters['step'])
        #cgrid = cgrid.reshape(1, ncpoints).repeat(nopoints, axis=0).reshape(nopoints, ncpoints, 1)
        if axis == 'x':
            cgrid = cgrid.reshape(1, ncpoints, 1)
        else:
            cgrid = cgrid.reshape(ncpoints, 1, 1)
        #import gc
        #gc.collect()
        #print np.shape(cgrid)
        return cgrid

    def _get_distances(self):
        pass

    def _get_kappa(self):
        """Now let's actually calculate the shear.

        The algorithm for all of these images is pretty simple.
          1) make a grid of (x,y) locations
          2) at each location, loop over all galaxies and calculate the tangential shear
             (and 45-degree shear) calculated from a rotation of e_1 and e_2.  First, the direction
             connecting the location and the position is calculated, and the sin(2*phi) and
             cos(2*phi) terms are  calculated. Then
                cos2phi = (dx*2 - dy**2) / squared_radius
                sin2phi = 2.0 * dx * dy / squared_radius
                etan = -1.0*( e1 * cos2phi + e2 * sin2phi)
                ecross = ( e2 * cos2phi - e1 * sin2phi)
          3) The rotated ellipticities are multiplied by a Weight, which depends on the type of map.
          4) The sum of the weighted ellipticities is divided by the sum of the weights to
             provide an output
        """
        # Compute all distance. We get a cube of distances. For each point of the grid, we have an
        # array of distances to all the sources of the catalog. This is a 3d array.
        dx = self.data['xsrc'].reshape(1, 1, len(self.data['xsrc'])) - self._get_axis_3dgrid(axis='x')
        print np.shape(dx)
        #dx = self.data['xsrc'].reshape(1, 1, len(self.data['xsrc']))  # first reshape to 3d
        #print np.shape(dx)
        #dx -= self._get_axis_3dgrid(axis='x')  # then subtract the x 3d grid
        dy = self.data['ysrc'].reshape(1, 1, len(self.data['ysrc'])) - self._get_axis_3dgrid(axis='y')
        #r2 = dx**2 - dy**2
        #cos2phi = (dx**2 - dy**2) / r2
        #sin2phi = 2.0 * dx * dy / r2
        self.dx = dx
        self.dy = dy
        
        return
        # Compute rotated ellipticities
        etan = -1.0 * (self.data['sch1'] * (dx**2 - dy**2) / (dx**2 + dy**2) + \
                       self.data['sch2'] * 2.0 * dx * dy / (dx**2 + dy**2))
        ecross = -1.0 * (self.data['sch2'] * (dx**2 - dy**2) / (dx**2 + dy**2) - \
                         self.data['sch1'] * 2.0 * dx * dy / (dx**2 + dy**2))

        # Transform the cube of distances into a cube of integer, it will serve a indexes for weight
        int_radius = np.array(np.sqrt(dx**2 + dy**2).tolist(), dtype=int)

        # Apply this indexes filter to get the right weigth array for each pixel of the grid
        # The output is also a 3D array (nxpoints, nypoints, len(x))
        weights = {cmap: self.weights[cmap][int_radius] for weight in self.weights}

        self.maps = {}
        for cmap in self.weights:
            sumw = np.sum(weights[cmap], axis=2)
            if '45' in cmap:
                self.maps[cmap] = np.sum(weights[cmap] * etan, axis=2) / sumw
            else:
                self.maps[cmap] = np.sum(weights[cmap] * ecross, axis=2) / sumw

    def plot_maps(self):
        pass

    def save_maps(self, format='fits'):
        # Now write the files out as fits files
        # To modify header, use this syntax: hdu.header["cd1_1"]=XXXX
        for cmap in self.maps:
            hdu = pyfits.PrimaryHDU(self.maps[cmap])
            hdu.writeto("%s.fits" % cmap, clobber=True)


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
    print "Get the main filter"
    filt = (abs(cat[e1]) < 1.2) & (abs(cat[e2] < 1.2))
    if 'z_flag_pdz' in cat.keys():
        print "Add pdz filter"
        filt &= cat['z_flag_pdz']
    return [cat[k][filt] for k in ['x_Src', 'y_Src', e1, e2]]


def aperture_mass_maturi_filter(tanhx, **kwargs):
    """Maturi et al filter for aperture mass."""
    tanha = kwargs.get("tanha", 6.0)
    tanhb = kwargs.get('tanhb', 150.0)
    tanhc = kwargs.get('tanhc', 50.0)
    tanhd = kwargs.get('tanhd', 47.0)
    tanhxc = kwargs.get('tanhxc', 0.1)
    return np.tanh(tanhx / tanhxc) / ((tanhx / tanhxc) * \
                                         (1 + np.exp(tanha - tanhb * tanhx) + \
                                          np.exp(tanhc * tanhx-tanhd)))


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
    print "Loading the data"
    cat = load_data(datafile)

    # Get the data from the catalog
    print "Getting the needed arrays"
    x, y, e1, e2 = get_data(cat)

    # Compute the size of the image in pixel
    print "Compute the size"
    sizex = max(x) - min(x)
    sizey = max(y) - min(y)

    # Determine the size and number of points in the image
    print "Determine the size and number of points in the image"
    nxpoints = int(sizex / step)
    nypoints = int(sizey / step)

    # this sets up the weights for the invlens algorithm
    # store the squares of rinner, router so we don't have to keep computing!
    ri2 = rinner*rinner
    ro2 = router*router
    rmax = np.sqrt(sizex*sizex+sizey*sizey)
    # no objects are ever separated by more than rmax, so we never need to
    # store cutoff weights for r> rmax as an array so we don't have to calculate it on the fly.
    irmax = int(rmax+0.5)
    # create an empty weight array for each method.  We want emty weights so that
    # odd values do not get imported by accident
    wt = np.zeros(irmax)
    wtmat = np.zeros(irmax)
    wtpot = np.zeros(irmax)
    wtint = np.zeros(irmax)
    wtapmass = np.zeros(irmax)

    # The Aperture mass radius from Schneider et al. 1998 is set to 3000
    apr2 = aprad * aprad

    # now populate all the weight arrays
    print "Now populate all the weight arrays"
    for i in range(irmax):
        r2 = i * i
        icut = 1 - np.exp( -r2 / (2.0 * ri2))
        ocut = np.exp( -r2 / (2.0 * ro2))
        icut2 = aperture_mass_maturi_filter(np.sqrt(r2 / (rinner * rinner)))
        wt[i] = icut * ocut / r2
        wtmat[i] = icut2
        wtpot[i] = icut * ocut / np.sqrt(r2)
        wtint[i] = icut * ocut
        wtapmass[i] = 6.0/np.pi * (1.0 - r2/apr2)*(r2/apr2) 
        if (r2 > apr2):
            wtapmass[i] = 0
        if (r2 > (theta0*theta0)):
            wtmat[i] = 0

    wt[0]=0
    wtpot[0]=0
    wtint[0]=0
    wtmat[0]=0

    # now let's actually calculate the shear.  The algorithm for all of these images is pretty simple.
    # 1) make a grid of (x,y) locations
    # 2) at each location, loop over all galaxies and calculate the tangential shear
    # (and 45-degree shear) calculated from a rotation of e_1 and e_2.  First, the direction
    # connecting the location and the position is calculated, and the sin(2*phi) and cos(2*phi) terms
    # are  calculated. then
    #   etan = -1.0*( e1 * cos2phi + e2 * sin2phi)
    #   ecross = ( e2 * cos2phi - e1 * sin2phi)
    # 3) The rotated ellipticities are multiplied by a Weight, which depends on the type of map.  
    # 4) The sum of the weighted ellipticities is divided by the sum of the weights to
    # provide an output

    
    # first, loop over  x and y points in output map
    mapnames = ["invlens", "invlens45", "maturi", "maturi45", "apmass", "apmass45",
                "potential", "potential45", "integralshear", "integralshear45"]
    maps = {key: np.zeros((nypoints,nxpoints)) for key in mapnames}
    print "Now loop over  x and y points in output map"
    for nxp in range(nxpoints):
        print nxp+1, nxpoints
        xp = min(x) + (nxp+0.5)*step
        for nyp in range(nypoints):
            yp = min(y) + (nyp+0.5)*step
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
            maps['invlens'][nyp][nxp] = np.sum(invlens) /  np.sum(wt[nr])
            maps['invlens45'][nyp][nxp] = np.sum(inv45) / np.sum(wt[nr])
            maps['maturi'][nyp][nxp] = np.sum(mat) / np.sum(wtmat[nr])
            maps['maturi45'][nyp][nxp] = np.sum(mat45) / np.sum(wtmat[nr])
            maps['apmass'][nyp][nxp] = np.sum(apmass) / np.sum(wtapmass[nr])
            maps['apmass45'][nyp][nxp] = np.sum(apmass45) / np.sum(wtapmass[nr])
            maps['potential'][nyp][nxp] = np.sum(pot) / np.sum(wtpot[nr])
            maps['potential45'][nyp][nxp] = np.sum(pot45) / np.sum(wtpot[nr])
            maps['integralshear'][nyp][nxp] = np.sum(integ) / np.sum(wtint[nr])
            maps['integralshear45'][nyp][nxp] = np.sum(int45) / np.sum(wtint[nr])
            #end x loop
        #end y loop
    #
    # Now write the files out as fits files
    # To modify header, use this syntax: hdu.header["cd1_1"]=XXXX
    for cmap in maps:
        hdu = pyfits.PrimaryHDU(maps[cmap])
        hdu.writeto("%s.fits" % cmap, clobber=True)
