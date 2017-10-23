"""Kappa analysis."""


from __future__ import print_function
import sys
import numpy as np
import pylab as pl
import seaborn
try:
    import numba
except ImportError:
    print("WARNING: optional module numba cannot be imported.")
import astropy.io.fits as pyfits
try:
    import numexpr
except ImportError:
    print("WARNING: optional module numexpr cannot be imported.")
from . import utils as cutils


seaborn.set(style='ticks')

class Kappa(object):

    """Kappa analysis of a shear map."""

    def __init__(self, xsrc, ysrc, sch1, sch2, **kwargs):
        """Kappa analysis of a shear map.

        :param list xsrc: x coordinates
        :param list xsrc: y coordinates
        :param list sch1: schear 1
        :param list sch2: schear 2
        :param list filt: filter to apply

        Possible kwargs:

        :param int step: Step of the grid
        :param float rinner: Inner cut in pixel
        :param float router: Outer cut in pixel
        :param float theta0: Maturi outer cut in pixel
        :param float aprad: Aperture mass outer cut in pixel
        """
        assert len(xsrc) == len(ysrc) == len(sch1) == len(sch2)

        # numba?
        self.use_numba = kwargs.get("numba", False) and 'numba' in sys.modules
        # numexpr?
        self.use_numexpr = kwargs.get("numexpr", False) and 'numexpr' in sys.modules
        if self.use_numba:
            for func in to_jit:
                apply_decorator(func, numba.jit)
            for func in to_vectorize:
                apply_decorator(func, numba.vectorize)
        elif self.use_numexpr:
            numexpr.set_num_threads(numexpr.detect_number_of_threads())

        # Make sure all list are actually numpy arrays
        xsrc, ysrc, sch1, sch2 = [np.array(x.tolist()) for x in [
            xsrc, ysrc, sch1, sch2]]

        # Do we have to filter?
        if kwargs.get('filt', None) is not None:
            assert len(kwargs['filt']) == len(xsrc)
            self._idata = {'xsrc': xsrc, 'ysrc': ysrc,
                           'sch1': sch1, 'sch2': sch2,
                           'flag': np.array(kwargs.get('filt'))}
            print("Filtering data")
            print(" - Before cut: %i sources" % len(xsrc))
            xsrc, ysrc, sch1, sch2 = [arr[np.array(kwargs.get('filt'))]
                                      for arr in [xsrc, ysrc, sch1, sch2]]
            print(" - After cut: %i sources" % len(xsrc))
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
        self._get_weights()
        self._get_kappa(xsampling=kwargs.get('xsampling', 10))
        self.save_maps()

    def _get_weights(self):
        """Set up the weights for the invlens algorithm.

        No objects are ever separated by more than rmax, so we never need to store cutoff weights
        for r > rmax as an array so we don't have to calculate it on the fly.
        """
        r2 = np.arange(int(self.parameters['rmax'] + 0.5))**2
        icut = 1 - np.exp(-r2 / (2.0 * self.parameters['rinner']**2))
        ocut = np.exp(-r2 / (2.0 * self.parameters['router']**2))
        icut2 = aperture_mass_maturi_filter(
            np.sqrt(r2 / self.parameters['rinner']**2))
        # Is the following formula the right one? 6.0 / pi * ... or (6 / pi) * ...?
        wtapmass = 6.0 / np.pi * (1.0 - r2 / self.parameters['aprad']**2) * \
            (r2 / self.parameters['aprad']**2)

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
                self.weights[weight][np.argwhere(
                    r2 > self.parameters['theta0']**2)] = 0
            if weight.startswith('apmass'):
                self.weights[weight][np.argwhere(
                    r2 > self.parameters['aprad']**2)] = 0
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
            reshapep = (self.parameters['nxpoints'], 1)
        else:
            cmin = min(self.data['ysrc'])
            carange = np.arange(self.parameters['nypoints']) + 0.5
            reshapep = (self.parameters['nypoints'], 1)
        cgrid = (cmin + carange * self.parameters['step']).reshape(reshapep)
        return cgrid

    def _get_kappa(self, xsampling=100000):
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
        dx = self.data['xsrc'].reshape(
            1, len(self.data['xsrc'])) - self._get_axis_3dgrid(axis='x')
        dy = self.data['ysrc'].reshape(
            1, len(self.data['ysrc'])) - self._get_axis_3dgrid(axis='y')

        if self.use_numba:
            print("Using numba to slightly speed up the process!")
        if self.use_numexpr:
            print("Using numexpr to slightly speed up the process!")
        # also loop over the x axis to pack them into arrays of 'xsampling' items
        xarange = [i for i in range(len(dx)) if not i % xsampling] + [len(dx)]
        dxs = [dx[xarange[jj]:xarange[jj + 1]]
               for jj in range(len(xarange[:-1]))]
        # Now loop over the y axis (explode memory otherwise)
        pbar = cutils.progressbar(len(dy) * (len(xarange) - 1))
        for ii, dyy in enumerate(dy):
            etan, ecross, int_radius = [], [], []
            for jj, dxx in enumerate(dxs):
                if self.use_numba:
                    cetan, cecross, cint_radius = get_params(dxx, dyy,
                                                             self.data['sch1'],
                                                             self.data['sch2'])
                    etan.extend(cetan)
                    ecross.extend(cecross)
                    int_radius.extend(cint_radius)
                elif self.use_numexpr:
                    dxxs = numexpr.evaluate("dxx**2")
                    dyys = numexpr.evaluate("dyy**2")
                    square_radius = numexpr.evaluate("dxxs + dyys")
                    cos2phi = numexpr.evaluate("(dxxs - dyys) / square_radius")
                    sin2phi = numexpr.evaluate("2.0 * dxx * dyy / square_radius")
                    # Compute rotated ellipticities
                    sch1 = self.data['sch1']
                    sch2 = self.data['sch2']
                    etan.extend(numexpr.evaluate("- (sch1 * cos2phi + sch2 * sin2phi)"))
                    ecross.extend(numexpr.evaluate("- (sch2 * cos2phi - sch1 * sin2phi)"))
                    # Transform cube of distances into a cube of integers
                    # (will serve as indexes for weight)
                    int_radius.extend(np.array(numexpr.evaluate("sqrt(square_radius)"), dtype=int))
                else:
                    dxxs, dyys = dxx**2, dyy**2
                    square_radius = dxxs + dyys
                    cos2phi = (dxxs - dyys) / square_radius
                    sin2phi = 2.0 * dxx * dyy / square_radius
                    # Compute rotated ellipticities
                    etan.extend(- (self.data['sch1'] *
                                   cos2phi + self.data['sch2'] * sin2phi))
                    ecross.extend(- (self.data['sch2'] *
                                     cos2phi - self.data['sch1'] * sin2phi))
                    # Transform cube of distances into a cube of integers
                    # (will serve as indexes for weight)
                    int_radius.extend(
                        np.array(np.sqrt(square_radius), dtype=int))
                pbar.update(ii + jj + (len(xarange) - 2) * ii + 1)

            # Apply this indexes filter to get the right weigth array for each pixel of the grid
            # The output is also a 3D array (nxpoints, nypoints, len(x))
            weights = {cmap: self.weights[cmap][int_radius]
                       for cmap in self.weights}
            for cmap in self.weights:
                if cmap not in self.maps:
                    self.maps[cmap] = []
                sumw = np.sum(weights[cmap], axis=1)
                cell = ecross if '45' in cmap else etan
                self.maps[cmap].append(
                    np.sum(weights[cmap] * cell, axis=1) / sumw)
        pbar.finish()

    def plot_maps(self, clust_coord=None, wcs=None, figsize=(10, 12)):
        """Plot the "kappa" maps."""
        if not hasattr(self, 'maps'):
            raise IOError("WARNING: No maps computed yet.")
        print("Plotting the following maps:")
        for cmap in self.maps.keys():
            print(" - ", cmap)
            fig = pl.figure(figsize=figsize)
            ax = fig.add_subplot(111)
            ax.set_xlabel(xlabel='X-coord (pixel)', fontsize=16)
            ax.set_ylabel(ylabel='Y-coord (pixel)', fontsize=16)
            extent = (min(self.data['xsrc']) + 0.5, max(self.data['xsrc']) + 0.5,
                      min(self.data['ysrc']) + 0.5, max(self.data['ysrc']) + 0.5)
            themap = ax.imshow(self.maps[cmap], origin='lower', zorder=0,
                               cmap=pl.cm.afmhot, extent=extent)
            cb = fig.colorbar(themap, pad=0.15 if wcs is not None else 0.05)
            cb.set_label(cmap, fontsize=20)
            cb.ax.tick_params(labelsize=14)
            ax.scatter(self.data['xsrc'] - 0.5, self.data['ysrc'] - 0.5,
                       s=3, color='b', zorder=1, alpha=0.4)
            if clust_coord is not None:
                if clust_coord[0] >= extent[0] and clust_coord[0] <= extent[1] and \
                   clust_coord[1] >= extent[2] and clust_coord[1] <= extent[3]:
                    ax.plot(clust_coord[0], clust_coord[1],
                            color='g', ms=25, mew=4, marker='+')
                elif wcs is None:
                    print("WARNING: You must provide coordinates in degree + the wcs" +
                          " or coordinates in pixel units.")
                    print(
                        "         Use wcs = kappa.load_wcs('data.hdf5') to get the wcs.")
                else:
                    xsrc, ysrc = cutils.skycoord_to_pixel(clust_coord, wcs)
                    ax.plot(xsrc, ysrc, color='g', ms=25, mew=4, marker='+')
            ax.set_xlim(xmin=min(self.data['xsrc']),
                        xmax=max(self.data['xsrc']) + 1)
            ax.set_ylim(ymin=min(self.data['ysrc']),
                        ymax=max(self.data['ysrc']) + 1)
            ax.tick_params(axis='both', labelsize=14)
            if wcs is not None:
                ra = cutils.pixel_to_skycoord(ax.get_xticks(),
                                              [np.mean(self.data['ysrc'])] *
                                              len(ax.get_xticks()),
                                              wcs).ra.value
                dec = cutils.pixel_to_skycoord([np.mean(self.data['xsrc'])] * len(ax.get_yticks()),
                                               ax.get_yticks(),
                                               wcs).dec.value
                ax2 = ax.twiny()
                ax2.set_xlim(ax.get_xlim())
                ax2.set_xticks(ax.get_xticks())
                ax2.set_xticklabels(["%.2f" % r for r in ra], fontsize=14)
                ax2.set_xlabel("RA (deg)", fontsize=16)
                ax2 = ax.twinx()
                ax2.set_ylim(ax.get_ylim())
                ax2.set_yticks(ax.get_yticks())
                ax2.set_yticklabels(["%.2f" % r for r in dec], fontsize=14)
                ax2.set_ylabel("DEC (deg)", fontsize=16)
            fig.savefig(cmap + ".png")
        pl.show()

    def plot_quiver(self, clust_coord=None, wcs=None, figsize=(10, 12)):
        """Quiver plot"""
        if not hasattr(self, 'maps'):
            raise IOError("WARNING: No maps computed yet.")
        print("Quiver plot for the following maps:")
        for cmap in self.maps.keys():
            if cmap.endswith("45"):
                continue
            print(" - ", cmap)
            fig = pl.figure(figsize=figsize)
            ax = fig.add_subplot(111)
            ax.set_xlabel(xlabel='X-coord (pixel)', fontsize=16)
            ax.set_ylabel(ylabel='Y-coord (pixel)', fontsize=16)

            def rebin2(inArray, shape):
                # rebin of 2D array into shape
                sh = shape[0], inArray.shape[0] // shape[0], shape[1], inArray.shape[1] // shape[1]
                return inArray.reshape(sh).mean(-1).mean(1)
            # rebin2(np.array(self.maps[cmap]),[64,64])
            gamma1 = np.array(self.maps[cmap])
            # rebin2(np.array(self.maps[cmap+'45']),[64,64])
            gamma2 = np.array(self.maps[cmap + '45'])
            gamma_mag = (gamma1**2 + gamma2**2)**0.5
            gamma_phi = 0.5 * np.arctan2(gamma2, gamma1)

            gammax = gamma_mag * np.cos(gamma_phi)
            gammay = gamma_mag * np.sin(gamma_phi)
            ax.quiver(gammax, gammay)
            ax.set_title(cmap)
            # extent = (min(self.data['xsrc']) + 0.5, max(self.data['xsrc']) + 0.5,
            #          min(self.data['ysrc']) + 0.5, max(self.data['ysrc']) + 0.5)
            # themap = ax.imshow(self.maps[cmap], origin='lower', zorder=0,
            #                   cmap=pl.cm.afmhot, extent=extent)
            #cb = fig.colorbar(themap, pad=0.15 if wcs is not None else 0.05)
            #cb.set_label(cmap, fontsize=20)
            # cb.ax.tick_params(labelsize=14)
            # ax.scatter(self.data['xsrc'] - 0.5, self.data['ysrc'] - 0.5,
            #           s=3, color='b', zorder=1, alpha=0.4)
            # if clust_coord is not None:
            #    if clust_coord[0] >= extent[0] and clust_coord[0] <= extent[1] and \
            #       clust_coord[1] >= extent[2] and clust_coord[1] <= extent[3]:
            #        ax.plot(clust_coord[0], clust_coord[1], color='g', ms=25, mew=4, marker='+')
            #    elif wcs is None:
            #        print("WARNING: You must provide coordinates in degree + the wcs" +
            #              " or coordinates in pixel units.")
            #        print("         Use wcs = kappa.load_wcs('data.hdf5') to get the wcs.")
            #    else:
            #        xsrc, ysrc = cutils.skycoord_to_pixel(clust_coord, wcs)
            #        ax.plot(xsrc, ysrc, color='g', ms=25, mew=4, marker='+')
            #ax.set_xlim(xmin=min(self.data['xsrc']), xmax=max(self.data['xsrc']) + 1)
            #ax.set_ylim(ymin=min(self.data['ysrc']), ymax=max(self.data['ysrc']) + 1)
            #ax.tick_params(axis='both', labelsize=14)
            # if wcs is not None:
            #    ra = cutils.pixel_to_skycoord(ax.get_xticks(),
            #                                 [np.mean(self.data['ysrc'])] * len(ax.get_xticks()),
            #                                 wcs).ra.value
            #    dec = cutils.pixel_to_skycoord([np.mean(self.data['xsrc'])] * len(ax.get_yticks()),
            #                                  ax.get_yticks(),
            #                                  wcs).dec.value
            #    ax2 = ax.twiny()
            #    ax2.set_xlim(ax.get_xlim())
            #    ax2.set_xticks(ax.get_xticks())
            #    ax2.set_xticklabels(["%.2f" % r for r in ra], fontsize=14)
            #    ax2.set_xlabel("RA (deg)", fontsize=16)
            #    ax2 = ax.twinx()
            #    ax2.set_ylim(ax.get_ylim())
            #    ax2.set_yticks(ax.get_yticks())
            #    ax2.set_yticklabels(["%.2f" % r for r in dec], fontsize=14)
            #    ax2.set_ylabel("DEC (deg)", fontsize=16)
            #fig.savefig(cmap + ".png")
        pl.show()

    def save_maps(self):
        """Save the maps in a fits files."""
        # Now write the files out as fits files
        # To modify header, use this syntax: hdu.header["cd1_1"]=XXXX
        for cmap in self.maps:
            hdu = pyfits.PrimaryHDU(self.maps[cmap])
            hdu.writeto("%s.fits" % cmap, overwrite=True) # clobber=True is deprecated


to_jit = []
to_vectorize = []


def numba_jit(func):
    """Apply jit decorator."""
    to_jit.append(func)
    return func


def numba_vectorize(func):
    """Applu vectorize decorator."""
    to_vectorize.append(func)
    return func


def apply_decorator(func, decorator):
    """Apply a given decorator to a function."""
    globals()[func.func_name] = decorator(func)


@numba_jit
def get_params(dx, dy, sch1, sch2):
    """Compute tangential and cross ellipticities and filter selection for weight array."""
    dxs, dys = squared_array(dx), squared_array(dy)
    square_radius = sum_arrays(dxs, dys)
    cos2phi = compute_cos2phi(dxs, dys, square_radius)
    sin2phi = compute_sin2phi(dx, dy, square_radius)
    # Compute rotated ellipticities
    etan = compute_etan(sch1, sch2, cos2phi, sin2phi)
    ecross = compute_ecross(sch1, sch2, cos2phi, sin2phi)
    # Transform cube of distances into a cube of integers
    # (will serve as indexes for weight)
    int_radius = np.array(sqrt_array(square_radius), dtype=int)
    return etan, ecross, int_radius


@numba_vectorize
def squared_array(x):
    """Squared an array."""
    return x**2


@numba_vectorize
def sqrt_array(x):
    """Squared root of an array."""
    return x**(1. / 2)


@numba_vectorize
def sum_arrays(x, y):
    """Sum two arrays."""
    return x + y


@numba_vectorize
def compute_cos2phi(dxs, dys, square_radius):
    """Compuet cos 2 phi."""
    return (dxs - dys) / square_radius


@numba_vectorize
def compute_sin2phi(dx, dy, square_radius):
    """Compute sin 2 phi."""
    return 2.0 * dx * dy / square_radius


@numba_vectorize
def compute_etan(sch1, sch2, cos2phi, sin2phi):
    """Compute tangential ellipticity."""
    return - (sch1 * cos2phi + sch2 * sin2phi)


@numba_vectorize
def compute_ecross(sch1, sch2, cos2phi, sin2phi):
    """Compute cross ellipticity."""
    return - (sch2 * cos2phi - sch1 * sin2phi)


def load_data(datafile):
    """Load the needed deepCoadd_meas catalog."""
    return cutils.read_hdf5(datafile, path='deepCoadd_meas', dic=False)


def load_wcs(datafile):
    """Load the WCS."""
    return cutils.load_wcs(cutils.read_hdf5(datafile, 'wcs', dic=False))


def get_cat(datafile, **kwargs):
    """Get a clean catalog (mostly for test purpose here).

    kwargs list can contain:

    :param string ell1: default is 'ext_shapeHSM_HsmShapeRegauss_e1'
    :param string ell2: default 'ext_shapeHSM_HsmShapeRegauss_e2'
    """
    ell1 = kwargs.get('ell1', 'ext_shapeHSM_HsmShapeRegauss_e1')
    ell2 = kwargs.get('ell1', 'ext_shapeHSM_HsmShapeRegauss_e2')
    cat = load_data(datafile)
    return cat[(abs(cat[ell1]) < 1.2) & (abs(cat[ell2] < 1.2) & (cat['filter'] == 'i'))]


def aperture_mass_maturi_filter(tanhx, **kwargs):
    """Maturi et al filter for aperture mass."""
    tanha = kwargs.get("tanha", 6.0)
    tanhb = kwargs.get('tanhb', 150.0)
    tanhc = kwargs.get('tanhc', 50.0)
    tanhd = kwargs.get('tanhd', 47.0)
    tanhxc = kwargs.get('tanhxc', 0.1)
    return np.tanh(tanhx / tanhxc) / ((tanhx / tanhxc) *
                                      (1 + np.exp(tanha - tanhb * tanhx) +
                                       np.exp(tanhc * tanhx - tanhd)))
