###############################
# Utilities to calculate NFW Halo Properties
###############################

import copy
import numpy as np
import scipy.optimize
from scipy.integrate import quad

#############################

std_G = 4.3e-9 # Newton's const   in Mpc (km/s)^2 M_sol^{-1
v_c = 299792.458 #km/s

#############################

class ComovingDistMemoization(object):

    def __init__(self, cosmology, memotable=None):

        if memotable is None:
            memotable = {}

        self.memotable = memotable
        self.cosmology = cosmology

    def __call__(self, z):


        if z in self.memotable:
            return self.memotable[z]

        def integrand(z):

            return 1./np.sqrt(self.cosmology.hubble2(z))

        y, err = quad(integrand, 0, z)

        dist = self.cosmology.v_c * y  #to get proper units, ie to put in the hubble length

        self.memotable[z] = dist

        return dist

##############################

class Cosmology(object):

    def __init__(self, omega_m=0.3, omega_l=0.7, h=0.7, w=-1, omega_r=0., G=std_G):

        self.omega_m = omega_m
        self.omega_l = omega_l
        self.omega_r = omega_r
        self.h = h
        self.w = w
        self.G = G
        self.v_c = v_c

        self.comovingdist = ComovingDistMemoization(self)

    def __copy__(self):

        return Cosmology(omega_m=self.omega_m, omega_l=self.omega_l,
                         h=self.h, w=self.w, omega_r=self.omega_r, G=self.G)

    def get_H0(self):
        return 100 * self.h  #km/s/MPC

    H0 = property(get_H0)

    def hubble2(self, z):

        inv_a = 1.+z
        return (self.omega_r * inv_a**4 + self.omega_m * inv_a**3 +
                self.omega_l * (inv_a**(3 * (1+self.w))) +
                (1 - self.omega_m - self.omega_l - self.omega_r) * inv_a**2) * self.H0**2

    def Ez(self, z):

        return np.sqrt(self.hubble2(z))/self.H0

    def get_hubble_length(self):
        return self.v_c / self.H0

    hubble_length = property(get_hubble_length)


    def _integral_1pzOverEz3(self, z_min, z_max=np.inf):

        def integrand(z):
            return (1.0 + z) / (self.Ez(z))**3

        y, err = quad(integrand, z_min, z_max)
        return y

    def UnnormedGrowthFactor(self, z):

        return 5.0 / 2.0 * self.omega_m * self.Ez(z) * self._integral_1pzOverEz3(z)

    def GrowthFactor(self, z):

        return self.UnnormedGrowthFactor(z)/self.UnnormedGrowthFactor(0.)

    def rho_crit(self, z):

        return 3. * self.hubble2(z)/(8 * np.pi * self.G)

    def angulardist(self, z, z2=None):

        if z2 is None:
            return self.comovingdist(z) / (1+z)

        return (self.comovingdist(z2) - self.comovingdist(z)) / (1+z2)

    def beta(self, z, zcluster):

        Ds = np.array([self.angulardist(zi) for zi in z])
        Dls = np.array([self.angulardist(zcluster, zi) for zi in z])

        Dls_over_Ds = np.zeros_like(Dls)
        Dls_over_Ds[Ds > 0] = Dls[Ds > 0] / Ds[Ds > 0]
        Dls_over_Ds[Dls <= 0] = 0

        return Dls_over_Ds


    def beta_s(self, z, zcluster):

        betainf = self.beta([1e6], zcluster)

        beta_s = self.beta(z, zcluster) / betainf

        return beta_s



###################################

class CosmologyFixedException(Exception):
    pass


class CosmologySingleton(object):

    def __init__(self, startCosmology=None, isMutable=True):
        self.isMutable = isMutable


        if startCosmology is None:
            self._cosmology = Cosmology()
        else:
            self._cosmology = copy.copy(startCosmology)

        self.comovingdist = ComovingDistMemoization(self._cosmology)

    def get_cosmology(self):
        return copy.copy(self._cosmology)

    def set_cosmology(self, newcosmo):
        if not self.isMutable:
            raise CosmologyFixedException
        self._cosmology = copy.copy(newcosmo)
        self.comovingdist = ComovingDistMemoization(self._cosmology)

    cosmology = property(get_cosmology, set_cosmology)

    def __getattr__(self, name):
        return getattr(self._cosmology, name)

###################################

std_cosmology = CosmologySingleton(isMutable=False)
    
global_cosmology = CosmologySingleton()


##################################

def deltaC(c, delta=200.):

    c = float(c)

    return (delta/3.) * c**3/(np.log(1+c) - (c/(1+c)))


###################################

def density(r, rs, c, z, cosmology=global_cosmology):

    x = r/rs

    delta_c = deltaC(c)

    rho_c = cosmology.rho_crit(z)

    return delta_c * rho_c/(x * ((1+x)**2))


##################################

def delta_vir(z, cosmology=global_cosmology):


    #Bryan & Norman 1998

    omega_z = cosmology.omega_m * cosmology.rho_crit(0.) * ((1+z)**3)/cosmology.rho_crit(z)

    x = omega_z - 1

    d_v = 18 * np.pi**2 - 82 * x - 39 * x**2 

    return d_v


###################################

def massInsideR(rs, c, z, R, cosmology=global_cosmology):

    x = R/rs
    delta_c = deltaC(c)
    rho_c = cosmology.rho_crit(z)

    return (np.log(1+x) - (x/(1+x))) * 4 * np.pi * delta_c * rho_c * rs**3

##############

def RsMassInsideR(mass, c, z, R, cosmology=global_cosmology):

    delta_c = deltaC(c)
    rho_c = cosmology.rho_crit(z)

    def f(x):

        xp = R/x

        return (np.log(1+xp) - (xp/(1+xp))) * 4 * np.pi * delta_c * rho_c * x**3 - mass

    try:
        rs = scipy.optimize.brenth(f, 0.01, 10.)

    except ValueError, e:
        print '!!!!!!!!!!!'
        print mass, c, f(0.01), f(10.)
        raise e

    return rs

###################################

def massInsideR_amp(rs, amp, z, R, cosmology=global_cosmology):

    Dl = cosmology.angulardist(z)

    betainf = cosmology.beta([1e6], z)[0]

    sigma_c_4pi = cosmology.v_c**2 / (cosmology.G * Dl * betainf)

    deltac_rhoc_rs3_4pi = sigma_c_4pi * amp * rs**2

    x = R / rs

    return (np.log(1+x) - (x/(1+x))) * deltac_rhoc_rs3_4pi

###################################

def rdelta2rs(rdelta, c, delta):

    delta_c = deltaC(c)

    # x = r_delta / rs
    def f(x):
        return 3 * delta_c * (np.log(1+x) - (x/(1+x)))/x**3 - delta

    x0 = scipy.optimize.brenth(f, 0.1, 20)

    return rdelta / x0


###################################

def rdelta(rs, c, delta):

    delta_c = deltaC(c)

    # x = r_delta / rs
    def f(x):
        return 3 * delta_c * (np.log(1 + x) - (x / (1 + x))) / x**3 - delta

    x0 = scipy.optimize.brenth(f, 0.1, 20)

    return x0 * rs

####################################

def Mdelta(rs, c, z, delta, cosmology=global_cosmology):

    r_delta = rdelta(rs, c, delta)

    rho_c = cosmology.rho_crit(z)

    return delta * rho_c * (4 * np.pi/3) * r_delta**3


######################################

def rdeltaConstM(mdelta,z, delta, cosmology=global_cosmology):

    rho_c = cosmology.rho_crit(z)

    rdelta = (3 * mdelta/(4 * delta * np.pi * rho_c))**(1./3.)

    return rdelta

#######################################


def rscaleConstM(mdelta, c200, z, delta, cosmology=global_cosmology):

    rho_c = cosmology.rho_crit(z)

    rdelta = (3 * mdelta/(4 * delta * np.pi * rho_c))**(1./3.)

    if delta == 200.:
        return rdelta / c200

    return rdelta2rs(rdelta, c200, delta)



###############################






