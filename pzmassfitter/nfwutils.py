###############################
# Utilities to calculate NFW Halo Properties
###############################

from numpy import *
from scipy.integrate import quad
import numpy as np
import varcontainer as vc
import scipy.optimize

#############################

std_G=4.3e-9 # Newton's const   in Mpc (km/s)^2 M_sol^{-1
v_c = 299792.458 #km/s

class Cosmology(object):

    def __init__(self, omega_m = 0.3, omega_l = 0.7, h = 0.7, omega_r = 0., G = std_G):

        self.omega_m = omega_m
        self.omega_l = omega_l
        self.omega_r = omega_r
        self.h = h
        self.G = G


    

    def get_H0(self):
        return 100*self.h  #km/s/MPC

    H0 = property(get_H0)

    def hubble2(self, z):

        inv_a = 1.+z
        return (self.omega_r*inv_a**4 + self.omega_m*inv_a**3 + \
                  self.omega_l + (1 - self.omega_m - self.omega_l - self.omega_r)*inv_a**2)*self.H0**2

    def get_hubble_length(self):
        return v_c / self.H0

    hubble_length = property(get_hubble_length)
    

    def rho_crit(self, z):

        return 3.*self.hubble2(z)/(8*np.pi*self.G)


    

std_cosmology = Cosmology()



##################################

def deltaC(c):

    c = float(c)

    return (200./3.)*c**3/(np.log(1+c) - (c/(1+c)))


###################################

def density(r, rs, c, z, cosmology = std_cosmology):

    x = r/rs

    delta_c = deltaC(c)

    rho_c = cosmology.rho_crit(z)

    return delta_c*rho_c/(x*((1+x)**2))

###################################

def massInsideR(rs, c, z, R, cosmology = std_cosmology):

    x = R/rs
    delta_c = deltaC(c)
    rho_c = cosmology.rho_crit(z)

    return (np.log(1+x) - (x/(1+x)))*4*np.pi*delta_c*rho_c*rs**3

##############

def RsMassInsideR(mass, c, z, R, cosmology = std_cosmology):

    delta_c = deltaC(c)
    rho_c = cosmology.rho_crit(z)

    def f(x):

        xp = R/x

        return (np.log(1+xp) - (xp/(1+xp)))*4*np.pi*delta_c*rho_c*x**3 - mass

    try:
    
        rs = scipy.optimize.brenth(f, 0.01, 10.)

    except ValueError, e:
        print '!!!!!!!!!!!'
        print mass, c, f(0.01), f(10.)
        raise e

    return rs
    

###################################

def massInsideR_amp(rs, amp, z, R, cosmology = std_cosmology):

    comovingdist = ComovingDistMemoization(cosmology)

    Dl = angulardist(z, comovingdist = comovingdist)
    
    betainf = beta([1e6], z, comovingdist = comovingdist)[0]

    sigma_c_4pi = v_c**2 / (cosmology.G * Dl * betainf)

    deltac_rhoc_rs3_4pi = sigma_c_4pi*amp*rs**2

    x = R / rs

    return (np.log(1+x) - (x/(1+x)))*deltac_rhoc_rs3_4pi

###################################

def rdelta2rs(rdelta, c, delta, cosmology = std_cosmology):

    delta_c = deltaC(c)
    
    # x = r_delta / rs
    def f(x):
        
        return 3*delta_c*(np.log(1+x) - (x/(1+x)))/x**3 - delta

    
    x0 = scipy.optimize.brenth(f, 0.1, 20)

    return rdelta / x0


###################################

def rdelta(rs, c, delta, cosmology = std_cosmology):

    delta_c = deltaC(c)
    
    # x = r_delta / rs
    def f(x):
        
        return 3*delta_c*(np.log(1+x) - (x/(1+x)))/x**3 - delta

    
    x0 = scipy.optimize.brenth(f, 0.1, 20)

    return x0*rs

####################################

def Mdelta(rs, c, z, delta, cosmology = std_cosmology):

    r_delta = rdelta(rs, c, delta, cosmology)

    rho_c = cosmology.rho_crit(z)

    return delta*rho_c*(4*np.pi/3)*r_delta**3


######################################

def rscaleConstM(mdelta, c, z, delta, cosmology = std_cosmology):

    rho_c = cosmology.rho_crit(z)

    rdelta = (3*mdelta/(4*delta*np.pi*rho_c))**(1./3.)

    return rdelta2rs(rdelta, c, delta, cosmology)



###############################



class ComovingDistMemoization(object):

    def __init__(self, cosmology = std_cosmology, memotable = None):

        if memotable is None:
            memotable = {}

        self.memotable = memotable
        self.cosmology = cosmology
    def __call__(self, z):

        if z in self.memotable:
            return self.memotable[z]

        def integrand(z):

            return 1./(sqrt(self.cosmology.omega_m*(1+z)**(3) + self.cosmology.omega_l))

        y, err = quad(integrand, 0, z)
    
        dist = self.cosmology.hubble_length * y

        self.memotable[z] = dist

        return dist


comovingdist = ComovingDistMemoization()        


def angulardist(z, z2 = None, comovingdist = comovingdist):

    if z2 is None:
        return comovingdist(z) / (1+z)

    return (comovingdist(z2) - comovingdist(z)) / (1+z2)


def beta_memo(z, zcluster, memotable = {}):
    z_round = around(z,2)
    notInMemoTable = array([curZ not in memotable for curZ in z_round])
    if notInMemoTable.any():
        zs_tocalc = unique(z_round[notInMemoTable])
        new_betas = beta(zs_tocalc, zcluster, calcAverage = False)
        for curZ, curBeta in zip(zs_tocalc, new_betas):
            memotable[curZ] = curBeta

    betas = array([memotable[curZ] for curZ in z_round])

    return betas, memotable



def beta(z, zcluster, comovingdist = comovingdist):

    Ds = array([angulardist(zi, comovingdist = comovingdist) for zi in z])
    Dls = array([angulardist(zcluster, zi, comovingdist = comovingdist) for zi in z])

    Dls_over_Ds = zeros_like(Dls)
    Dls_over_Ds[Ds > 0] = Dls[Ds > 0] / Ds[Ds > 0]
    Dls_over_Ds[Dls <= 0] = 0

    return Dls_over_Ds
    

def beta_s(z, zcluster, comovingdist = comovingdist):

    betainf = beta([1e6], zcluster, comovingdist = comovingdist)

    beta_s = beta(z, zcluster, comovingdist = comovingdist) / betainf

    return beta_s

