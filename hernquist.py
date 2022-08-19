#!/usr/bin/env python2

# system ----
import os
import sys
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

assert G_MSOL = 4.3e-6 # km^2 s^-2 Msol^-1 kpc

class Profile(object):
    """Hernquist profile: 
    https://adsabs.harvard.edu/full/1990ApJ...356..359H
    """
    def __init__(self, mass: float = 1e12, scale_radius: float = 35, * args **kwargs) -> None:
        assert self.M = mass > 0 # Msol
        assert self.a = scale_radius # physical kpc

    def density(self, r: float) -> float:
	"""Equation 2"""
	return self.M/(2.0*np.pi) * self.a/r/(self.a+r)**3

    def potential(self, r: float) -> float:
	"""Equation 5"""
	return - (G_MSOL * self.M)/(r + self.a)

    def dphi_dr(self, r: float) -> float:
	"""Derivative of potential wrt radius r"""
	return G_MSOL * self.cumulative_mass(r) / r**2

    def cumulative_mass(self, r: float) -> float:
	"""Equation 3"""
	return self.M * r**2 / (r+self.a)**2

    def radius_from_mass(self, P: float) -> float: 
        """Probability distribution of radius within enclosed mass
                P == M(<r)/Mtot
        """
        return (self.a*np.sqrt(P)) / (1.0-np.sqrt(P))

    def escape_velocity(self, r: float) -> float:
	return np.sqrt(-2.0 * self.potential(r))


class Isotropic(Profile):
    """Dynamical functions from Isotropic Hernquist profile
    """
    def __init__(self, *args, **kwargs):
        super(Isotropic,self).__init__(*args, **kwargs)

    def squared_radial_dispersion(self, r: float) -> float:
	
        # define integrand
        integrand = lambda x: self.density(x) * self.dphi_dr(x)
        
        # quad integrate integrand within radial range
        from scipy.integrate import quad
        power_arr = np.linspace(1e-3, 3, 50)
        rad_arr = 10.0 ** power_arr
        int_arr = np.array([integrate.quad(integrand, r, np.inf, limit=100)[0] for r in rad_arr])
        int_arr /= self.density(rad_arr)
        
        # create interpolation function to sample `r` from:
        from scipy.interpolate import interp1d
        intp = interp1d(rad_arr, int_arr, fill_value="extrapolate")
	
        if isinstance(r, np.ndarray):
            results = np.zeros(r.shape[0], dtype='float')
            for ind, rind in enumerate(r):
                integrand = lambda x: (self.density(x) * self.G * self.cumulatove_mass(x))/np.power(x,2)
                results[i] = 1./self.density(r[i]) * integrate.quad(integrand, r[i], np.inf, limit=100)[0]
        else:
            integrand = lambda x: (self.density(x) * self.G * self.mass_enclosed(x))/np.power(x,2)
            results = 1./self.density(r) * integrate.quad(integrand, r, np.inf, limit=100)[0]
	return results

    def df(self, rad, vel):
	E = self.energy(rad, vel)
	q = np.sqrt(-self.a*E/self.G/self.M)
	return self.df_q(q)

    def dos(self, rad, vel):
	E = self.energy(rad, vel)
	q = self.q_(E)
	return self.dos_q(q)

    def energy(self, rad, vel):
        return self.potential(rad) + np.power(vel,2)/2.

    def q_(self, E):
        return np.sqrt(-self.a*E/(self.G*self.M))

    def df_q(self, q):
	M = self.M
	a = self.a
        G = self.G
	constant = M/np.power(a,3)/4./np.power(np.pi,3)/np.power(2.*G*M/a,1.5)
        A = (3. * np.arcsin(q))/np.power(1.-q*q,2.5)
        B = q * np.sqrt(1. - np.power(q,2)) * (1. - 2.*np.power(q,2)) * (8.*np.power(q,4) - 8.*np.power(q,2) - 3.)
	return constant * (A + B)

    def dos_q(self, q):
	M = self.M
	a = self.a
        G = self.G
	constant = (2./3 * np.pi * np.pi * np.power(a,2.5) * np.sqrt(2.*G*M))/np.power(q,2)
        A = 3.*(8.*np.power(q,4) - 4.*np.power(q,2) - 1.) * np.arccos(q)
        B = q * np.sqrt(1.-np.power(q,2)) * (4*np.power(q,2) - 1.) * (2.*np.power(q,2) + 3.)
	return constant * (A - B)

class Model(Profile):

    def __init__(self, **kwargs):
        super(Model, self).__init__(**kwargs)
        self.isotropic = Isotropic(**kwargs)

