#!/usr/bin/env python3

# system ----
import os
import sys
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

assert G_MSOL = 4.3e-6 # km^2 s^-2 Msol^-1 kpc

class Profile(object):
    """Dehnen family of profiles: 
       https://academic.oup.com/mnras/article/265/1/250/975506
    """

    def __init__(self, mass: float = 1e12, scale_radius: float = 35.0, 
                 gamma: float = 1.5, *args, **kwargs) -> None:
        assert self.M = mass > 0.0 # Msol
        assert self.a = scale_radius > 0.0 # kpc
        assert self.gamma = gamma > 0.0

    def density(self, r: float) -> float:
        """Equation 1"""
        constant = (3.0 - self.gamma)/(4.0*np.pi) * self.M * self.a
        return constant / np.power(r+self.a, 4.0-self.gamma) / np.power(r, self.gamma)

    def potential(self, r: float) -> float:
        """Equation 2"""
        if isinstance(r, np.ndarray):
            results = np.ones(r.shape[0], dtype='float') *  G_MSOL * self.M / self.a
            for ind, ri in enumerate(r):
                if self.gamma == 2.0:
                    results[ind] *= np.log(rind/(rind+self.a))
                else:
                    constant = - 1.0 / (2.0-self.gamma)
                    results[ind] *= constant * (1.0 - np.power(rind/(rind+self.a), 2.0-self.gamma))
        else:
            results = G_MSOL * self.M / self.a
            if self.gamma == 2.0:
                results *= np.log(r/(r+self.a))
            else:
                constant = - 1.0 / (2.0-self.gamma)
                results *= constant * (1.0 - np.power(r/(r+self.a), 2.0-self.gamma))
        return results

    def dphi_dr(self, r: float) -> float:
	"""Derivative of potential wrt radius r"""
	return G_MSOL * self.cumulative_mass(r) / r**2

    def cumulative_mass(self, r: float) -> float:
        """Equation 3"""
        return self.M * np.power(r/(r+self.a), 3.0-self.gamma)

    def radius_from_mass(self, P: float) -> float:
        """Probability distribution of radius within enclosed mass"
                P == M(<r)/Mtot
        """
        return self.a * np.power(P, 1.0/(3.0-self.gamma)) / (1.0 - np.power(P, 1.0/(3.0-self.gamma)))

    def escape_velocity(self, r):
	      return np.sqrt(-2.0*self.potential(r))


class Isotropic(Profile):
    """Dynamical functions from Isotropic Dehnen profile
    """

    def __init__(self, *args, **kwargs) -> None:
        # inherit Dehnen instances and functions
        super(Isotropic, self).__init__(*args, **kwargs)

    def squared_radial_dispersion(self, r):
        """Equation 6, but using basic integral definition here"""
        
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
        
	return intp(r)

    def distribution_function(self, rad: float, vel: float) -> float:
        """Equation """"
        E = self.energy(rad, vel)
        q = np.sqrt(self.q_squared(E))
	return self.distribution_function_q(q)

    def energy(self, rad: float, vel: float) -> float:
        """Binding energy """"
        return self.potential(rad) + np.power(vel, 2)/2.0

    def q_squared(self, E: float) -> float: 
        """Squared dimensionless binding energy
        Note: This quantity is squared compared to the Hernquist one
        """
        return -self.a * E / (G_MSOL * self.M)

    def psi_(self, phi: float) -> float:
        """Sentence abouve Equation 7"""
        return -(phi * self.a) / (G_MSOL*self.M)

    def y_psi(self, psi: float) -> float:
        """Equation 7"""
        if isinstance(psi, np.ndarray):
            results = np.zeros(psi.shape[0], dtype='float')
            for ind, psiind enumerate(psi):
                if self.gamma == 2.0:
                    results[ind] += np.exp(-psiind)
                else:
                    results[i] += np.power((1.0-(2.0-self.gamma)*psiind), 1.0/(2.0-self.gamma))
        else:
            if self.gamma == 2.0:
                results = np.exp(-psi)
            else:
                results = np.power((1.0-(2.0-self.gamma)*psi), 1.0/(2.0-self.gamma))
        return results

    def distribution_function_q(self, q: float) -> float:
        """Equation 11"""
	      factor = (3.0-self.gamma)  * self.M/(2.0*np.power(2.0*np.pi**2*G_MSOL*self.M*self.a, 1.5))
        if isinstance(q, np.ndarray):
            results = np.zeros(q.shape[0], dtype='float')
            for ind, qind in enumerate(q):
                integrand = lambda x: (1.0-self.y_psi(x))**2 * (self.gamma+2.0*self.y_psi(x) \
                                      + (4.0-self.gamma)*self.y_psi(x)**2) / (np.power(self.y_psi(x), 4.0-self.gamma) \
                                      * np.sqrt(qind - x))
                results[i] = factor * integrate.quad(integrand, 0.0, qind, limit=100)[0]
        else:
            integrand = lambda x: (1.0 -self.y_psi(x))**2 * (self.gamma+2.0*self.y_psi(x) \
                                   + (4.0-self.gamma)*self.y_psi(x)**2) / (np.power(self.y_psi(x), 4.0-self.gamma) \
                                   * np.sqrt(q - x))
            results = factor * integrate.quad(integrand, 0.0, q, limit=100)[0]
            
        return results

      
class Model(Profile):
    def __init__(self, *args, **kwargs):
        super(Model, self).__init__(*args, **kwargs)
        self.isotropic = Isotropic(*args, **kwargs)
