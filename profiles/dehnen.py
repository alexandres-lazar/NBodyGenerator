#!/usr/bin/env python3

# system ----
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

G_MSOL = 4.3e-6  # km^2 s^-2 Msol^-1 kpc

from cosmology import SetCosmology

class Profile(SetCosmology):
    """Dehnen family of profiles:
       https://academic.oup.com/mnras/article/265/1/250/975506
    """

    def __init__(self, mass: float = 1e12, scale_radius: float = 35.0,
                 concentration: float = None, gamma: float = 1.5,
                 *args, **kwargs) -> None:
        super().__init__()
        self.M = mass           # Msol
        self.rs = scale_radius  # kpc
        self.c = concentration
        if self.c is None and self.rs is None:
            raise Exception("!!! Either the concentration or the scale radius must be defined !!!")
        if self.c is not None and self.rs is not None:
            raise Exception("!!! Either define the concentration or the scale radius!!!")
        self.gamma = gamma

    def halo_radius(self):
        """Halo radius from specified definition. Will only focus on z=0
        with Bryan and Norman (1998) definition.
        """
        return SetCosmology().virial_radius(0.0, self.M)

    def density(self, r: float) -> float:
        """Equation 1"""
        constant = (3.0 - self.gamma)/(4.0*np.pi) * self.M * self.rs
        return constant / np.power(r+self.rs, 4.0-self.gamma) / np.power(r, self.gamma)

    def scale_radius(self) -> float:
        if isinstance(self.rs, float):
            return self.rs
        elif self.rs is None:
            rvir = self.halo_radius()
            return rvir / self.c

    def concentration(self) -> float:
        if isinstance(self.c, float):
            return self.c
        elif self.c is None:
            rvir = self.halo_radius()
            return rvir / self.rs

    def potential(self, r: float) -> float:
        """Equation 2"""
        factor = G_MSOL * self.M / self.rs
        if isinstance(r, np.ndarray):
            results = np.ones(r.shape[0], dtype='float') * factor
            for ind, ri in enumerate(r):
                if self.gamma == 2.0:
                    results[ind] *= np.log(rind/(rind+self.rs))
                else:
                    constant = - 1.0 / (2.0-self.gamma)
                    results[ind] *= constant * (1.0 - np.power(rind/(rind+self.rs), 2.0-self.gamma))
        else:
            results = factor
            if self.gamma == 2.0:
                results *= np.log(r/(r+self.rs))
            else:
                constant = - 1.0 / (2.0-self.gamma)
                results *= constant * (1.0 - np.power(r/(r+self.rs), 2.0-self.gamma))
        return results

    def dphi_dr(self, r: float) -> float:
        """Derivative of potential wrt radius r"""
        return G_MSOL * self.cumulative_mass(r) / r**2

    def cumulative_mass(self, r: float) -> float:
        """Equation 3"""
        return self.M * np.power(r/(r+self.rs), 3.0-self.gamma)

    def radius_from_mass(self, P: float) -> float:
        """Probability distribution of radius within enclosed mass"
                P == M(<r)/Mtot
        """
        sigma_inv = 1.0 / (3.0-self.gamma)
        prob_factor = P**sigma_inv / (1.0 - P**sigma_inv)
        return self.rs * prob_factor

    def escape_velocity(self, r: float) -> float:
        return np.sqrt(-2.0*self.potential(r))


class Isotropic(Profile):
    """Dynamical functions from Isotropic Dehnen profile
    """
    def __init__(self, *args, **kwargs) -> None:
        # inherit Dehnen instances and functions
        super(Isotropic, self).__init__(*args, **kwargs)

    def squared_radial_dispersion(self, r: float) -> float:
        """Equation 6, but using basic integral definition here"""

        # define integrand
        integrand = lambda x: self.density(x) * self.dphi_dr(x)

        # quad integrate integrand within radial range
        power_arr = np.linspace(1e-3, 3, 50)
        rad_arr = 10.0 ** power_arr
        int_arr = np.array([integrate.quad(integrand, r, np.inf, limit=100)[0] for r in rad_arr])
        int_arr /= self.density(rad_arr)

        # create interpolation function to sample `r` from:
        intp = interp1d(rad_arr, int_arr, fill_value="extrapolate")
        return intp(r)

    def distribution_function(self, rad: float, vel: float) -> float:
        """Equation 11"""
        E = self.energy(rad, vel)
        q = np.sqrt(self.q_squared(E))
        return self.distribution_function_q(q)

    def energy(self, rad: float, vel: float) -> float:
        """Binding energy"""
        return self.potential(rad) + np.power(vel, 2)/2.0

    def q_squared(self, E: float) -> float:
        """Squared dimensionless binding energy
        Note: This quantity is squared compared to the Hernquist one
        """
        return -self.rs * E / (G_MSOL * self.M)

    def psi_(self, phi: float) -> float:
        """Sentence abouve Equation 7"""
        return -(phi * self.rs) / (G_MSOL*self.M)

    def y_psi(self, psi: float) -> float:
        """Equation 7"""
        if isinstance(psi, np.ndarray):
            results = np.zeros(psi.shape[0], dtype='float')
            for ind, psiind in enumerate(psi):
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
        factor = (3.0-self.gamma) * self.M \
                 / (2.0*np.power(2.0*np.pi**2*G_MSOL*self.M*self.rs, 1.5))
        if isinstance(q, np.ndarray):
            results = np.zeros(q.shape[0], dtype='float')
            for ind, qind in enumerate(q):
                integrand = lambda x: (1.0-self.y_psi(x))**2 * (self.gamma+2.0*self.y_psi(x) \
                                      + (4.0-self.gamma)*self.y_psi(x)**2) \
                                      / (np.power(self.y_psi(x), 4.0-self.gamma) \
                                      * np.sqrt(qind - x))
                results[i] = factor * quad(integrand, 0.0, qind, limit=100)[0]
        else:
            integrand = lambda x: (1.0 -self.y_psi(x))**2 * (self.gamma+2.0*self.y_psi(x) \
                                   + (4.0-self.gamma)*self.y_psi(x)**2) \
                                   / (np.power(self.y_psi(x), 4.0-self.gamma) \
                                   * np.sqrt(q - x))
            results = factor * quad(integrand, 0.0, q, limit=100)[0]
        return results


class Model(Profile):
    def __init__(self, *args, **kwargs):
        super(Model, self).__init__(*args, **kwargs)
        self.isotropic = Isotropic(*args, **kwargs)
