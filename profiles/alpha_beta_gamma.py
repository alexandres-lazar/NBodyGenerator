#!/usr/bin/env python3

# system ----
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

G_MSOL = 4.3e-6  # km^2 s^-2 Msol^-1 kpc

class Profile(object):
    """The alpha-beta-gamma profile:
    https://adsabs.harvard.edu/full/1990ApJ...356..359H
    """
    def __init__(self, mass: float = 1e12, scale_radius: float = 35,
                 cutoff_radius: float = None, alpha: float = 1.0,
                 beta: float = 3.0, gamma: float = 1.0,
                 *args, **kwargs) -> None:
        self.M = mass           # Msol
        self.rs = scale_radius  # physical kpc
        self.alpha = alpha
        self.beta = beta
        if gamma > 3.0:
            raise Exception("!!! gamma must be less than 3 for the mass \
                            not to diverge at the center !!!")
        else:
            self.gamma = gamma
        if beta > 3.0:
            self.rcut = np.inf
        elif (beta <= 3.0):
            if cutoff_radius is not None:
                self.rcut == cutoff_radius  # physical kpc
            else:
                self.rcut == 10.0 * self.rs  # physicl kpc
                # raise Exception("with beta <= 3, cutoff_radius mus be specified")
        self.rdec = 0.3 * self.rcut  # physical kpc

    def density(self, r: float) -> float:
        if r <= self.rcut:
            result = density_inner_cutoff(rhos, r)
        elif r > self.rcut:
            result = density_outer_cutoff(rhos, r)
        return result

    def density_inner_cutoff(self, r: float) -> float:
        """Profile iff r <= self.rcut)"""
        x = r/self.rs
        denom = x**self.gamma * np.power(1.0+x**self.alpha,
                                         (self.beta-self.gamma)/self.alpha)
        result = self.rhos() / denom
        return result

    def density_outer_cutoff(self, r: float) -> float:
        """Profile iff r > self.rcut)"""
        A = (r/self.rcut)**self.delta()
        B = np.exp(-(r-self.rcut)/self.rdec)
        return self.density_inner_cutoff(self.rcut) * A * B

    def delta(self):
        xcut = self.rcut/self.rs
        return self.rcut/self.rdec - (self.gamma
               + self.beta*xcut**self.alpha) (1.0+xcut**self.alpha)

    def rhos(self):
        return self.M / (4.0*np.pi*self.rs**3) \
              / (self.integral_M() + self.integral_Mcutoff())

    def integral_M(self):
        """   """
        q = self.rcut / self.rs
        integrand = lambda x: np.power(x, 2.0-self.gamma) \
                        / np.power(1.0+x**self.alpha, (self.beta-self.gamma)/self.alpha)
        return quad(integrand, 0.0, q, limit=100)[0]

    def integral_Mcutoff(self):
        """  """
        q = self.rcut / self.rs
        factor_inv = self.rs**3 * q**self.gamma \
                    * np.power(1.0 + q**self.alpha, (self.beta-self.gamma)/self.alpha)
        integrand = lambda x: x**2 (x/self.rcut)**self.delta() \
                              * np.exp(-(x-self.rcut)/self.rdec)
        return quad(integrand, self.rcut, np.inf, limit=100)[0]

    def potential(self, r: float) -> float:
        rmax = 100 * self.rcut
        power_arr = np.linspace(-2, np.log10(rmax), 100)
        rad_arr = 10.0 ** power_arr

        # integrate weighted density profile
        integrand = lambda x: self.cumulative_mass(x) / x**2
        integral = np.array([quad(integrand, r, np.infty, limit=100) for r in rad_arr])

        intp = interp1d(rad_arr, integral, fill_value='extrapolate')
        return -G_MSOL * intp(r)

    def dphi_dr(self, r: float) -> float:
        """Derivative of potential wrt radius r"""
        return G_MSOL * self.cumulative_mass(r) / r**2

    def cumulative_mass(self, r: float) -> float:
        rmax = 100 * self.rcut
        power_arr = np.linspace(-2, np.log10(rmax), 100)
        rad_arr = 10.0 ** power_arr

        # integrate weighted density profile
        integrand = lambda x: x**2 * self.density(x)
        integral = np.array([quad(integrand, 0.0, r, limit=100) for r in rad_arr])

        # compute interpolation
        intp = interp1d(rad_arr, integral, fill_value='extrapolate')
        return 4.0 * np.pi * intp(r)

    def escape_velocity(self, r: float) -> float:
        return np.sqrt(-2.0 * self.potential(r))


