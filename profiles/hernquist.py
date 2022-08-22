#!/usr/bin/env python3

# system ----
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

G_MSOL = 4.3e-6 # km^2 s^-2 Msol^-1 kpc

from cosmology import SetCosmology

class Profile(SetCosmology):
    """Hernquist profile:
    https://adsabs.harvard.edu/full/1990ApJ...356..359H
    """
    def __init__(self, mass: float = 1e12, concentration: float = None,
                       scale_radius: float = None, *args, **kwargs) -> None:
        super().__init__()
        self.M = mass           # Msol
        self.rs = scale_radius  # kpc
        self.c = concentration
        if self.c is None and self.rs is None:
            raise Exception("!!! Either the concentration or the scale radius must be defined !!!")
        if self.c is not None and self.rs is not None:
            raise Exception("!!! Either define the concentration or the scale radius!!!")

    def halo_radius(self):
        """Halo radius from specified definition. Will only focus on z=0
        with Bryan and Norman (1998) definition.
        """
        return SetCosmology().virial_radius(0.0, self.M)

    def density(self, r: float) -> float:
        return self.M*self.rs / (2.0*np.pi) / r / (r+self.rs)**3

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
        """Equation 5"""
        return - (G_MSOL * self.M)/(r + self.rs)

    def dphi_dr(self, r: float) -> float:
        """Derivative of potential wrt radius r"""
        return G_MSOL * self.cumulative_mass(r) / r**2

    def cumulative_mass(self, r: float) -> float:
        """Equation 3"""
        return self.M * r**2 / (r+self.rs)**2

    def radius_from_mass(self, P: float) -> float:
        """Probability distribution of radius within enclosed mass
                P == M(<r)/Mtot
        """
        return (self.rs*np.sqrt(P)) / (1.0-np.sqrt(P))

    def escape_velocity(self, r: float) -> float:
        return np.sqrt(-2.0 * self.potential(r))


class Isotropic(Profile):
    """Dynamical functions from Isotropic Hernquist profile
    """
    def __init__(self, *args, **kwargs):
        super(Isotropic, self).__init__(*args, **kwargs)

    def squared_radial_dispersion(self, r: float, integrator: str = None) -> float:

        rmax = np.log10(1e3) # can be changed to whatever
        power_arr = np.linspace(1e-3, 3, 50)
        rad_arr = 10.0 ** power_arr

        if integrator == "quadrature":
            # define integrand
            integrand = lambda x: self.density(x) * self.dphi_dr(x)
            # quad integrate integrand within radial range
            from scipy.integrate import quad
            int_arr = np.array([integrate.quad(integrand, r, np.inf, limit=100)[0] for r in rad_arr])
            int_arr /= self.density(rad_arr)

            # create interpolation function to sample `r` from:
            from scipy.interpolate import interp1d
            intp = interp1d(rad_arr, int_arr, fill_value="extrapolate")
            result = intp(r)

        elif integrator == "trapz":
            pass

        elif integrator is None:
            factor = G_MSOL * self.M  / 12.0 / self.rs
            A = 12.0 * r * (r+self.rs)**3/self.rs**4 * np.log((r+self.rs)/r)
            B = - r / (r+self.rs)
            C = 25.0 + 52.0*(r/self.rs) + 42.0*(r/self.rs)**2 + 12.0*(r/self.rs)**3
            result = factor * (A - B*C)

        return result

    def distribution_function(self, rad: float, vel: float) -> float:
        """Equation 17 as function of radius and velocity"""
        E = self.energy(rad, vel)
        q = np.sqrt(self.q_squared(E))
        return self.distribution_function_q(q)

    def energy(self, rad: float, vel: float) -> float:
        """Binding energy"""
        return self.potential(rad) + vel**2/2.0

    def density_of_states(self, rad: float, vel: float) -> float:
        """Equation 23 as function of radius and velocity"""
        E = self.energy(rad, vel)
        q = np.sqrt(self.q_squared(E))
        return self.density_of_states_q(q)

    def q_squared(self, E: float) -> float:
        """Squared dimensionless binding energy
        Note: This quantity is squared compared to the Hernquist one
        """
        return -self.rs * E / (G_MSOL * self.M)

    def distribution_function_q(self, q: float) -> float:
        """Equatin 17"""
        constant = self.M / (self.rs*np.pi)**3 / 4.0 / np.power(2.0*G_MSOL*self.M/self.rs, 1.5)
        A = 1.0 / np.power(1.0-q**2, 2.5)
        B = 3.0 * np.arcsin(q)
        C = q * np.sqrt(1.0 - q**2) * (1.0 - 2.0*q**2) * (8.0*q**4 - 8.0*q**2 - 3.0)
        return constant * A*(B + C)

    def density_of_states_q(self, q: float) -> float:
        """Equation 23"""
        constant = (2.0/3.0 * np.pi**2 * self.rs**2.5 * np.sqrt(2.0*G_MSOL*self.M)) / q**5
        A = 3.0 * (8.0*q**4 - 4.0*q**2 + 1.0) * np.arccos(q)
        B = q * np.sqrt(1.0-q**2) * (4.0*q**2 - 1.0) * (2.0*q**2 + 3.0)
        return constant * (A - B)


class Model(Profile):
    def __init__(self, **kwargs):
        super(Model, self).__init__(**kwargs)
        self.isotropic = Isotropic(**kwargs)

