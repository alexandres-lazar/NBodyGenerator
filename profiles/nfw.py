#!/usr/bin/env python3

# system ----
import numpy as np
from scipy.special import spence
from scipy.special import lambertw
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

G_MSOL = 4.3e-6  # km^2 s^-2 Msol^-1 kpc

from cosmology import SetCosmology

class Profile(SetCosmology):
    """The NFW profile:
    https://arxiv.org/abs/astro-ph/9508025
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
        rs = self.halo_radius() / self.c
        x = r / rs
        rhos = self.normalization()
        return rhos / x / (1.0+x)**2

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

    def normalization(self):
        rvir = self.halo_radius()
        if isinstance(self.rs, float):
            rs = self.rs
            c = rvir / rs
        else:
            c = self.c
            rs = rvir / c
        return self.M / (4.0*np.pi*rs**3) / self.fx(c)

    def fx(self, x: float) -> float:
        """Result us factor of mass integration"""
        return np.log(1.0+x) - x/(1.0+x)

    def gx(self, x: float) -> float:
        """Inverse of fx"""
        return 1.0 / self.fx(x)
    """
    def potential_(self, r: float) -> float:
        rvir = self.halo_radius()
        s = r / rvir
        if isinstance(self.rs, float):
            rs = self.rs
        else:
            rs = rvir / self.c
        cs = s * self.c
        factor = -self.M * G_MSOL * self.gx(self.c) / rvir
        return factor * np.log(1.0+cs) / s
    """
    def potential(self, r: float) -> float:
        rvir = self.halo_radius()
        rhos = self.normalization()
        if isinstance(self.rs, float):
            rs = self.rs
        else:
            rs = rvir / self.c
        factor = 4.0 * np.pi * rs**2 * rhos * G_MSOL
        x = r / rs
        return -factor * np.log(1.0+x) / x

    def dphi_dr(self, r: float) -> float:
        """Derivative of potential wrt radius r"""
        return G_MSOL * self.cumulative_mass(r) / r**2

    def cumulative_mass(self, r: float) -> float:
        """Equation 14 from LM00"""
        rvir = self.halo_radius()
        if isinstance(self.rs, float):
            rs = self.rs
        else:
            rs = rvir / self.c
        rhos = self.normalization()
        x = r / rs
        return 4.0 * np.pi * rs**3 * rhos * self.fx(x)

    def escape_velocity(self, r: float) -> float:
        return np.sqrt(-2.0 * self.potential(r))

    def radius_from_mass(self, P: float) -> float:
        """Implement Lambert W function (zeroeth branch) to sample"""
        if self.c is not None:
            c = self.c
        else:
            c = self.halo_radius() / self.rs

        exp_factor = - (P*self.fx(c) + 1.0)
        x = - np.exp(exp_factor)
        w0 = lambertw(x, k=0)
        w0_real = np.real(w0)
        factor = - self.halo_radius() / c
        return factor * (1.0 + 1.0/w0_real)


class Isotropic(Profile):
    """Dynamical functions from Isotropic NFW profile
    """
    def __init__(self, *args, **kwargs):
        super(Isotropic, self).__init__(*args, **kwargs)

    def squared_radial_dispersion(self, r: float) -> float:
        """Equation 14 from LM00"""
        rvir = self.halo_radius()
        vvir_sq = self.M * G_MSOL / rvir  # km/s
        s = r / rvir

        if isinstance(self.c, float):
            c = self.c
        else:
            c = rvir / self.rs
        cs = c * s

        A = 0.5 * self.c**2 * self.gx(self.c) * s * (1.0+cs)**2
        B = np.pi**2 - np.log(cs) - 1.0/cs - 1.0/(1.0+cs)**2 - 6.0/(1.0+cs)
        C = (1.0 + 1.0/cs**2 - 4.0/cs - 2.0/(1.0+cs)) * np.log(1.0+cs)
        D = 3.0 * np.log(1.0+cs)**2 + 6.0*spence(-cs)
        return vvir_sq * A * (B + C + D)

    def distribution_function(self, rad: float, vel: float) -> float:
        """Equation A2 in Widrow 2000 for Model III"""
        fE = self.fancy_E(rad, vel)

        # parameters of the fit
        F0, q = 9.1968e-2, -2.7419
        lambda_ = 5.0/2.0
        p1, p2, p3, p4 = 0.3620, -0.5639, -0.0859, -0.4912
        P = p1*fE + p2*fE**2 + p3*fE**3 + p4*fE**4

        # fitting function
        A = F0 * np.power(fE, 3.0/2.0) * (1.0-fE)**(-lambda_)
        B = (-np.log(fE)/(1.0-fE))**q * np.exp(P)
        return A * B

    def fancy_E(self, rad: float, vel: float) -> float:
        rho0 = self.normalization()
        rvir = self.halo_radius()

        if isinstance(self.rs, float):
            rs = self.rs
        else:
            rs = rvir / self.c

        E = self.energy(rad, vel)
        energy_diff = E - self.potential_infty()
        denom = 4.0 * np.pi * G_MSOL * rho0 * rs**2
        return - energy_diff / denom

    def energy(self, rad: float, vel: float) -> float:
        """Binding energy"""
        return self.potential(rad) + vel**2/2.0

    def potential_infty(self) -> float:
        """Potential out to infinity, from Widrow 2000,
           solve for this algebraicly
        """
        return 0.0


class Model(Profile):
    def __init__(self, **kwargs):
        super(Model, self).__init__(**kwargs)
        self.isotropic = Isotropic(**kwargs)
