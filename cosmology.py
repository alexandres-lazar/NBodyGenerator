#!/usr/bin/env python3

import sys
import numpy as np
from scipy.integrate import quad

if sys.version_info < (3, 6):
    sys.exit("!!! Please use Python 3.6+ to execute script!!!")

G_MSOL = 4.3e-6        # kpc/Msol (km/s)^2
SPEED_OF_LIGHT = 3e+5  # km/s

class SetCosmology(object):

    def __init__(self, h0: float = 0.6774, Om0: float = 0.3089, Ol0: float = 0.6911,
                       Ob0: float = 0.0486, Or0: float = 0.0,
                       *args, **kwargs) -> None:
        self.h0 = float(h0)
        self.Om0 = float(Om0)
        self.Ol0 = float(Ol0)
        self.Ob0 = float(Ob0)
        self.Or0 = float(Or0)

    def scale_factor(self, z: float) -> float:
        return 1.0 / (z+1.0)

    def redshift(self, a: float) -> float:
        return 1.0/a - 1.0

    def H(self, z: float) -> float:
        H0 = 100 * self.h0  # km/s/Mpc
        omegaMz = (self.Om0) * (1.0+z)**3
        omegaLz = self.Ol0
        return H0 * np.sqrt(omegaMz + omegaLz)  # km/s/Mpc

    def Ea(self, a: float) -> float:
        return (self.Om0) * np.power(a, -3) + self.Ol0

    def Ez(self, z: float) -> float:
        return (self.Om0) * (1.0+z)**3 + self.Ol0

    def delta_vir(self, z: float) -> float:
        """The virial overdensity from Bryan & Norman (1998)"""
        a = self.scale_factor(z)
        Omega = self.Om0*(1.0 + z)**3 / (self.Om0*(1.0+z)*3 + self.Ol0)
        x = Omega - 1.0
        return (18.0*np.pi**2 + 82.0*x - 39.0*x**2)

    def rho_c(self, z: float) -> float:
        a = self.scale_factor(z)
        hubble = self.H(z) * 1e-3             # km/s/kpc
        return 3.0*hubble**2 / (8.0*np.pi*G_MSOL)  # Msol/kpc

    def virial_radius(self, z: float, Mvir: float) -> float:
        Mvir = np.float64(Mvir)
        r_cubed = 3.0*Mvir / (4.0*np.pi*self.delta_vir(z)*self.rho_c(z))
        return np.power(r_cubed, 1.0/3.0)  # physical kpc

    def angular_distance(self, z: float) -> float:
        a = self.scale_factor(z)
        H0 = 100 * self.h0  # km/s/Mpc
        func = lambda x: 1.0 / np.power(x, 2) / self.Ea(x)
        return a * SPEED_OF_LIGHT/H0 * quad(func, a, 1.0)[0] * 1000.0  # kpc/arcsec

    def angular_phy_size(self, z: float, theta: float) -> float:
        return self.angular_distance(z) * theta  * 4.84814e-6 # kpc/arcsecond --> kpc
