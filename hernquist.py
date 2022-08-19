#!/usr/bin/env python2

# system ----
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.integrate as integrate
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

# local ----

class Profile(object):

    def __init__(self, mass=1e12, scale_radius=35, **kwargs):
        self.M = mass # Msol
        self.a = scale_radius # kpc
        self.G = 4.3e-6 # km^2 s^-2 Msol^-1 kpc

    def density(self, r):
	return self.M/(2.*np.pi) * self.a/(r*np.power(self.a+r,3))

    def potential(self, r):
	return -(self.G * self.M)/(r + self.a)

    def mass_enclosed(self, r):
	return self.M * np.power(r,2)/np.power(r+self.a,2)

    def radius_from_mass(self, P):
        # P == M(<r)/Mtot
        return (self.a*np.sqrt(P))/(1. - np.sqrt(P))

    def escape_velocity(self, r):
	return np.sqrt(-2. * self.potential(r))

class Isotropic(Profile):

    def __init__(self, **kwargs):
        super(Isotropic,self).__init__(**kwargs)

    def squared_radial_dispersion(self, r):
	G = self.G
        if(type(r) is np.ndarray):
            results = np.zeros(r.size, dtype='float')
            for i in range(0, r.size):
                integrand = lambda x: (self.density(x) * self.G * self.mass_enclosed(x))/np.power(x,2)
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

