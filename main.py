import os
import sys
import numpy as np

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

# CONSTANTS AND CONVERISONS
# -------------------------
G = 4.3e-9 * 1e+3  # (km/s)^2 Msol^-1 kpc

class main():
    #def __init__(self, **kwargs):
    #    self.__dict__.update(kwargs)
    def __init__(self, M=1e+12, a=35., N=1e+6):
        self.M = float(M) # Total mass of halo in solar mass
        self.a = float(a) # Scale radius of profile in kpc
        self.N = int(N)   # Total number of particles to generate

    def generate_nbody(self, *args, **kwargs):
        a = self.a
        N = self.N

        print "Generating the positions and velocities for your profile..."
        pos = np.zeros((N,3))
        vel = np.zeros((N,3))
        for _ in range(N):
            # Scheme of sampling cartesian positions from
            # the cumulative mass profile
            x = np.random.uniform()
            if(x > 0.):
                # ***
                # This needs to grab the inverse of the
                # profiles cumulative mass function
                r = a * np.sqrt(x)/(1.-np.sqrt(x))
            elif(x == 0.):
                r = 0.
            theta = np.random.uniform() * (2.*np.pi)
            varphi = np.random.uniform() * np.pi
            pos_x = r * np.sin(theta) * np.cos(varphi)
            pos_y = r * np.sin(theta) * np.sin(varphi)
            pos_z = r * np.cos(theta)
            pos_vec =  [pos_x, pos_y, pos_z]
            pos[_] += pos_vec
            pos_mag = np.sqrt(pos_vec[0]*pos_vec[0] + pos_vec[1]*pos_vec[1] + pos_vec[2]*pos_vec[2])

            # ***
            # This needs to grab the potential from whatever
            # profile selected
            pot = self.potential(pos_mag)
            vmax = np.sqrt(-2. * pot)
            # ***
            # This needs to grab the analytical form
            # of the profiles distribution function
            fmax = vmax * vmax * self.distribution_function(rad, vel=0.0)
            # ***
            # MBK said this was wrong -- should fix
            # Rejection Method
            counter = 0
            fail = 0
            while(counter < 1):
                wf = np.random.uniform() * fmax
                wv = np.random.uniform() * vmax
                if(wf <= (wv * wv * self.distribution_function(rad, vel=wv))):
                    vel_mag = wv
                    counter += 1
                else:
                    failure += 1
            #v_mag = self.generate_vel(r_mag)
            theta = np.random.uniform() * 2*np.pi
            varphi = np.random.uniform() * np.pi
            vel_x = vel_mag * np.sin(theta) * np.cos(varphi)
            vel_y = vel_mag * np.sin(theta) * np.sin(varphi)
            vel_z = vel_mag * np.cos(theta)
            vel_vec = [vel_x, vel_y, vel_z]
            vel[_] += vel_vec

        print "Generation complete, now saving the positions and velocities into separate files..."
        # write scheme to save the positions into a .txt file and the velocities into another .txt file
        return [pos, vel]
