#!/usr/bin/env python2

# system ----
import os
import sys
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

# local ----
import _all_

class PhaseSpace(_all_.Profiles):

    def __init__(self, Nparticles=1e+6, profile='hernquist', distribution='isotropic', **kwargs):
        super(PhaseSpace,self).__init__(**kwargs)
        self.distribution = distribution
        self.N = Nparticles
        self.profile = _all_.Profiles(**kwargs).__getattribute__(profile)

    def all(self):
        N = int(self.N)
        print "Generating the positions and velocities for your profile..."
        pos_arr = np.zeros((N,3))
        vel_arr = np.zeros((N,3))
        for _ in range(0,N):
            posN = self.position()
            pos_arr[_] += posN
            vel_arr[_] += self.velocity(self.magnitude(posN))
        print "Generation complete, now saving the positions and velocities into separate files..."
        #np.savetxt('positions.txt', pos_arr, delimiter=' ')
        #np.savetxt('velocities.txt', vel_arr, delimiter=' ')
        np.save('positions.npy', pos_arr)
        np.save('velocities.npy', vel_arr)

    def position(self):
        x = np.random.uniform()
        if(x > 0.):
            r = self.profile.radius_from_mass(x)
        elif(x == 0.):
            r = 0.
        theta = np.random.uniform() * (2.*np.pi)
        varphi = np.random.uniform() * np.pi
        pos_x = r * np.sin(theta) * np.cos(varphi)
        pos_y = r * np.sin(theta) * np.sin(varphi)
        pos_z = r * np.cos(theta)
        pos_vec = [pos_x, pos_y, pos_z]
        return pos_vec

    def velocity(self, rad): # This needs to be symbiotic with the positions attribute
        vmax = self.profile.escape_velocity(rad)
        fmax = vmax * vmax * self.profile.__getattribute__(self.distribution).df(rad, vel=0.0)
        counter = 0
        failure = 0
        while(counter < 1): # Rejection Method
            wf = np.random.uniform() * fmax
            wv = np.random.uniform() * vmax
            if(wf <= (wv * wv * self.profile.__getattribute__(self.distribution).df(rad, vel=wv))):
                vel_mag = wv
                counter += 1
            else:
                failure += 1
        theta = np.random.uniform() * 2*np.pi
        varphi = np.random.uniform() * np.pi
        vel_x = vel_mag * np.sin(theta) * np.cos(varphi)
        vel_y = vel_mag * np.sin(theta) * np.sin(varphi)
        vel_z = vel_mag * np.cos(theta)
        vel_vec = [vel_x, vel_y, vel_z]
        return vel_vec

    def magnitude(self, vec):
        return np.sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2])

