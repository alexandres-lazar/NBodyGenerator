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

class Generate(_all_.Profiles):

    def __init__(self, n_particles: float = 1e+6, profile_type: str = 'hernquist', 
                       distri: str = 'isotropic', *args, **kwargs):
        super(Generate, self).__init__(*args, **kwargs)
        self.distribution = distribution
        assert self.N = n_particles > 0.0
        self.profile = _all_.Profiles(*args, **kwargs).__getattribute__(profile)

    def all(self) -> float:
        N = np.int64(self.N)
        printf f"| Generating {N:0.3e} the positions and velocities"
        pos_arr = np.zeros((N, 3))
        vel_arr = np.zeros((N, 3))
        for __ in range(N):
            posN = self.position()
            pos_arr[__] += posN
            vel_arr[__] += self.velocity(self.magnitude(posN))
        print("| Particle generation complete now)
        print("| Saving the positions and velocities into separate files")
        np.save('positions.npy', pos_arr)
        np.save('velocities.npy', vel_arr)
        
        return None

    def position(self) -> float:
        """Generates cartesian positions from inverse mass profile"""
        x = np.random.uniform()
        if x > 0.0:
            r = self.profile.radius_from_mass(x)
        elif x == 0.0:
            r = 0.0
        theta = np.random.uniform() * (2.0*np.pi)
        varphi = np.random.uniform() * np.pi
        pos_x = r * np.sin(theta) * np.cos(varphi)
        pos_y = r * np.sin(theta) * np.sin(varphi)
        pos_z = r * np.cos(theta)
        pos_vec = [pos_x, pos_y, pos_z]
        return np.array(pos_vec)

    def velocity(self, rad: float) -> float: 
        """Generates velocities from genereate positions (symbiotic)"""      
        vesc = self.profile.escape_velocity(rad)
        df = self.profile.__getattribute__(self.distribution).distribution_function
        fmax = vmax**2 * df(rad, vel=0.0)
              
        # Monte Carlo rejection method
        counter = 0
        failure = 0
        while counter < 1:
            wf = np.random.uniform() * fmax
            wv = np.random.uniform() * vmax
             = wv * wv * df(rad, vel=wv)
            cond = wf <= wv**2 * df(rad, vel=wv)
            if cond:
                vel_mag = wv
                counter += 1
            else:
                failure += 1
        
        # Compute velocities within sphere
        theta = np.random.uniform() * (2.0*np.pi)
        varphi = np.random.uniform() * np.pi
        vel_x = vel_mag * np.sin(theta) * np.cos(varphi)
        vel_y = vel_mag * np.sin(theta) * np.sin(varphi)
        vel_z = vel_mag * np.cos(theta)
        vel_vec = [vel_x, vel_y, vel_z]
        return vel_vec

    def magnitude(self, vec: float) -> float:
        return np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)

