#!/usr/bin/env python3

# system ----
import h5py
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

# local ----
import all_models

class Generate(all_models.Profiles):
    """Generates the initial conditions from a analytical profile
    """
    def __init__(self, mass: float = 1e12, concentration: float = None,
                       scale_radius: float = None, gamma: float = 1.5,
                       n_particles: float = 1e6,
                       profile: str = 'hernquist', distr: str = 'isotropic',
                       *args, **kwargs):
        #print(concentration, scale_radius)
        super(Generate, self).__init__(concentration=concentration,
                                       scale_radius=scale_radius,
                                       *args, **kwargs)
        self.tot_mass = mass # Msol
        self.N = n_particles
        self.part_mass = mass / n_particles  # Msol
        self.distr = distr
        self.profile_name = profile
        self.gamma = gamma
        if profile == 'dehnen':
            self.profile = all_models.Profiles(concentration=concentration,
                                               scale_radius=scale_radius,
                                               gamma=gamma,
                                               *args, **kwargs).__getattribute__(profile)
        else:
            self.profile = all_models.Profiles(concentration=concentration,
                                               scale_radius=scale_radius,
                                               *args, **kwargs).__getattribute__(profile)
        if scale_radius is not None:
            self.rs = scale_radius
            self.conc = self.profile.concentration()
            self.saveName = f"{profile}_{distr}_n{np.log10(n_particles):.0f}_m{np.log10(mass):.0f}_rs{self.rs:.0f}"
        elif concentration is not None:
            self.conc = concentration
            self.rs = self.profile.scale_radius()
            self.saveName = f"{profile}_{distr}_n{np.log10(n_particles):.0f}_m{np.log10(mass):.0f}_c{self.conc:.0f}"
        if profile == "dehnen":
            self.saveName += f"_gamma{self.gamma:.2f}"

    def all(self) -> float:
        N = np.int64(self.N)
        print(f">>> Particle generation for {self.distr} {self.profile_name} halo with...")
        print(f"::: mass {self.tot_mass:.2e} Msol")
        print(f"::: concentration: {self.conc:.2f}")
        print(f"::: scale_radius: {self.rs:.2e} kpc")
        if self.profile_name == "nfw":
            print(f"::: scale density: {self.profile.normalization():.2e} Msol/kpc")
        if self.profile_name == "dehnen":
            print(f"::: gamma slope: {self.gamma:.2f}")
        print(f"| Generating coordinates and velocities for {N:0.3e} particles")
        mass_arr = np.ones(N) * self.part_mass
        pos_arr = np.zeros((N, 3))
        vel_arr = np.zeros((N, 3))
        for __ in range(N):
            posN = self.position()
            pos_arr[__] += posN
            vel_arr[__] += self.velocity(self.magnitude(posN))
        print("| Particle generation complete")
        print(f": Particle masses: {self.part_mass:0.3e} Msol")
        print(f": Coordinate range: {pos_arr.min():0.3f} -- {pos_arr.max():0.3f} kpc")
        print(f": Velocity range: {vel_arr.min():0.3f} -- {vel_arr.max():0.3f} km/s")

        print("| Saving the positions and velocities into separate files")
        outPath = "results"
        with h5py.File(f'{outPath}/{self.saveName}.hdf5', 'w') as h5:
            p1 = h5.create_group('PartType1')
            p1.create_dataset('Masses', data=mass_arr)
            p1.create_dataset('Coordinates', data=pos_arr)
            p1.create_dataset('Velocities', data=vel_arr)
        """
        np.save(f"{outPath}/positions.npy", pos_arr)
        np.save(f"{outPath}/velocities.npy", vel_arr)
        """
        return None

    def position(self) -> float:
        """Generates cartesian positions from inverse mass profile"""
        p = np.random.uniform()
        if p > 0.0:
            r = self.profile.radius_from_mass(p)
        elif p == 0.0:
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
        df = self.profile.__getattribute__(self.distr).distribution_function
        fmax = vesc**2 * df(rad, vel=0.0)

        # Monte Carlo rejection method
        counter = 0
        failure = 0
        while counter < 1:
            wf = np.random.uniform() * fmax
            wv = np.random.uniform() * vesc
            cond = (wf <= wv**2 * df(rad, vel=wv))
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

