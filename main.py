#!/usr/bin/env python3

# system ---
from timeit import default_timer as timer

# local ----
import phase_space

# ---------------------------------------------------------------------------

def main() -> None:
    start = timer()
    print("-" * 75)
    halo_mass = 1e12  # Msol
    npart = 1e6
    c = 12
    #prof_list = ['nfw', 'herquist']
    prof_list = ['nfw']

    for prof_name in prof_list:
        phase_space.Generate(mass=halo_mass,
                             concentration=c,
                             n_particles=npart,
                             profile=prof_name).all()
        print("-" * 75)

    print(f"Wall-clock time from execuation: {timer() - start} sec")
    print("-" * 75)

# ---------------------------------------------------------------------------

if __name__ == '__main__':
    main()
