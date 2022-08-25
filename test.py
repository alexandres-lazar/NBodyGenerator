#!/usr/bin/env python3

# system ---
from timeit import default_timer as timer

# local ----
import phase_space

# ---------------------------------------------------------------------------

def main() -> None:
    start = timer()
    print("-" * 75)

    npart = 1e4

    """
    kwargs = {'mass': 1e12, 'n_particles': npart, 'concentration': 12}
    phase_space.Generate(profile='nfw', **kwargs).all()

    kwargs = {'mass': 1e12, 'n_particles': npart, 'scale_radius': 35}
    phase_space.Generate(profile='hernquist', **kwargs).all()
    """

    kwargs = {'mass': 1e12, 'n_particles': npart, 'scale_radius': 35,
              'gamma': 1.5}
    phase_space.Generate(profile='dehnen', **kwargs).all()

    print("-" * 75)

    print(f"Wall-clock time from execuation: {timer() - start} sec")
    print("-" * 75)

# ---------------------------------------------------------------------------

if __name__ == '__main__':
    main()
