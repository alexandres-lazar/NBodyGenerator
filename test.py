#!/usr/bin/env python3

# system ----
from timeit import default_timer as timer

# local ----
import phase_space

start = timer()
phase_space.Generate(n_particles=1e+6).all()
print(f"Wall-clock time from execuation: {timer() - start} sec")
