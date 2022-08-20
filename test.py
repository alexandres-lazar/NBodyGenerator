#!/usr/bin/env python3

# system ----
from timeit import default_timer as timer

# local ----
import phase_space

start = timer()
print("-" * 75)
phase_space.Generate(n_particles=1e+6).all()
print("-" * 75)
print(f"Wall-clock time from execuation: {timer() - start} sec")
print("-" * 75)
