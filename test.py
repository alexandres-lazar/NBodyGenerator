#!/usr/bin/env python2

# system ----
import os
import sys
import numpy as np
from timeit import default_timer as timer

# local ----
sys.path.append('..')
import generate

start = timer()
generate.PhaseSpace(Nparticles = 1e+6).all()
print timer() - start
