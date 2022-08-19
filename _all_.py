#!/usr/bin/env python3

# system ----
import os
import sys
import numpy as np
import scipy.integrate as integrate
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

# local ----
import hernquist
import dehnen

class Profiles(object):

    def __init__(self, distribution: str ='isotropic', *args, **kwargs) -> None:
        self.hernquist = hernquist.Model( *args, **kwargs)
        self.dehnen = dehnen.Model(*args, **kwargs)
