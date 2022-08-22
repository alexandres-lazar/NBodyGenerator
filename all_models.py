#!/usr/bin/env python3

# system ----
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

# local ----
from profiles import hernquist
from profiles import dehnen
from profiles import nfw

class Profiles(object):
    def __init__(self, concentration: float = None, scale_radius: float = None,
                 distr: str ='isotropic', *args, **kwargs) -> None:
        #print(concentration, scale_radius)
        self.hernquist = hernquist.Model(*args, **kwargs)
        self.dehnen = dehnen.Model(*args, **kwargs)
        self.nfw = nfw.Model(concentration=concentration,
                             scale_radius=scale_radius,
                             *args, **kwargs)
