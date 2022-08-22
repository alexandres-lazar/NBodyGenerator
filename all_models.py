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
                 gamma: float = 1.5, distr: str ='isotropic', *args, **kwargs) -> None:
        self.hernquist = hernquist.Model(concentration=concentration,
                             scale_radius=scale_radius,
                             *args, **kwargs)
        self.dehnen = dehnen.Model(concentration=concentration,
                             scale_radius=scale_radius,
                             gamma=gamma,
                             *args, **kwargs)
        self.nfw = nfw.Model(concentration=concentration,
                             scale_radius=scale_radius,
                             *args, **kwargs)
