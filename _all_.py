#!/usr/bin/env python3

# system ----
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

# local ----
import hernquist
import dehnen

class Profiles(object):

    def __init__(self, distr: str ='isotropic', *args, **kwargs) -> None:
        self.hernquist = hernquist.Model( *args, **kwargs)
        self.dehnen = dehnen.Model(*args, **kwargs)
