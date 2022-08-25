# NBodyGenerator
#### _Generates the positions and velocities of a N-body (gravitational) halos utilizing the Collisionless Boltzmann Equation._ 

This was inspired by my final project in graduate Gravitational Dynamics class at UT Austin :metal:. For the main methods used, look to my submitted homework in `_docs/methods.pdf`.

* Works with python version >3.6 

* Simply specify the halo mass [Msol] and concenctration (halo radius computed from Bryan and Norman (1998) definition at redshift zero) or the scale radius [kpc] within `main.py` and run with `python test.py` in the command line. Outputs are stored in `results` directory.

* `phase_space.py` contains the computational work to generate particles:
    * Coordinates are randomly sampled from inverted cumulative mass profiles.
    * Velocities are sampled from distribution function using a Monte Carlo Rejection method.

---

Profiles available are:
* [Hernquist Halo](https://adsabs.harvard.edu/full/1990ApJ...356..359H):
   * Only isotropic dynamics.

* [Navarro, Frenk, and White](https://arxiv.org/pdf/astro-ph/9508025.pdf) (NFW) Halo: 
    * Analytical functions worked out in [Lokas and Mamon (2001)](https://arxiv.org/pdf/astro-ph/0002395.pdf).
    * Radial sampling using real value of [Lambert W function](https://en.wikipedia.org/wiki/Lambert_W_function) transformation of inverted mass profile.
    * Distribution functions approximated using Appendix of (Widrow (2000))[https://arxiv.org/pdf/astro-ph/0003302.pdf].
      * Only isotropic profile encoded.

* __[Work in progress]___ Dehnen (+Jaffe and Hernquist) Halo (Isotropic):
    * _Warning_: Quadrature integration is used solve distribution function makes this computaitonally exhausting.

---

