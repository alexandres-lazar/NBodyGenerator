# NBodyGenerator
_Generates the positions and velocities of a N-body (gravitational) halos utilizing the Collisionless Boltzmann Equation._ 

* Works with python version >3.6 

* Simply specify the halo mass [Msol] and concenctration (halo radius computed from Bryan and Norman (1998) definition at redshift zero) or the scale radius [kpc] within `main.py` and run with `python main.py` in the command line. Outputs are stored in `results` directory.

* `phase_space.py` contains the computational work to generate particles:
    * Coordinates are randomly sampled from inverted cumulative mass profiles.
    * Velocities are sampled from distribution function using a Monte Carlo Rejection method.

---

Profiles available are:
* Hernquist Halo (Isotropic):

* Navarro, Frenk, and White (NFW) Halo (Isotropic): 
    * Analytical functions worked our here.
    * Radial sampling using real value of Lambert W function transformation of mass profile.
    * Distribution functions approximated using Appendix of Widrow (2000).

* __[Work in progress]___ Dehnen (+Jaffe and Hernquist) Halo (Isotropic):
    * _Warning_: Quadrature integration is used solve distribution function makes this computaitonally exhausting.

---

