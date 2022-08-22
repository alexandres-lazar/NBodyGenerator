# NBodyGenerator
Generates the positions and velocities of a N-body (gravitational) halos utilizing the Collisionless Boltzmann Equation. 

---

Profiles available are:
* Hernquist Halo (Isotropic):

* Dehnen (+Jaffe and Hernquist) Halo (Isotropic):

* Navarro, Frenk, and White (NFW) Halo (Isotropic): 
    * Analytical functions worked our here.
    * Radial sampling using real value of Lambert W function transformation of mass profile.
    * Distribution functions approximated using Appendix of Widrow (2000).

---

Simply specify the halo mass and concenctration (halo radius computed from Bryan and Norman (1998) definition at redshift zero) within `main.py` and run. 
