Assumptions and limitations
===========================

PackLab provides complementary models, not interchangeable ways to obtain the
same result. Record the physical inputs and choose the workflow whose
assumptions match the question.

Common inputs
-------------

Radii, number fractions, volume fraction, domain dimensions, and boundary
conditions define the physical system. Use explicit units for every dimensional
input. A radius distribution is not merely a numerical setting: it changes the
mixture composition and its structure.

RSA
---

Random sequential adsorption produces non-overlapping configurations by
irreversible deposition. Accepted spheres never move, so early random choices
can block later placements. RSA is appropriate for deposition-like packings and
for workflows that require explicit sphere centres. It is not an equilibrium
hard-sphere sampler.

Check the stopping criterion, random seed, realised packing fraction, and box
size. A configuration can be valid yet still be too small or too sparse for a
stable bulk correlation estimate.

Metropolis Monte Carlo
----------------------

Metropolis moves target an equilibrium hard-sphere system at fixed particle
count and volume. A valid RSA packing is a convenient initial state, but its
deposition history must be forgotten through burn-in before collecting data.
Assess stationarity, autocorrelation, acceptance rate, and sensitivity to box
size; a fixed sweep count alone is not evidence of equilibration.

Percus--Yevick
--------------

The analytical workflow is a Percus--Yevick approximation for an equilibrium
hard-sphere mixture. It supplies pair correlations and structure factors
quickly, but does not create centres, represent deposition history, or remove
finite-size effects from a simulation. Verify numerical resolution by refining
the wavenumber grid when real-space correlations are important.

Scattering
----------

Scattering is optional and uses PyMieSim for the single-particle optical
amplitudes. A mixture phase function additionally depends on the supplied
structure tensor. When that tensor comes from Percus--Yevick, the scattering
result is structure-aware but remains an analytical equilibrium prediction,
not a scattering calculation from an RSA configuration.

Comparisons
-----------

When comparing workflows, match radii, composition, total volume fraction,
boundary conditions, and the definition and resolution of the observable.
Agreement between equilibrated Metropolis samples and Percus--Yevick tests an
analytical approximation. Differences between RSA and Percus--Yevick are often
expected because they describe different physical processes.
