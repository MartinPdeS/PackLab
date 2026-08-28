:orphan:

.. _theory:

Theoretical background
======================

PackLab offers complementary descriptions of sphere packings. The Monte-Carlo
tools generate explicit configurations, while the analytical tools evaluate a
Percus--Yevick approximation. They answer related, but not identical,
questions and are often most useful together.

Random sequential adsorption
----------------------------

Random sequential adsorption (RSA) proposes one sphere at a time. A proposal
is retained only when it does not overlap an accepted sphere. The process
therefore produces a physically valid, history-dependent configuration rather
than an equilibrium packing. In PackLab, ``maximum_attempts`` limits rejected
proposals and ``target_packing_fraction`` provides an optional earlier stopping
condition.

For spheres with radii :math:`r_i` in a domain of volume :math:`V`, the packing
fraction is

.. math::

   \phi = \frac{1}{V}\sum_i \frac{4\pi r_i^3}{3}.

The radius sampler is consequently part of the physical model: a
monodisperse simulation and a broad polydisperse simulation at the same
nominal packing fraction need not have the same structure.

Periodic boundaries
-------------------

With periodic boundaries enabled, opposite faces of the simulation box are
identified. A sphere close to one face can therefore overlap a sphere close
to the opposite face. This reduces boundary artefacts when estimating bulk
properties such as a pair-correlation function. Use non-periodic boundaries
when the edges themselves are part of the system being modelled.

Pair correlations and structure factors
---------------------------------------

The Monte-Carlo result can estimate the radial pair-correlation function
:math:`g(r)` from its accepted centres. The analytical namespace calculates a
Percus--Yevick solution and its associated structure factor over a user-supplied
wave-vector grid. The analytical model is efficient for parameter sweeps;
explicit RSA configurations are preferable when finite-size effects, a chosen
radius distribution, or the actual sphere centres matter.
