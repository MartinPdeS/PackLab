:orphan:

.. _theory:

Theoretical background
======================

PackLab offers three complementary hard-sphere workflows. Percus--Yevick (PY)
is an analytical equilibrium reference; random sequential adsorption (RSA)
constructs explicit, irreversible deposition configurations; and Metropolis
Monte Carlo (MC) samples an equilibrium system through particle moves. They
can share radii, composition, and volume fraction, but they do not represent
the same physical process.

Technique guides
----------------

Use the :doc:`analytical` guide for the PY Python workflow and the
:doc:`monte_carlo` guide for RSA packing generation and Metropolis MC moves.
The :doc:`references` page collects the underlying literature and release
citation information.

.. toctree::
   :hidden:

   analytical
   monte_carlo
   references

Percus--Yevick hard-sphere mixtures
-----------------------------------

The analytical namespace evaluates the multicomponent PY approximation for an
equilibrium hard-sphere mixture. Given particle radii, number fractions, and a
total volume fraction, it derives species densities and computes partial pair
correlations :math:`g_{ij}(r)` and structure factors :math:`S_{ij}(k)`.

PY is useful when a rapid equilibrium reference is needed for parameter sweeps
or a structure-aware scattering calculation without generating sphere centres.
It is an approximation: its physical accuracy depends on density, composition,
and size ratio. Its numerical accuracy also depends on the wavenumber grid
used to recover real-space correlations. These are separate questions. The
automatic grid is a useful default, while explicit grids allow a resolution
study for a particular distance range.

Random sequential adsorption
----------------------------

Random sequential adsorption (RSA) proposes one sphere at a time. A proposal
is retained only when it does not overlap an accepted sphere. Once accepted,
the sphere never moves; a rejected proposal is discarded. Early accepted
spheres can therefore block positions permanently, making the final packing
depend on the random proposal order. RSA produces a physically valid,
history-dependent deposition configuration, not an equilibrium hard-sphere
sample. In PackLab, ``maximum_attempts`` limits rejected proposals and
``target_packing_fraction`` provides an optional earlier stopping condition.

For spheres with radii :math:`r_i` in a domain of volume :math:`V`, the packing
fraction is

.. math::

   \phi = \frac{1}{V}\sum_i \frac{4\pi r_i^3}{3}.

The radius sampler is consequently part of the physical model: a
monodisperse simulation and a broad polydisperse simulation at the same
nominal packing fraction need not have the same structure.

Metropolis Monte Carlo
----------------------

Metropolis MC starts from any valid non-overlapping configuration, keeps every
particle radius and class label fixed, and proposes symmetric single-particle
displacements. A proposed move is accepted when it remains in the domain and
creates no overlap. Repeated moves target the equilibrium distribution at the
chosen particle count and volume; unlike RSA, accepted particles can move and
rearrange.

The initial states retain memory of their starting configuration, so discard a
``burn-in`` interval before collecting measurements. Successive MC states may
also be similar; their autocorrelation indicates how many sweeps should
separate approximately independent samples. Finally, repeat a calculation for
larger boxes or particle counts when bulk correlations are required, because a
finite simulation box can alter the measured result.

Periodic boundaries
-------------------

With periodic boundaries enabled, opposite faces of the simulation box are
identified. A sphere close to one face can therefore overlap a sphere close
to the opposite face. This reduces boundary artefacts when estimating bulk
properties such as a pair-correlation function. Use non-periodic boundaries
when the edges themselves are part of the system being modelled.

Pair correlations and structure factors
---------------------------------------

RSA and MC configurations can estimate the radial pair-correlation function
:math:`g_{ij}(r)` from accepted centres. The analytical namespace calculates
the corresponding PY pair correlations and structure factors over a
user-supplied wavenumber grid. Match radii, number fractions, total volume
fraction, boundary conditions, and the definition of each observable before
comparing results.

Agreement between an equilibrated MC calculation and PY can test how closely
the approximation represents a specified equilibrium mixture. RSA--PY
differences are expected even with matched inputs, because irreversible RSA
deposition and an equilibrium hard-sphere fluid are different models. Use
explicit configurations when their centres, finite size, or deposition history
are the quantities of interest; use PY when a fast analytical equilibrium
reference is the useful object.
