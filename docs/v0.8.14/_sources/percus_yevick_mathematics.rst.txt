.. _percus-yevick-mathematics:

Percus--Yevick: mathematics and assumptions
===========================================

Percus--Yevick (PY) is a closure approximation for the correlations of a
homogeneous, isotropic equilibrium hard-sphere mixture. It produces functions
of separation and wavenumber rather than particle centres or a deposition
history. See :doc:`analytical` for the Python workflow,
:doc:`metropolis_mathematics` for explicit equilibrium sampling, and
:doc:`rsa_mathematics` for irreversible deposition.

Physical inputs and correlation functions
-----------------------------------------

Species :math:`i` has radius :math:`a_i>0`, number fraction :math:`x_i`, and
number density :math:`\rho_i`, with :math:`\sum_i x_i=1`. A **number fraction**
counts particles; a **volume fraction** weights their occupied volume. For
total volume fraction :math:`\eta`,

.. math::

   \rho=\frac{\eta}{\sum_i x_i(4\pi a_i^3/3)},\qquad
   \rho_i=x_i\rho,\qquad
   \eta_i=\rho_i\frac{4\pi a_i^3}{3},\qquad \sum_i\eta_i=\eta.

**Additive** spheres have contact distance :math:`\sigma_{ij}=a_i+a_j`.
Their potential is infinite below contact and zero outside. The equilibrium
fluid assumptions exclude walls, directional order, attractions, and
irreversible insertion. Mathematical validity of an input volume fraction
does not establish a stable fluid phase at that density.

The **partial pair correlation** :math:`g_{ij}(r)` gives the relative
likelihood of finding a species :math:`j` centre at separation :math:`r`
from species :math:`i`:

.. math::

   \mathrm{d}N_j=\rho_jg_{ij}(r)4\pi r^2\mathrm{d}r,\qquad
   h_{ij}(r)=g_{ij}(r)-1.

:math:`h` is the **total correlation**: it includes both direct and indirect
correlation through other particles. In a homogeneous fluid, :math:`g` tends
to one far away; the hard core requires :math:`g=0` below contact.

Ornstein--Zernike equation and the PY postulate
-----------------------------------------------

The **direct correlation** :math:`c_{ij}` separates the direct part from
chains mediated by other species. The mixture Ornstein--Zernike equation is

.. math::

   h_{ij}(\mathbf{r})=c_{ij}(\mathbf{r})+
   \sum_\ell\rho_\ell\int_{\mathbb{R}^3}
   c_{i\ell}(\mathbf{r}-\mathbf{r}')h_{\ell j}(\mathbf{r}')\,
   \mathrm{d}^3\mathbf{r}'.

This relation alone contains two unknown functions. A **closure** supplies
another relation; PY postulates

.. math::

   c_{ij}(r)=\bigl(1+h_{ij}(r)-c_{ij}(r)\bigr)
   \bigl(\exp[-\beta u_{ij}(r)]-1\bigr).

For hard spheres this is equivalent to

.. math::

   h_{ij}(r)=-1\quad(r<\sigma_{ij}),\qquad
   c_{ij}(r)=0\quad(r>\sigma_{ij}).

The first condition enforces impenetrability; the second is the approximation.
The total correlation outside contact generally remains nonzero because the
Ornstein--Zernike convolution propagates indirect correlations. A closure is
not a deposition or particle-motion rule.

Fourier transforms and PackLab's matrix normalization
-----------------------------------------------------

**Wavenumber** :math:`k` has units of inverse length. For an isotropic function,
PackLab's three-dimensional radial transform convention is

.. math::

   \widetilde{h}_{ij}(k)=4\pi\int_0^\infty
   r^2h_{ij}(r)\operatorname{sinc}(kr)\,\mathrm{d}r,\qquad
   \operatorname{sinc}(z)=\frac{\sin z}{z},\quad\operatorname{sinc}(0)=1.

Here :math:`\widetilde h` has units of volume. PackLab returns the
**density-normalized**, dimensionless matrices

.. math::

   C_{ij}(k)=\sqrt{\rho_i\rho_j}\,\widetilde c_{ij}(k),\qquad
   H_{ij}(k)=\sqrt{\rho_i\rho_j}\,\widetilde h_{ij}(k).

The result field ``Cpy`` is :math:`C`, and ``H`` is :math:`H`, not the
unnormalized transform :math:`\widetilde h`. Multiplying the Fourier-space
Ornstein--Zernike equation by the density factors gives

.. math::

   H=C+CH,\qquad (I-C)H=C,\qquad S=I+H=(I-C)^{-1}.

:math:`I` is the identity matrix and :math:`S` is the corresponding partial
**structure-factor matrix**. It describes density fluctuations in reciprocal
space. ``H`` and ``Cpy`` have shape ``(species, species, wavenumber)``;
``g`` and ``h`` use distance as their last axis. Form :math:`S` from ``H``
by adding the species identity at each wavenumber. An experimentally weighted
mixture scattering intensity also needs species amplitudes and their weights;
it is not simply every element of :math:`S` summed together.

The inverse transform for positive-density species is

.. math::

   h_{ij}(r)=\frac{1}{2\pi^2\sqrt{\rho_i\rho_j}}
   \int_0^\infty k^2H_{ij}(k)\operatorname{sinc}(kr)\,\mathrm{d}k.

PackLab evaluates this numerically on a finite wavenumber grid. It explicitly
enforces the hard core in the reported ``g``; ``h`` retains the numerical
inverse transform, so ``g = 1 + h`` can fail inside the core. For pairs
involving a zero-density species, the inverse transform sets ``h`` to zero
and ``g`` is zero inside contact and one outside; this convention does not
compute infinite-dilution tracer correlations.

Resolution, periodic terminology, and interpretation
----------------------------------------------------

Use ``wavenumber="auto"`` for a first calculation. For a chosen real-space
resolution :math:`\delta r` and largest requested distance :math:`r_{\max}`,
a useful grid design is

.. math::

   k_{\max}=\frac{\pi}{\delta r},\qquad
   \Delta k\lesssim\frac{2\pi}{m r_{\max}},

where :math:`m` is the number of samples per oscillation of the inverse kernel
(default 12). Increasing :math:`k_{\max}` resolves sharper features;
decreasing :math:`\Delta k` resolves oscillations over longer distances.
These are independent requirements. Finite cutoff causes ringing near the
contact discontinuity; coarse spacing causes integration errors. Retain
PackLab's inadequate-sampling warning and test convergence by increasing the
cutoff and refining the spacing separately. A smooth-looking curve alone is
not an accuracy test. Compare distances away from contact as well as the
contact region, where limits from either side differ.

**Periodic images** and **minimum-image distances** apply to the explicit
simulation box (see :doc:`rsa_mathematics`). The PY solver assumes an infinite
bulk mixture; its wavenumber grid is an integration grid, not a list of
periodic box modes. It has no wrapped centres, winding numbers, or finite
particle count. Match radii, densities, and normalization when comparing with
Metropolis, and check equilibration and box size on the simulation side.
RSA need not agree with PY because its deposition process is irreversible.
Grid refinement reduces numerical error within PY; it does not remove the
physical error of the closure approximation.

Further reading
---------------

Percus and Yevick introduced the approximation in
`Analysis of Classical Statistical Mechanics by Means of Collective Coordinates
<https://doi.org/10.1103/PhysRev.110.1>`_ (1958). For an example of mixture
structure-factor normalization and scattering weights, see
`Measurable structure factors of dense dispersions containing polydisperse
optically inhomogeneous particles <https://journals.iucr.org/j/issues/2024/05/00/ge5154/>`_
(2024). The equations above specify the normalization used by PackLab.
