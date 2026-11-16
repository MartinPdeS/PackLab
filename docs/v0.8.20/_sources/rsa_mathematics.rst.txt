.. _rsa-mathematics:

RSA Monte Carlo: mathematics and assumptions
============================================

Random sequential adsorption (RSA), also called random sequential addition,
builds an explicit packing by random insertion. Here, **Monte Carlo** means
that proposals are drawn randomly. It does not imply equilibrium sampling.
The defining RSA postulate is irreversible deposition: an accepted sphere
stays at its centre forever. See :doc:`monte_carlo` for the Python workflow,
:doc:`metropolis_mathematics` for movable spheres, and
:doc:`percus_yevick_mathematics` for the analytical equilibrium reference.

Model and insertion rule
------------------------

A sphere has centre :math:`\mathbf{x}_n` and radius :math:`a_n>0`. The
orthorhombic box has side lengths :math:`L_x,L_y,L_z` and volume
:math:`V=L_xL_yL_z`. Spheres are rigid, additive, and impenetrable: two
surfaces can touch but their interiors cannot intersect. The exclusion
(or contact) distance between centres is :math:`a_n+a_m`.

For a proposed radius :math:`a` and centre :math:`\mathbf{x}`, the acceptance
indicator is

.. math::

   A(\mathbf{x},a)=B(\mathbf{x},a)
   \prod_{n=1}^{N}\Theta\!\left(d(\mathbf{x},\mathbf{x}_n)-a-a_n\right),

where :math:`\Theta(z)=1` for :math:`z\geq0` and zero otherwise.
:math:`B` checks boundary admissibility, and :math:`d` is the appropriate
centre distance. Acceptance adds one sphere; rejection leaves the packing
unchanged. Every acceptance changes the exclusion region for later trials,
so the result depends on proposal history. Gravity, attraction, particle
deformation, and rearrangement are absent.

With walls, :math:`a\leq x_\alpha\leq L_\alpha-a` for each Cartesian
coordinate :math:`\alpha`. PackLab draws wall-bounded centres from this
admissible interval, rather than counting centres outside it as trials.
``containment_padding`` adds an optional wall margin;
``minimum_center_separation_addition`` adds an optional pair clearance.
These settings change the geometric constraints from ideal touching spheres.

Periodic images and the minimum-image convention
------------------------------------------------

**Periodic boundaries** identify opposite faces of a box. A **periodic
image** is a translated copy of a particle, with centre

.. math::

   \mathbf{x}_n+(m_xL_x,m_yL_y,m_zL_z),\qquad m_\alpha\in\mathbb{Z}.

Images describe a repeated configuration; they are not additional independent
particles in the stored box. Wrapping a centre into the primary box uses
:math:`x_\alpha\bmod L_\alpha`.

The **minimum-image convention** selects the closest translated copy when
computing a separation. For an orthorhombic box,

.. math::

   \Delta x_\alpha^{\mathrm{MI}}=
   \Delta x_\alpha-L_\alpha\operatorname{round}
   \!\left(\frac{\Delta x_\alpha}{L_\alpha}\right),\qquad
   d_{\mathrm{MI}}=\left\|\Delta\mathbf{x}^{\mathrm{MI}}\right\|.

For example, centres at :math:`0.1\,\mathrm{m}` and :math:`9.9\,\mathrm{m}`
along an axis of length :math:`10\,\mathrm{m}` are separated by
:math:`0.2\,\mathrm{m}` across the boundary. A particle must also avoid its
**own** images: :math:`2a\leq\min(L_x,L_y,L_z)`. Equality permits contact.
RSA rejects a trial violating this condition; Metropolis initialization and
packing archive validation reject such a configuration.

Radius probabilities and realised composition
---------------------------------------------

The **proposal distribution** describes radii attempted, whereas the
**realised distribution** describes radii accepted. Larger spheres can be
harder to place, so these distributions may differ. With
``enforce_radii_distribution=True``, PackLab retries a drawn radius at new
centres until it succeeds or a stopping condition is reached. This changes
the deposition protocol; it does not create equilibrium samples.

``NormalRadiusSampler`` draws an underlying normal variable
:math:`X\sim\mathcal{N}(\mu,\sigma^2)`, conditions on :math:`X>0`, then
clips at :math:`a_{\max}`. For :math:`\sigma>0`, with :math:`\Phi` the
standard normal cumulative distribution function, its density below the cap
is

.. math::

   f(a)=\frac{\exp[-(a-\mu)^2/(2\sigma^2)]}
   {\sigma\sqrt{2\pi}\,\Phi(\mu/\sigma)},\quad 0<a<a_{\max}.

There is an atom (a finite probability at one radius) at the upper cap of
mass :math:`[1-\Phi((a_{\max}-\mu)/\sigma)]/\Phi(\mu/\sigma)`.
This is **conditioning** at zero and **clipping** at the upper bound, not
truncation at both ends. The mean parameter must be non-negative; when
:math:`\sigma=0`, it must be positive and the radius is
:math:`\min(\mu,a_{\max})`. Optional binning maps draws to bin midpoints;
``to_bins`` integrates the same probabilities, including the upper atom.
The resulting mean and standard deviation generally differ from the input
parameters. Earlier zero-clipped normal sampling could produce zero radii;
positive conditioning changes that distribution and seeded trajectories.

Packing fraction, correlations, and stopping
--------------------------------------------

The occupied volume fraction is

.. math::

   \eta=\frac{1}{V}\sum_{n=1}^{N}\frac{4\pi a_n^3}{3}.

For a homogeneous mixture, the **partial pair correlation** :math:`g_{ij}`
describes the expected number of species :math:`j` centres in a shell around
a species :math:`i` centre:

.. math::

   \mathrm{d}N_j=\rho_jg_{ij}(r)4\pi r^2\mathrm{d}r.

Here :math:`\rho_j` is number density, not volume fraction. PackLab counts
unordered pairs, contributes both ordered class orientations to the histogram,
and accounts for the finite number of same-class pairs. If
:math:`T=N(N-1)/2` pairs are available and :math:`M<T` are sampled uniformly
without replacement, each has inclusion probability :math:`M/T`; multiplying
sampled counts by :math:`T/M` gives an unbiased estimate of full counts.
``maximum_pairs=0`` examines all pairs. ``random_seed`` controls capped pair
sampling and defaults to deterministic zero; this differs from simulation
seeds, where zero requests a random seed. Pair-distance extraction caps the
number examined before filtering by radius, so it can return fewer distances.

The correlation range ends at half the shortest box side, within which a
periodic spherical shell is unambiguous. With walls, bulk shell normalization
has boundary bias; box-size checks are necessary. Sampling noise decreases
with more pairs and independent packings. Reordering particles leaves full
pair counts unchanged; capped estimates agree statistically, not necessarily
for the same seed.

``maximum_attempts`` limits the outer insertion loop;
``maximum_consecutive_rejections``, ``maximum_spheres``, and
``target_packing_fraction`` provide other stopping conditions. When radius
retries are enabled, one outer iteration can make several placement attempts.
A rejection limit is not a proof of **saturation**, meaning no admissible
insertion remains. Report the seed, stopping rule, realised composition,
volume fraction, and box dimensions. RSA is neither close packing nor an
equilibrium fluid calculation.

Further reading
---------------

Zhang and Torquato discuss RSA and rigorous saturation detection in
`Precise algorithm to generate random sequential addition of hard hyperspheres
at saturation <https://doi.org/10.1103/PhysRevE.88.053312>`_ (2013).
PackLab uses finite proposal budgets rather than their saturation-certifying
algorithm.
