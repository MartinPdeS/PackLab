.. _metropolis-mathematics:

Metropolis Monte Carlo: mathematics and assumptions
===================================================

Metropolis sampling moves an existing non-overlapping packing at fixed
particle count, radii, composition, and box volume. Its postulate is the
canonical equilibrium distribution for rigid hard spheres. An RSA packing
can supply the starting centres, but equilibration must be assessed before
interpreting samples as equilibrium measurements. See :doc:`monte_carlo`
for API usage and :doc:`rsa_mathematics` for irreversible deposition.

Hard-sphere potential and equilibrium measure
---------------------------------------------

For species radii :math:`a_i,a_j`, the additive hard-sphere potential is

.. math::

   u_{ij}(r)=\begin{cases}
   +\infty,&r<a_i+a_j,\\
   0,&r\geq a_i+a_j.
   \end{cases}

For a configuration :math:`\mathbf{X}`, the potential energy :math:`U` sums
these pair interactions and any wall exclusions. Its configurational
probability density is

.. math::

   \pi(\mathbf{X})=Z^{-1}\exp[-\beta U(\mathbf{X})],\qquad
   \beta=(k_BT)^{-1}.

:math:`Z` normalizes probability over allowed configurations. Consequently,
all valid configurations have the same density, and overlapping states have
zero probability. Temperature cancels from this purely hard-sphere acceptance
rule. Packing fraction stays fixed because particle radii and volume do.

Symmetric proposals and acceptance
----------------------------------

PackLab proposes one displacement per particle per **sweep**, in particle
index order. Each component is sampled independently from
:math:`\mathrm{Uniform}[-d,d]`, where :math:`d` is
``maximum_displacement``. Thus :math:`d` is a coordinate bound, not a bound on
the vector length, whose maximum is :math:`\sqrt{3}d`.

A **symmetric proposal** has equal forward and reverse probability density.
The single-particle Metropolis rule is then

.. math::

   P_{\mathrm{accept}}(\mathbf{X}\to\mathbf{X}')=
   \min\{1,\exp[-\beta(U(\mathbf{X}')-U(\mathbf{X}))]\}.

Starting from a valid packing, this is one if the new centre satisfies walls
and all pair exclusions, and zero otherwise. A rejected move retains the
previous centre. Rejected states must still contribute to time-indexed
measurements; measuring only accepted moves changes the sampled distribution.
Each particle update obeys **detailed balance** (equal equilibrium probability
flow forward and backward); their sequential sweep preserves the equilibrium
measure, although the composed sweep need not itself obey detailed balance.

Reaching equilibrium also requires adequate exploration, or **mixing**.
Local moves can mix very slowly in crowded states. A valid configuration or
high acceptance rate alone does not establish equilibration. Small proposals
are often accepted but explore slowly; large proposals can be mostly rejected.

Periodic geometry and accumulated displacement
----------------------------------------------

A **periodic image** is a particle translated by integer box lengths. The
**minimum-image convention** finds the nearest such copy for overlap checks;
its formula and a worked example are in :doc:`rsa_mathematics`. Spheres must
satisfy :math:`2a\leq\min(L_x,L_y,L_z)` to avoid their own periodic images.
Stored centres are wrapped into the primary box after accepted periodic moves.

Wrapped centres are suitable for pair distances, but lose the number of times
a particle has crossed a face. A **winding** counts these traversals.
**Unwrapped displacement** retains every accepted proposal vector:

.. math::

   \mathbf{D}_n(s)=\sum_{\ell\leq s}
   I_{n\ell}\,\delta\mathbf{x}_{n\ell},\qquad
   \mathrm{MSD}(s;s_0)=\frac{1}{N}\sum_n
   \left\|\mathbf{D}_n(s)-\mathbf{D}_n(s_0)\right\|^2.

:math:`I_{n\ell}` is one for an accepted move and zero for a rejection.
``simulator.cumulative_displacements`` stores :math:`\mathbf{D}_n` as a length
quantity of shape ``(N, 3)`` and clears it on ``reset``.
``run_metropolis_diagnostics`` measures mean squared displacement (MSD)
relative to the start of that call, even when the simulator has already run.
Tracking the actual proposal preserves travel even across several box lengths
in one accepted move. Minimum-image subtraction of two endpoints would lose
this information and artificially bound the MSD.

MSD here measures movement in Monte Carlo sweeps. Metropolis proposals are
numerical sampling steps, not Newtonian or Brownian time evolution. A physical
diffusion coefficient requires a separately justified dynamical model.
Because MSD relative to the beginning of the run typically changes with time,
it is a mobility diagnostic rather than a stationary equilibrium observable.
By default the diagnostic report contains MSD and move statistics, and its
observable and uncertainty fields are ``None``. Supply an ``observable``
callable to request uncertainty estimates for a scalar quantity chosen to be
approximately stationary after burn-in. Choosing an observable does not itself
establish stationarity or equilibration.

Burn-in, autocorrelation, and uncertainty
-----------------------------------------

**Burn-in** is an initial interval excluded from measurements because the
chain still remembers its starting state. **Stationarity** means the measured
observable's distribution no longer changes with sampling time. Check several
observables, different initial states, and longer runs; no fixed burn-in count
guarantees stationarity.

For a stationary scalar observable :math:`A_t`, define the lag correlation
and integrated autocorrelation time (in retained-sample intervals) by

.. math::

   \rho_A(\ell)=\frac{\operatorname{Cov}(A_t,A_{t+\ell})}
   {\operatorname{Var}(A_t)},\qquad
   \tau_{\mathrm{int}}=\frac12+\sum_{\ell=1}^{\infty}\rho_A(\ell),\qquad
   N_{\mathrm{eff}}\simeq\frac{M}{2\tau_{\mathrm{int}}}.

:math:`M` is the number of retained observations, and
:math:`N_{\mathrm{eff}}` is the **effective sample size** for that observable.
It is not the particle count. PackLab truncates the correlation sum at its
first non-positive lag, so estimates can be unreliable for short or complex
trajectories. **Block averaging** divides observations into complete blocks;
if block means are approximately independent, the standard error is their
sample standard deviation divided by the square root of the block count.
PackLab's interval uses a normal quantile; few blocks can give poor coverage.
If every retained value is identical, autocorrelation and uncertainty are
reported as unavailable: the trace alone cannot distinguish a truly constant
quantity from a chain that has failed to explore its distribution.

For example, after constructing a periodic ``simulator``, the cosine component
of its longest-wavelength density mode can be monitored as follows:

.. code-block:: python

   import numpy as np
   from PackLab import monte_carlo

   box_length_m = simulator.domain.length_x.to("meter").magnitude

   def density_mode(configuration):
       x = configuration.positions.to("meter").magnitude[:, 0]
       return float(np.cos(2 * np.pi * x / box_length_m).mean())

   report = monte_carlo.run_metropolis_diagnostics(
       simulator, number_of_sweeps=10_000, sample_interval=10,
       burn_in_sweeps=2_000, observable=density_mode,
       observable_name="density_mode_x", block_size=50,
   )
   print(report.effective_sample_size, report.confidence_interval)

This is an illustration, not an established burn-in or block size: inspect
the observable, compare independent chains, and check block-size sensitivity.
Statistics use only regularly spaced post-burn-in samples; a final shorter
interval remains in the recorded trajectory but is excluded from them.
The report does not prove equilibration.
Periodic boundaries reduce wall effects but preserve finite-box effects;
compare larger boxes and independent seeded trajectories before drawing bulk
conclusions. Equilibrated Metropolis and PY comparisons test an approximation;
RSA and PY can differ for physical reasons even with matched inputs.

Further reading
---------------

The original algorithm and periodic construction appear in Metropolis,
Rosenbluth, Rosenbluth, Teller, and Teller,
`Equation of State Calculations by Fast Computing Machines
<https://doi.org/10.1063/1.1699114>`_ (1953).
