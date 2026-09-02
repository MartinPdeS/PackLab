Validation and comparison
=========================

These examples make three complementary validation checks:

* compare an ensemble of explicit, finite RSA packings with a matching
  Percus--Yevick reference, including the RSA standard error;
* compare samples from a fixed-volume Metropolis hard-sphere chain with the
  same Percus--Yevick reference;
* compare PackLab's binary-mixture Percus--Yevick solution with digitised
  Percus--Yevick curves from Tsang et al. (2001);
* verify convergence of the numerical Percus--Yevick inverse transform by
  refining its wavenumber grid.

RSA is an irreversible, history-dependent packing process. Metropolis samples
the equilibrium hard-sphere ensemble only after adequate burn-in and sampling.
Percus--Yevick is an analytical equilibrium reference, so agreement with RSA
is informative but is not an expectation of exact equality.
