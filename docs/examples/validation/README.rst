Validation and comparison
=========================

These examples make complementary validation checks:

* compare an ensemble of explicit, finite RSA packings with a matching
  Percus--Yevick reference, including the RSA standard error, and repeat the
  comparison across packing fractions;
* compare samples from a fixed-volume Metropolis hard-sphere chain with the
  same Percus--Yevick reference;
* compare PackLab's binary-mixture Percus--Yevick solution with digitised
  Percus--Yevick curves from Tsang et al. (2001);
* compare PackLab's monodisperse Percus--Yevick solution at $f=0.2$ and
  $f=0.3$ with digitised curves from the same reference;
* verify convergence of the numerical Percus--Yevick inverse transform by
  refining its wavenumber grid across mixture sizes and volume fractions.

RSA is an irreversible, history-dependent packing process. Metropolis samples
the equilibrium hard-sphere ensemble only after adequate burn-in and sampling.
Percus--Yevick is an analytical equilibrium reference, so agreement with RSA
is informative but is not an expectation of exact equality.
