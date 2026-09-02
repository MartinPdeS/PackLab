Examples
========

The gallery is organised by the question each workflow answers:

* **Monte Carlo** creates explicit RSA packings and samples fixed-volume
  equilibrium hard-sphere configurations with Metropolis moves.
* **Analytical** evaluates Percus--Yevick mixture models and scattering
  quantities without generating a packing.
* **Validation** compares RSA and Metropolis estimates with their analytical
  counterparts.
* **Benchmarks** records reproducible runtime, scaling, and memory measurements
  for supported workflows.

Start with the minimal RSA example, then use the Metropolis and validation
examples when you need an equilibrium hard-sphere workflow or want to assess
agreement with the analytical reference. Use benchmarks to measure
computational cost for a clearly specified machine and workload; they do not
establish physical agreement between models.
