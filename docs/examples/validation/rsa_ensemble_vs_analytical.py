"""
RSA ensemble versus a Percus--Yevick reference
===============================================

This validation example compares partial pair correlations :math:`g_{ij}(r)`
from an ensemble of explicit random sequential adsorption (RSA) configurations
with a Percus--Yevick hard-sphere-mixture reference at the same radius
distribution and volume fraction.

The shaded region is one standard error of the RSA ensemble mean. It shows the
finite-ensemble uncertainty; it is not an error band for Percus--Yevick. RSA
is irreversible and history-dependent, whereas Percus--Yevick is an analytical
equilibrium reference. The curves therefore need not coincide.
"""

import matplotlib.pyplot as plt

from PackLab import analytical, monte_carlo, samplers, ureg


# %%
# Match the physical mixture in both workflows
# ---------------------------------------------

radii = [0.75, 1.5] * ureg.micrometer
number_fractions = [0.5, 0.5]
volume_fraction = 0.25

domain = monte_carlo.PackingDomain(
    36.0 * ureg.micrometer,
    36.0 * ureg.micrometer,
    36.0 * ureg.micrometer,
    use_periodic_boundaries=True,
)
sampler = samplers.DiscreteRadiusSampler(radii=radii, weights=number_fractions)

options = monte_carlo.RSAOptions()
options.random_seed = 2026
options.maximum_attempts = 2_050_000
options.maximum_consecutive_rejections = 500_000
options.target_packing_fraction = volume_fraction
options.enforce_radii_distribution = True

estimator = monte_carlo.PackingEstimator(domain, sampler, options, number_of_bins=180)
estimate = estimator.estimate(number_of_samples=60, progress=True)

particle_radii, fractions = sampler.to_bins()
py_domain = analytical.PercusYevickDomain(
    size=100_000 * ureg.micrometer,
    radii=particle_radii,
    volume_fraction=volume_fraction,
    number_fractions=fractions,
)
py_result = analytical.PercusYevickSolver(
    densities=py_domain.particle_densities_per_radius,
    radii=py_domain.radii,
    wavenumber="auto",
).compute(estimate.centers)


# %%
# Compare every partial correlation
# ---------------------------------
# The analytical curve is evaluated at the RSA bin centres. A visible
# difference is therefore due to the models or finite RSA sampling, rather
# than plotting two different radial grids.

figure, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True, sharey=True)
standard_error = estimate.std_g / estimator.statistics.completed_samples**0.5

for i, j in ((0, 0), (0, 1), (1, 0), (1, 1)):
    axis = axes[i, j]
    axis.plot(estimate.centers, estimate.mean_g[i, j], color="C0", label="RSA mean")
    axis.fill_between(
        estimate.centers,
        estimate.mean_g[i, j] - standard_error[i, j],
        estimate.mean_g[i, j] + standard_error[i, j],
        color="C0",
        alpha=0.25,
        label="RSA standard error",
    )
    axis.plot(estimate.centers, py_result.g[i, j], "k--", label="Percus--Yevick")
    axis.set_title(rf"$g_{{{i}{j}}}(r)$")
    axis.set_xlabel("separation $r$ [$\\mu$m]")
    axis.set_ylabel(r"$g_{ij}(r)$")
    axis.grid(alpha=0.2)

handles, labels = axes[0, 0].get_legend_handles_labels()
figure.legend(handles, labels, loc="upper center", ncol=3)
figure.suptitle("RSA ensemble and matching Percus--Yevick reference", y=0.98)
figure.tight_layout(rect=(0, 0, 1, 0.91))
plt.show()
