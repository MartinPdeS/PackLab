"""
Equilibrating an RSA configuration with Metropolis moves
=========================================================

Random sequential addition (RSA) creates a valid hard-sphere packing, but it
is irreversible and does not sample the equilibrium hard-sphere ensemble.
This example uses that RSA configuration only as an initial state for
single-particle Metropolis moves. The Metropolis sampler keeps every radius,
class label, and the particle count fixed; it only changes sphere centres.

For a production equilibrium calculation, check that observables have reached
stationarity and estimate their uncertainty from independent samples. The
acceptance rate printed below helps tune the proposal displacement.
"""

from PackLab import monte_carlo, samplers
from PackLab.units import ureg
import matplotlib.pyplot as plt

# %%
# Create a convenient non-overlapping initial state with RSA
# -----------------------------------------------------------

domain = monte_carlo.PackingDomain(
    5 * ureg.micrometer,
    5 * ureg.micrometer,
    5 * ureg.micrometer,
    use_periodic_boundaries=True,
)
radius_sampler = samplers.ConstantRadiusSampler(150 * ureg.nanometer, bins=1)

rsa_options = monte_carlo.RSAOptions()
rsa_options.random_seed = 12
rsa_options.maximum_attempts = 50_000
rsa_options.target_packing_fraction = 0.08
initial_result = monte_carlo.RSASimulator(domain, radius_sampler, rsa_options).run()


# %%
# Equilibrate with fixed-volume Metropolis sampling
# -------------------------------------------------

metropolis_options = monte_carlo.MetropolisOptions()
metropolis_options.random_seed = 34
metropolis_options.number_of_sweeps = 500
metropolis_options.maximum_displacement = 50 * ureg.nanometer

simulator = monte_carlo.MetropolisSimulator(
    domain,
    initial_result.sphere_configuration,
    metropolis_options,
)
equilibrated_result = simulator.run()

statistics = simulator.statistics
print(f"Metropolis acceptance rate: {statistics.acceptance_rate:.3f}")
print(f"Completed sweeps: {statistics.completed_sweeps}")
print(f"Packing fraction: {equilibrated_result.statistics.packing_fraction_geometry:.3f}")

figure = equilibrated_result.plot_slice_2d(show=False)
_ = figure.suptitle("Hard-sphere configuration after Metropolis equilibration")
plt.show()
