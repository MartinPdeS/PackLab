from PackLab.monte_carlo import Domain, Options, Simulator, samplers

from PackLab.binary.interface_estimator import Estimator
from TypedUnit.units import ureg
import matplotlib.pyplot as plt
import numpy as np

domain = Domain(
    length_x=6.0 * ureg.millimeter,
    length_y=6.0 * ureg.millimeter,
    length_z=6.0 * ureg.millimeter,
    use_periodic_boundaries=True
)

domain.scale(10)

radius_sampler = samplers.Discrete(
    radii=[1, 2] * ureg.millimeter,
    weights=[0.5, 0.5],
)

options = Options()
options.random_seed = 123
options.maximum_attempts = 2_500_000
options.maximum_consecutive_rejections = 500_000
options.target_packing_fraction = 0.3
options.minimum_center_separation_addition = 0.0
options.enforce_radii_distribution = True


estimator = Estimator(domain=domain, radius_sampler=radius_sampler, options=options, number_of_bins=100)

estimate_result = estimator.estimate(number_of_samples=40, maximum_pairs=10_000_000)
mean_g_array = np.asarray(estimate_result.mean_g)
std_g_array = np.asarray(estimate_result.std_g)

plt.figure()
for i in range(2):
    for j in range(2):
        mean_g = mean_g_array[i, j, :]
        std_g = std_g_array[i, j, :]

        plt.plot(mean_g, label=f"g_{i+1}{j+1}(r)")
        plt.fill_between(
            np.arange(len(mean_g)),
            mean_g - std_g,
            mean_g + std_g,
            alpha=0.3
        )

plt.show()
