"""
Comparing radius samplers
=========================

Radius samplers are explicit physical inputs to an RSA simulation. Their
``to_bins`` method converts a continuous distribution into radius classes and
number fractions for use by the analytical Percus--Yevick model.
"""

import matplotlib.pyplot as plt
import numpy as np

from PackLab import samplers
from PackLab.units import ureg


uniform = samplers.UniformRadiusSampler(
    minimum_radius=80 * ureg.nanometer,
    maximum_radius=160 * ureg.nanometer,
    bins=8,
)
log_normal = samplers.LogNormalRadiusSampler(
    median_radius=110 * ureg.nanometer,
    geometric_standard_deviation=1.25,
    maximum_radius_clip=250 * ureg.nanometer,
    bins=8,
)

figure, axes = plt.subplots(1, 2, figsize=(9, 3.5), sharey=True)
for axis, sampler, label in zip(axes, (uniform, log_normal), ("Uniform", "Log-normal")):
    _ = sampler.plot_histogram(ax=axis, color="tab:blue", width=7)
    axis.set_title(label)
figure.tight_layout()

assert np.isclose(sum(uniform.to_bins()[1]), 1.0)
assert np.isclose(sum(log_normal.to_bins()[1]), 1.0)
