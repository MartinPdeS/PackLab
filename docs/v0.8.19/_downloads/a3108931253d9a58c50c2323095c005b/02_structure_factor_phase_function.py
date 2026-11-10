"""
Mixture phase function with a structure factor
==============================================

Combine optical sphere amplitudes with the Percus--Yevick reciprocal-space
correlation tensor of a binary hard-sphere mixture. The resulting phase
function includes the analytical mixture structure; it does not represent an
RSA configuration.
"""

import numpy as np

from PackLab import analytical, scattering, ureg


# %%
# Calculate the Percus--Yevick mixture structure
# -----------------------------------------------

radii = np.array([75, 125]) * ureg.nanometer
domain = analytical.PercusYevickDomain(
    size=100 * ureg.micrometer,
    radii=radii,
    volume_fraction=0.12,
    number_fractions=np.array([0.6, 0.4]),
)
distances = np.linspace(0.0, 2.0, 240) * ureg.micrometer
structure = analytical.PercusYevickSolver(
    densities=domain.particle_densities_per_radius,
    radii=domain.radii,
    wavenumber="auto",
).compute(distances)


# %%
# Evaluate optical amplitudes for the same size classes
# ------------------------------------------------------

dataset = scattering.compute_scattering_amplitudes(
    wavelength=532 * ureg.nanometer,
    diameters=2 * radii,
    material=1.59,
    medium=1.33,
    phi=np.linspace(-90, 90, 181) * ureg.degree,
)
dataset.process()


# %%
# Project the structure-corrected phase function into two dimensions
# -------------------------------------------------------------------

phi, theta, phase_function = dataset.get_phase_function(
    densities=structure.densities,
    H=structure.H,
    wavenumber=structure.wavenumber,
    theta_points=120,
)
phase_function_values = (
    phase_function.magnitude if hasattr(phase_function, "magnitude") else phase_function
)
_ = scattering.plottings.plot_phase_function_2d_projection(
    phi=phi,
    theta=theta,
    phase_function=phase_function_values,
    projection="azimuth_average",
    normalize=True,
)
