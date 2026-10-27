"""
Independent scattering from a size mixture
==========================================

The independent scattering coefficient sums each cross section weighted by
its number density. This approximation omits interference from particle
correlations; see the structure-factor example for a correlated reference.
"""

import matplotlib.pyplot as plt
import numpy as np
from PackLab import scattering, ureg

diameters = np.array([150, 250]) * ureg.nanometer
dataset = scattering.compute_scattering_amplitudes(
    wavelength=532 * ureg.nanometer,
    diameters=diameters,
    material=1.59,
    medium=1.33,
    phi=np.linspace(-90, 90, 181) * ureg.degree,
)
# Keep the total number density fixed while changing composition.
total_density = 1e18 / ureg.meter**3
fractions = np.linspace(0, 1, 11)
coefficients = []
for fraction in fractions:
    densities = np.array([fraction, 1 - fraction]) * total_density
    coefficients.append(dataset.get_mu_independant(densities).to("1 / meter").magnitude)
figure, axis = plt.subplots()
axis.plot(fractions, coefficients, "o-")
axis.set(
    xlabel="number fraction of 150 nm spheres",
    ylabel="scattering coefficient (1/m)",
    title="Independent scattering at fixed total number density",
)
axis.grid(alpha=0.2)
figure.tight_layout()
