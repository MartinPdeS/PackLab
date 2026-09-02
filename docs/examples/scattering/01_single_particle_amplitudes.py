"""
Single-particle scattering amplitudes
====================================

Evaluate the far-field amplitudes of a dielectric sphere and inspect its
scattering cross section. This example uses PackLab's optional PyMieSim
integration and does not require a packing configuration.
"""

import matplotlib.pyplot as plt
import numpy as np

from PackLab import scattering, ureg


# %%
# Evaluate amplitude functions for two sphere diameters
# -----------------------------------------------------

diameters = np.array([150, 300]) * ureg.nanometer
dataset = scattering.compute_scattering_amplitudes(
    wavelength=532 * ureg.nanometer,
    diameters=diameters,
    material=1.59,
    medium=1.33,
    phi=np.linspace(-90, 90, 181) * ureg.degree,
)


# %%
# Plot the polarisation-resolved amplitudes for one diameter
# -----------------------------------------------------------

figure = dataset[0].plot()
_ = figure.suptitle("Amplitude functions for a 150 nm sphere")


# %%
# Compare scattering cross sections across the diameter grid
# -----------------------------------------------------------

dataset.process()
figure, axis = plt.subplots(figsize=(6, 4))
_ = axis.plot(
    diameters.to("nanometer").magnitude,
    dataset.Csca.to("nanometer**2").magnitude,
    "o-",
)
_ = axis.set(
    xlabel="sphere diameter (nm)",
    ylabel=r"scattering cross section (nm$^2$)",
    title="Single-particle scattering cross section",
)
axis.grid(alpha=0.25)
figure.tight_layout()
