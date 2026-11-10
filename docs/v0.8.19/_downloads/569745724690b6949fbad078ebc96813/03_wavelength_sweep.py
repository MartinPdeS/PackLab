"""
Single-sphere wavelength sweep
==============================

Compare scattering cross sections of a fixed sphere across wavelengths in the
same medium. This independent-particle optical calculation contains no
packing structure or multiple scattering.
"""

import matplotlib.pyplot as plt
import numpy as np
from PackLab import scattering, ureg

wavelengths = np.array([400, 450, 500, 550, 600, 650, 700])
cross_sections = []
for wavelength in wavelengths:
    dataset = scattering.compute_scattering_amplitudes(
        wavelength=wavelength * ureg.nanometer,
        diameters=np.array([200]) * ureg.nanometer,
        material=1.59,
        medium=1.33,
        phi=np.linspace(-90, 90, 181) * ureg.degree,
    )
    dataset.process()
    cross_sections.append(dataset.Csca.to("nanometer**2").magnitude[0])
figure, axis = plt.subplots()
axis.plot(wavelengths, cross_sections, "o-")
axis.set(
    xlabel="wavelength (nm)", ylabel="scattering cross section (nm²)", title="Fixed 200 nm sphere"
)
axis.grid(alpha=0.2)
figure.tight_layout()
