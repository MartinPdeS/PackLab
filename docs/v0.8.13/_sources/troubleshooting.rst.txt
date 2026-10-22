Troubleshooting
===============

Installation and imports
------------------------

Use the same Python interpreter for installation and execution::

    python -m pip install --upgrade PackLab
    python -c "import PackLab; print(PackLab.__version__)"

A missing native module often means that no compatible wheel was installed or
a source build did not finish. Confirm your Python version and platform before
changing the simulation. Source builds need the repository's compiler and native
library dependencies.

Plots on a server
-----------------

Use ``MPLBACKEND=Agg`` for Matplotlib on a headless machine. Save figures to files
when no interactive display is available.

Physical assumptions and resolution
-----------------------------------

Use dimensional quantities for lengths. For analytical calculations, start with
``wavenumber="auto"`` and retain Fourier-sampling warnings. For RSA, inspect the
proposal budget and accepted configuration; failure to reach a target fraction
does not imply equilibrium. Use deterministic seeds when comparing stochastic
runs. See :doc:`assumptions` and :doc:`theory`.

Optional scattering
-------------------

Install ``packlab[scattering]`` only when using the optical workflow. Keep its
PyMieSim dependency compatible with the versions declared by PackLab; the core
packing and analytical workflows do not require it.
