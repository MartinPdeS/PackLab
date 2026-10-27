:orphan:

API reference
=============

PackLab keeps compiled extensions alongside its public Python APIs. The
reference below is generated from the installed package, so signatures and
docstrings remain aligned with the released interfaces.

Start with one of these public namespaces:

* :mod:`PackLab.samplers` for radius distributions;
* :mod:`PackLab.monte_carlo` for random sequential adsorption (RSA);
* :mod:`PackLab.analytical` for Percus--Yevick mixture calculations;
* :mod:`PackLab.scattering` for the optional PyMieSim integration.

Only public classes and functions are listed below. Use the module-source link
next to each object when implementation details are useful.

Radius samplers
---------------

.. automodule:: PackLab.samplers
   :members:
   :member-order: bysource

Monte-Carlo hard-sphere workflows
---------------------------------

.. automodule:: PackLab.monte_carlo

Packing domain
~~~~~~~~~~~~~~

.. automodule:: PackLab.monte_carlo.domain
   :members:
   :member-order: bysource

RSA simulation
~~~~~~~~~~~~~~

.. automodule:: PackLab.monte_carlo.simulator
   :members:
   :member-order: bysource

Metropolis sampling
~~~~~~~~~~~~~~~~~~~

.. automodule:: PackLab.monte_carlo.metropolis
   :members:
   :member-order: bysource

Results and visualisation
~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: PackLab.monte_carlo.results
   :members:
   :member-order: bysource

Statistics and ensemble estimates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: PackLab.monte_carlo.statistics
   :members:
   :member-order: bysource

.. automodule:: PackLab.monte_carlo.estimator
   :members:
   :member-order: bysource

Analytical model
----------------

.. automodule:: PackLab.analytical

Mixture domain
~~~~~~~~~~~~~~

.. automodule:: PackLab.analytical.domain
   :members:
   :member-order: bysource

Wavenumber grids
~~~~~~~~~~~~~~~~

.. automodule:: PackLab.analytical.grid
   :members:
   :member-order: bysource

Percus--Yevick solver
~~~~~~~~~~~~~~~~~~~~~

.. automodule:: PackLab.analytical.solver
   :members:
   :member-order: bysource

Scattering tools
----------------

.. automodule:: PackLab.scattering

Data containers
~~~~~~~~~~~~~~~

.. automodule:: PackLab.scattering.data
   :members:
   :member-order: bysource

Amplitude generation
~~~~~~~~~~~~~~~~~~~~

.. automodule:: PackLab.scattering.model
   :members:
   :member-order: bysource

Phase-function plots
~~~~~~~~~~~~~~~~~~~~

.. automodule:: PackLab.scattering.plottings
   :members:
   :member-order: bysource
