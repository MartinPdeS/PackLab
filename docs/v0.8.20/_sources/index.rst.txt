PackLab documentation
=====================

Explore structure in three-dimensional hard-sphere systems. Compute a Percus–Yevick equilibrium reference, generate an irreversible RSA packing, or equilibrate a configuration with Metropolis moves.

**Version:** |version|

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Get started
      :link: getting_started
      :link-type: doc

      Install the package and run a first model.

   .. grid-item-card:: Choose a workflow
      :link: guide
      :link-type: doc

      Find the approach that matches your scientific question.

   .. grid-item-card:: Explore examples
      :link: gallery/index
      :link-type: doc

      Run worked examples and inspect their outputs.

   .. grid-item-card:: Look up the API
      :link: code
      :link-type: doc

      Find public classes, parameters and result helpers.

Before you interpret results
----------------------------

RSA is an irreversible deposition process. Metropolis moves sample an
equilibrium system, while Percus–Yevick supplies an approximate analytical
reference. Match physical inputs before comparing them; see :doc:`assumptions`.

.. toctree::
   :maxdepth: 2
   :hidden:

   Start <getting_started>
   Guide <guide>
   Examples <gallery/index>
   API <code>
   Resources <resources>
