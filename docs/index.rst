.. Straph documentation master file, created by
sphinx-quickstart on Tue Mar 16 18:19:07 2021.
You can adapt this file completely to your liking, but it should at least
contain the root `toctree` directive.

Straph's documentation
======================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

``Straph`` is an open source Python 3 package, under the licence Apache 2.0, for exploration
and analysis of real and artificial stream graphs. This library provides specific
data structures for representing different types of stream graphs, algorithms to compute
basic properties and measures, readers and writers for various data formats as
well as random generators.

``Straph`` can be used to teach stream graph theory or illustrate particular concepts.
Several Jupyter notebooks have been written in order to demonstrate its functionalities.
Straph can be used by users or developers that are not necessarily experts in
programming or in stream graph theory.

Stream Graph
============

Stream graph formalism aims at extending basic graph theory concepts where time and structure are equally important.
Stream graphs are particularly suited to model temporal networks with a highly internal dynamics.

We refer to latapy_2017_ and latapy_2019_ for a detailed presentation of stream graph theory.


Summary
=======

.. toctree::
   :maxdepth: 1

   installation.rst
   notebooks/Getting Started.ipynb
   tutorials.rst
   api_reference.rst


.. _latapy_2017: https://arxiv.org/abs/1710.04073
.. _latapy_2019: https://arxiv.org/abs/1906.04840
