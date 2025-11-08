.. default-domain:: mat
.. currentmodule:: OneDim.+ct

One-dimensional Reacting Flows
==============================
.. contents::
   :local:

Composite Domains
-----------------

CounterFlowDiffusionFlame
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: CounterFlowDiffusionFlame(left, flow, right, oxidizer)

.. _sec-matlab-flow-domains:

Flow Domains
------------

FreeFlow
^^^^^^^^
.. autoclass:: FreeFlow(phase[, name])

UnstrainedFlow
^^^^^^^^^^^^^^
.. autoclass:: UnstrainedFlow(phase[, name])

AxisymmetricFlow
^^^^^^^^^^^^^^^^
.. autoclass:: AxisymmetricFlow(phase[, name])

.. _sec-matlab-boundary-domains:

Boundaries
----------

Inlet1D
^^^^^^^
.. autoclass:: Inlet1D(phase[, name])

Outlet1D
^^^^^^^^
.. autoclass:: Outlet1D(phase[, name])

OutletRes1D
^^^^^^^^^^^
.. autoclass:: OutletRes1D(phase[, name])

SymmetryPlane1D
^^^^^^^^^^^^^^^
.. autoclass:: SymmetryPlane1D(phase[, name])

Surface1D
^^^^^^^^^
.. autoclass:: Surface1D(phase[, name])

ReactingSurface1D
^^^^^^^^^^^^^^^^^
.. autoclass:: ReactingSurface1D(phase[, name])


Base Classes
------------

Domain1D
^^^^^^^^
.. autoclass:: Domain1D

Boundary1D
^^^^^^^^^^
.. autoclass:: Boundary1D

Sim1D
^^^^^
.. autoclass:: Sim1D(domains)
