.. default-domain:: mat
.. currentmodule:: OneDim

One-dimensional Reacting Flows
==============================
.. contents::
   :local:

Composite Domains
-----------------

CounterflowDiffusionFlame
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

Inlet
^^^^^
.. autoclass:: Inlet(phase[, name])

Outlet
^^^^^^
.. autoclass:: Outlet(phase[, name])

OutletRes
^^^^^^^^^
.. autoclass:: OutletRes(phase[, name])

SymmPlane
^^^^^^^^^
.. autoclass:: SymmPlane(phase[, name])

Surface
^^^^^^^
.. autoclass:: Surface(phase[, name])

ReactingSurface
^^^^^^^^^^^^^^^
.. autoclass:: ReactingSurface(phase[, name])


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
