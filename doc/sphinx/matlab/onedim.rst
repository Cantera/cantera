.. default-domain:: mat
.. currentmodule:: +ct.+oneD

One-dimensional Reacting Flows
==============================

.. caution::
    The MATLAB toolbox is an experimental part of the Cantera API and may be changed
    without notice. It includes breaking changes from the legacy MATLAB API. While
    almost all features of the legacy MATLAB API are implemented, the toolbox does not
    include all functionality available for the C++ and Python interfaces.

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

Inlet
^^^^^
.. autoclass:: Inlet(phase[, name])

Outlet
^^^^^^
.. autoclass:: Outlet(phase[, name])

OutletRes
^^^^^^^^^
.. autoclass:: OutletRes(phase[, name])

SymmetryPlane
^^^^^^^^^^^^^
.. autoclass:: SymmetryPlane(phase[, name])

Surface
^^^^^^^
.. autoclass:: Surface(phase[, name])

ReactingSurface
^^^^^^^^^^^^^^^
.. autoclass:: ReactingSurface(phase[, name])


Base Classes
------------

Domain
^^^^^^
.. autoclass:: Domain

Boundary
^^^^^^^^
.. autoclass:: Boundary

Sim1D
^^^^^
.. autoclass:: Sim1D(domains)
