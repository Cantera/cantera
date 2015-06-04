.. py:currentmodule:: cantera

.. _sec-cython-onedim:

One-dimensional Reacting Flows
==============================

.. contents::
   :local:

Composite Domains
-----------------

FreeFlame
^^^^^^^^^
.. autoclass:: FreeFlame(gas, grid=None)

BurnerFlame
^^^^^^^^^^^
.. autoclass:: BurnerFlame(gas, grid=None)

CounterflowDiffusionFlame
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: CounterflowDiffusionFlame(gas, grid=None)

CounterflowPremixedFlame
^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: CounterflowPremixedFlame(gas, grid=None)

ImpingingJet
^^^^^^^^^^^^
.. autoclass:: ImpingingJet(gas, grid=None)

Flow Domains
------------

FreeFlow
^^^^^^^^
.. autoclass:: FreeFlow(thermo)
    :inherited-members:

AxisymmetricStagnationFlow
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: AxisymmetricStagnationFlow(thermo)
    :inherited-members:

Boundaries
----------

Inlet1D
^^^^^^^
.. autoclass:: Inlet1D(phase, name=None)

Outlet1D
^^^^^^^^
.. autoclass:: Outlet1D(phase, name=None)

OutletReservoir1D
^^^^^^^^^^^^^^^^^
.. autoclass:: OutletReservoir1D(phase, name=None)

SymmetryPlane1D
^^^^^^^^^^^^^^^
.. autoclass:: SymmetryPlane1D(phase, name=None)

Surface1D
^^^^^^^^^
.. autoclass:: Surface1D(phase, name=None)

ReactingSurface1D
^^^^^^^^^^^^^^^^^
.. autoclass:: ReactingSurface1D(phase, name=None)


Base Classes
------------
Domain1D
^^^^^^^^
.. autoclass:: Domain1D(name=None)

Boundary1D
^^^^^^^^^^
.. autoclass:: Boundary1D(phase, name=None)

Sim1D
^^^^^
.. autoclass:: Sim1D(domains)

FlameBase
^^^^^^^^^
.. autoclass:: FlameBase(domains, gas, grid=None)
