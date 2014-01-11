.. py:currentmodule:: cantera

One-dimensional Reacting Flows
==============================

Composite Domains
-----------------
.. autoclass:: FreeFlame(gas, grid=None)
.. autoclass:: BurnerFlame(gas, grid=None)
.. autoclass:: CounterflowDiffusionFlame(gas, grid=None)

Flow Domains
------------
.. autoclass:: FreeFlow(thermo)
    :inherited-members:
.. autoclass:: AxisymmetricStagnationFlow(thermo)
    :inherited-members:

Boundaries
----------
.. autoclass:: Inlet1D(phase, name=None)
.. autoclass:: Outlet1D(phase, name=None)
.. autoclass:: OutletReservoir1D(phase, name=None)
.. autoclass:: SymmetryPlane1D(phase, name=None)
.. autoclass:: Surface1D(phase, name=None)
.. autoclass:: ReactingSurface1D(phase, name=None)

Base Classes
------------
.. autoclass:: Domain1D(name=None)
.. autoclass:: Boundary1D(phase, name=None)
.. autoclass:: Sim1D(domains)
.. autoclass:: FlameBase(domains, gas, grid=None)
