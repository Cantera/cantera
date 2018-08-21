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
.. autoclass:: FreeFlame(gas, grid=None, width=None)

BurnerFlame
^^^^^^^^^^^
.. autoclass:: BurnerFlame(gas, grid=None, width=None)

CounterflowDiffusionFlame
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: CounterflowDiffusionFlame(gas, grid=None, width=None)

CounterflowPremixedFlame
^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: CounterflowPremixedFlame(gas, grid=None, width=None)

ImpingingJet
^^^^^^^^^^^^
.. autoclass:: ImpingingJet(gas, grid=None, width=None)

IonFreeFlame
^^^^^^^^^^^^
.. autoclass:: IonFreeFlame(gas, grid=None, width=None)

   .. autoattribute:: E
   .. autoattribute:: electric_field_enabled
   .. automethod:: solve

IonBurnerFlame
^^^^^^^^^^^^^^
.. autoclass:: IonBurnerFlame(gas, grid=None, width=None)

   .. autoattribute:: E
   .. autoattribute:: electric_field_enabled
   .. automethod:: solve

Flow Domains
------------

IdealGasFlow
^^^^^^^^^^^^
.. autoclass:: IdealGasFlow(thermo)
    :inherited-members:

IonFlow
^^^^^^^
.. autoclass:: IonFlow(thermo)

FreeFlow
^^^^^^^^
.. autoclass:: FreeFlow(thermo)

AxisymmetricStagnationFlow
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: AxisymmetricStagnationFlow(thermo)

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
