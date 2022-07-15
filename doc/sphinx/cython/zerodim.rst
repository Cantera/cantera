.. py:currentmodule:: cantera

.. _sec-cython-zerodim:

Zero-Dimensional Reactor Networks
=================================

.. contents::
   :local:

Defining Functions
------------------

.. autoclass:: Func1

.. autoclass:: TabulatedFunction

Base Classes
------------

ReactorBase
^^^^^^^^^^^
.. autoclass:: ReactorBase

FlowDevice
^^^^^^^^^^
.. autoclass:: FlowDevice

Reactor Networks
----------------

.. autoclass:: ReactorNet

Reactors
--------

Reservoir
^^^^^^^^^
.. autoclass:: Reservoir(contents=None, name=None)

Reactor
^^^^^^^
.. autoclass:: Reactor

IdealGasReactor
^^^^^^^^^^^^^^^
.. autoclass:: IdealGasReactor(contents=None, *, name=None, energy='on')

IdealGasMoleReactor
^^^^^^^^^^^^^^^^^^^
.. autoclass:: IdealGasMoleReactor(contents=None, *, name=None, energy='on')

ConstPressureReactor
^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ConstPressureReactor(contents=None, *, name=None, energy='on')

IdealGasConstPressureReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: IdealGasConstPressureReactor(contents=None, *, name=None, energy='on')

IdealGasConstPressureMoleReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: IdealGasConstPressureMoleReactor(contents=None, *, name=None, energy='on')

FlowReactor
^^^^^^^^^^^
.. autoclass:: FlowReactor(contents=None, *, name=None, energy='on')

ExtensibleReactor
^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleReactor(contents=None, *, name=None, energy='on')

ExtensibleIdealGasReactor
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleIdealGasReactor(contents=None, *, name=None, energy='on')

ExtensibleConstPressureReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleConstPressureReactor(contents=None, *, name=None, energy='on')

ExtensibleIdealGasConstPressureReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleIdealGasConstPressureReactor(contents=None, *, name=None, energy='on')

Walls
-----

Wall
^^^^
.. autoclass:: Wall(left, right, *, name=None, A=None, K=None, U=None, Q=None, velocity=None)
   :inherited-members:

Surfaces
--------

ReactorSurface
^^^^^^^^^^^^^^
.. autoclass:: ReactorSurface(kin=None, r=None, *, A=None)

Flow Controllers
----------------

MassFlowController
^^^^^^^^^^^^^^^^^^
.. autoclass:: MassFlowController
   :inherited-members:

Valve
^^^^^
.. autoclass:: Valve
   :inherited-members:

PressureController
^^^^^^^^^^^^^^^^^^
.. autoclass:: PressureController
   :inherited-members:

Preconditioners
---------------

AdaptivePreconditioner
^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: AdaptivePreconditioner
