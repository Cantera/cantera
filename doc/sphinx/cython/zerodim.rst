.. py:currentmodule:: cantera

.. _sec-cython-zerodim:

Zero-Dimensional Reactor Networks
=================================

.. contents::
   :local:

Defining Functions
------------------

.. autoclass:: Func1

Base Classes
------------

ReactorBase
^^^^^^^^^^^
.. autoclass:: ReactorBase(contents=None, name=None)

FlowDevice
^^^^^^^^^^
.. autoclass:: FlowDevice(upstream, downstream, *, name=None)

Reactor Networks
----------------

.. autoclass:: ReactorNet(reactors=())

Reactors
--------

Reservoir
^^^^^^^^^
.. autoclass:: Reservoir(contents=None, name=None)

Reactor
^^^^^^^
.. autoclass:: Reactor(contents=None, *, name=None, energy='on')

IdealGasReactor
^^^^^^^^^^^^^^^
.. autoclass:: IdealGasReactor(contents=None, *, name=None, energy='on')

ConstPressureReactor
^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ConstPressureReactor(contents=None, *, name=None, energy='on')

IdealGasConstPressureReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: IdealGasConstPressureReactor(contents=None, *, name=None, energy='on')

FlowReactor
^^^^^^^^^^^
.. autoclass:: FlowReactor(contents=None, *, name=None, energy='on')

Walls
-----

Wall
^^^^
.. autoclass:: Wall(left, right, *, name=None, A=None, K=None, U=None, Q=None, velocity=None, kinetics=(None,None))

WallSurface
^^^^^^^^^^^
.. autoclass:: WallSurface(wall, side)

Surfaces
--------

ReactorSurface
^^^^^^^^^^^^^^
.. autoclass:: ReactorSurface(kin=None, r=None, *, A=None)

Flow Controllers
----------------

MassFlowController
^^^^^^^^^^^^^^^^^^
.. autoclass:: MassFlowController(upstream, downstream, *, name=None, mdot=None)

Valve
^^^^^
.. autoclass:: Valve(upstream, downstream, *, name=None, K=None)

PressureController
^^^^^^^^^^^^^^^^^^
.. autoclass:: PressureController(upstream, downstream, *, name=None, master=None, K=None)
