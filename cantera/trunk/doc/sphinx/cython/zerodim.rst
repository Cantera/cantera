.. py:currentmodule:: cantera

Zero-Dimensional Reactor Networks
=================================

Defining Functions
------------------

.. autoclass:: Func1

Base Classes
------------

.. autoclass:: ReactorBase(contents=None, name=None)
.. autoclass:: FlowDevice(upstream, downstream, *, name=None)

Reactor Networks
----------------

.. autoclass:: ReactorNet(reactors=())

Reactors
--------

.. autoclass:: Reservoir(contents=None, name=None)
.. autoclass:: Reactor(contents=None, *, name=None, energy='on')
.. autoclass:: IdealGasReactor(contents=None, *, name=None, energy='on')
.. autoclass:: ConstPressureReactor(contents=None, *, name=None, energy='on')
.. autoclass:: IdealGasConstPressureReactor(contents=None, *, name=None, energy='on')
.. autoclass:: FlowReactor(contents=None, *, name=None, energy='on')

Flow Controllers
----------------

.. autoclass:: Wall(left, right, *, name=None, A=None, K=None, U=None, Q=None, velocity=None, kinetics=(None,None))
.. autoclass:: MassFlowController(upstream, downstream, *, name=None, mdot=None)
.. autoclass:: Valve(upstream, downstream, *, name=None, K=None)
.. autoclass:: PressureController(upstream, downstream, *, name=None, master=None, K=None)