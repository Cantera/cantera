.. default-domain:: mat
.. currentmodule:: Reactor.+ct

Zero-Dimensional Reactor Networks
=================================
.. contents::
   :local:

Defining Functions
------------------
.. autoclass:: Func1(typ, varargin)

Reactor Networks
----------------
.. autoclass:: ReactorNet(reactors)

.. _sec-matlab-reactors:

Reactors
--------

ReactorBase
^^^^^^^^^^^
.. autoclass:: ReactorBase(id)

Reservoir
^^^^^^^^^
.. autoclass:: Reservoir(phase[, name][, clone=true])

Reactor
^^^^^^^
.. autoclass:: Reactor(phase[, name][, clone=true])

IdealGasReactor
^^^^^^^^^^^^^^^
.. autoclass:: IdealGasReactor(phase[, name][, clone=true])

ConstPressureReactor
^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ConstPressureReactor(phase[, name][, clone=true])

IdealGasConstPressureReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: IdealGasConstPressureReactor(phase[, name][, clone=true])

FlowReactor
^^^^^^^^^^^
.. autoclass:: FlowReactor(phase[, name][, clone=true])

Surfaces
--------

ReactorSurface
^^^^^^^^^^^^^^
.. autoclass:: ReactorSurface(surf, reactors[, name][, clone=true])

Connectors
----------

Connector
^^^^^^^^^
.. autoclass:: Connector(typ, r1, r2[, name])

Wall
^^^^
.. autoclass:: Wall(l, r[, name])

FlowDevice
^^^^^^^^^^
.. autoclass:: FlowDevice(typ, upstream, downstream[, name])

MassFlowController
^^^^^^^^^^^^^^^^^^
.. autoclass:: MassFlowController(upstream, downstream[, name])

Valve
^^^^^
.. autoclass:: Valve(upstream, downstream[, name])
