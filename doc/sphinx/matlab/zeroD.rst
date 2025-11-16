.. default-domain:: mat
.. currentmodule:: +ct.+zeroD

Zero-Dimensional Reactor Networks
=================================

.. caution::
    The MATLAB toolbox is an experimental part of the Cantera API and may be changed
    without notice. It includes breaking changes from the legacy MATLAB API. While
    almost all features of the legacy MATLAB API are implemented, the toolbox does not
    include all functionality available for the C++ and Python interfaces.

.. contents::
   :local:

Defining Functions
------------------
.. autoclass:: Base.+ct.Func1(typ, varargin)

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
