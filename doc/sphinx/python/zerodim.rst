.. py:currentmodule:: cantera

.. _sec-python-zerodim:

Zero-Dimensional Reactor Networks
=================================

.. contents::
   :local:

Defining Functions
------------------

.. autoclass:: Func1

.. autoclass:: Tabulated1

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

.. _sec-python-reactors:

Reactors
--------

Reservoir
^^^^^^^^^
.. autoclass:: Reservoir(contents=None, name=None)

Reactor
^^^^^^^
.. autoclass:: Reactor(contents=None, *, name=None, energy='on', node_attr=None, group_name="")

MoleReactor
^^^^^^^^^^^
.. autoclass:: MoleReactor(contents=None, *, name=None, energy='on', node_attr=None, group_name="")

IdealGasReactor
^^^^^^^^^^^^^^^
.. autoclass:: IdealGasReactor(contents=None, *, name=None, energy='on', node_attr=None, group_name="")

IdealGasMoleReactor
^^^^^^^^^^^^^^^^^^^
.. autoclass:: IdealGasMoleReactor(contents=None, *, name=None, energy='on', node_attr=None, group_name="")

ConstPressureReactor
^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ConstPressureReactor(contents=None, *, name=None, energy='on', node_attr=None, group_name="")

ConstPressureMoleReactor
^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ConstPressureMoleReactor(contents=None, *, name=None, energy='on', node_attr=None, group_name="")

IdealGasConstPressureReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: IdealGasConstPressureReactor(contents=None, *, name=None, energy='on', node_attr=None, group_name="")

IdealGasConstPressureMoleReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: IdealGasConstPressureMoleReactor(contents=None, *, name=None, energy='on', node_attr=None, group_name="")

FlowReactor
^^^^^^^^^^^
.. autoclass:: FlowReactor(contents=None, *, name=None, energy='on', node_attr=None, group_name="")

ExtensibleReactor
^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleReactor(contents=None, *, name=None, energy='on', node_attr=None, group_name="")

ExtensibleIdealGasReactor
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleIdealGasReactor(contents=None, *, name=None, energy='on', node_attr=None, group_name="")

ExtensibleConstPressureReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleConstPressureReactor(contents=None, *, name=None, energy='on', node_attr=None, group_name="")

ExtensibleIdealGasConstPressureReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleIdealGasConstPressureReactor(contents=None, *, name=None, energy='on', node_attr=None, group_name="")

ExtensibleMoleReactor
^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleMoleReactor(contents=None, *, name=None, energy='on', node_attr=None, group_name="")

ExtensibleIdealGasMoleReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleIdealGasMoleReactor(contents=None, *, name=None, energy='on', node_attr=None, group_name="")

ExtensibleConstPressureMoleReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleConstPressureMoleReactor(contents=None, *, name=None, energy='on', node_attr=None, group_name="")

ExtensibleIdealGasConstPressureMoleReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleIdealGasConstPressureMoleReactor(contents=None, *, name=None, energy='on', node_attr=None, group_name="")

Walls
-----

Wall
^^^^
.. autoclass:: Wall(left, right, *, name=None, A=None, K=None, U=None, Q=None, velocity=None, edge_attr=None)
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

Drawing Reactor Networks
------------------------
These functions provide the implementation behind the ``draw`` methods of the
corresponding classes.

.. autofunction:: cantera.drawnetwork.draw_reactor

.. autofunction:: cantera.drawnetwork.draw_reactor_net

.. autofunction:: cantera.drawnetwork.draw_surface

.. autofunction:: cantera.drawnetwork.draw_flow_controllers

.. autofunction:: cantera.drawnetwork.draw_walls
