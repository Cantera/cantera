.. py:currentmodule:: cantera

.. _sec-python-zerodim:

Zero-Dimensional Reactor Networks
=================================

.. contents::
   :local:

Defining Functions
------------------

.. autoclass:: Func1(callable)

.. autoclass:: Tabulated1(time, fval, method='linear')

Base Classes
------------

ReactorBase
^^^^^^^^^^^
.. autoclass:: ReactorBase()

FlowDevice
^^^^^^^^^^
.. autoclass:: FlowDevice()

Reactor Networks
----------------

.. autoclass:: ReactorNet(reactors=())

.. _sec-python-reactors:

Reactors
--------

Reservoir
^^^^^^^^^
.. autoclass:: Reservoir(phase, name=None)

Reactor
^^^^^^^
.. autoclass:: Reactor(phase, *, clone=None, name=None, energy='on', node_attr=None, group_name="")

MoleReactor
^^^^^^^^^^^
.. autoclass:: MoleReactor(phase, *, clone=None, name=None, energy='on', node_attr=None, group_name="")

IdealGasReactor
^^^^^^^^^^^^^^^
.. autoclass:: IdealGasReactor(phase, *, clone=None, name=None, energy='on', node_attr=None, group_name="")

IdealGasMoleReactor
^^^^^^^^^^^^^^^^^^^
.. autoclass:: IdealGasMoleReactor(phase, *, clone=None, name=None, energy='on', node_attr=None, group_name="")

ConstPressureReactor
^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ConstPressureReactor(phase, *, clone=None, name=None, energy='on', node_attr=None, group_name="")

ConstPressureMoleReactor
^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ConstPressureMoleReactor(phase, *, clone=None, name=None, energy='on', node_attr=None, group_name="")

IdealGasConstPressureReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: IdealGasConstPressureReactor(phase, *, clone=None, name=None, energy='on', node_attr=None, group_name="")

IdealGasConstPressureMoleReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: IdealGasConstPressureMoleReactor(phase, *, clone=None, name=None, energy='on', node_attr=None, group_name="")

FlowReactor
^^^^^^^^^^^
.. autoclass:: FlowReactor(phase, *, clone=None, name=None, energy='on', node_attr=None, group_name="")

ExtensibleReactor
^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleReactor(phase, *, clone=None, name=None, energy='on', node_attr=None, group_name="")

ExtensibleIdealGasReactor
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleIdealGasReactor(phase, *, clone=None, name=None, energy='on', node_attr=None, group_name="")

ExtensibleConstPressureReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleConstPressureReactor(phase, *, clone=None, name=None, energy='on', node_attr=None, group_name="")

ExtensibleIdealGasConstPressureReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleIdealGasConstPressureReactor(phase, *, clone=None, name=None, energy='on', node_attr=None, group_name="")

ExtensibleMoleReactor
^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleMoleReactor(phase, *, clone=None, name=None, energy='on', node_attr=None, group_name="")

ExtensibleIdealGasMoleReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleIdealGasMoleReactor(phase, *, clone=None, name=None, energy='on', node_attr=None, group_name="")

ExtensibleConstPressureMoleReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleConstPressureMoleReactor(phase, *, clone=None, name=None, energy='on', node_attr=None, group_name="")

ExtensibleIdealGasConstPressureMoleReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleIdealGasConstPressureMoleReactor(phase, *, clone=None, name=None, energy='on', node_attr=None, group_name="")

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
.. autoclass:: ReactorSurface(phase, r=None, clone=None, name="(none)", *, A=None)

ExtensibleReactorSurface
^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleReactorSurface(phase, r=None, clone=None, name="(none)", *, A=None)

ExtensibleMoleReactorSurface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ExtensibleMoleReactorSurface(phase, r=None, clone=None, name="(none)", *, A=None)

Flow Controllers
----------------

MassFlowController
^^^^^^^^^^^^^^^^^^
.. autoclass:: MassFlowController(upstream, downstream, *, name=None, mdot=1.0, edge_attr=None)
   :inherited-members:

Valve
^^^^^
.. autoclass:: Valve(upstream, downstream, *, name=None, K=1.0, edge_attr=None)
   :inherited-members:

PressureController
^^^^^^^^^^^^^^^^^^
.. autoclass:: PressureController(upstream, downstream, *, name=None, primary=None, K=1.0, edge_attr=None)
   :inherited-members:

Drawing Reactor Networks
------------------------
These functions provide the implementation behind the ``draw`` methods of the
corresponding classes.

.. autofunction:: cantera.drawnetwork.draw_reactor

.. autofunction:: cantera.drawnetwork.draw_reactor_net

.. autofunction:: cantera.drawnetwork.draw_surface

.. autofunction:: cantera.drawnetwork.draw_flow_controllers

.. autofunction:: cantera.drawnetwork.draw_walls
