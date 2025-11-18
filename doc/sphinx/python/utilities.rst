.. py:currentmodule:: cantera

Utilities
=========

.. contents::
   :local:

YAML Input/Output
-----------------

AnyMap
^^^^^^

.. autoclass:: AnyMap()

YamlWriter
^^^^^^^^^^
.. autoclass:: YamlWriter()

Unit Conversions
----------------

UnitSystem
^^^^^^^^^^
.. autoclass:: UnitSystem(units)
   :no-undoc-members:

Units
^^^^^
.. autoclass:: Units(name)
   :no-undoc-members:


.. _sec-python-jacobians:

Jacobians & Preconditioners
---------------------------

SystemJacobian
^^^^^^^^^^^^^^
.. autoclass:: SystemJacobian()

AdaptivePreconditioner
^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: AdaptivePreconditioner()

BandedJacobian
^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: BandedJacobian()

EigenSparseDirectJacobian
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: EigenSparseDirectJacobian()


.. _sec-python-global-funcs:

Global Functions
----------------

.. autofunction:: add_data_directory
.. autofunction:: get_data_directories
.. autofunction:: appdelete
.. autofunction:: use_sparse

.. autofunction:: make_deprecation_warnings_fatal
.. autofunction:: suppress_deprecation_warnings
.. autofunction:: suppress_thermo_warnings
.. autofunction:: use_legacy_rate_constants
.. autofunction:: debug_mode_enabled
.. autofunction:: print_stack_trace_on_segfault

.. autofunction:: extension(name: str)

Exceptions
----------

CanteraError
^^^^^^^^^^^^
.. autoclass:: CanteraError
   :no-undoc-members:


ThermoModelMethodError
^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: ThermoModelMethodError
   :no-undoc-members:
