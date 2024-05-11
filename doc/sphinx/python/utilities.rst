.. py:currentmodule:: cantera

Utilities
=========

.. contents::
   :local:

YAML Input/Output
-----------------

AnyMap
^^^^^^

.. autoclass:: AnyMap

YamlWriter
^^^^^^^^^^
.. autoclass:: YamlWriter

Unit Conversions
----------------

UnitSystem
^^^^^^^^^^
.. autoclass:: UnitSystem
   :no-undoc-members:

Units
^^^^^
.. autoclass:: Units
   :no-undoc-members:

.. _sec-python-global-funcs:

Global Functions
----------------

.. autofunction:: add_directory
.. autofunction:: get_data_directories
.. autofunction:: import_phases
.. autofunction:: appdelete
.. autofunction:: use_sparse

.. autofunction:: make_deprecation_warnings_fatal
.. autofunction:: suppress_deprecation_warnings
.. autofunction:: suppress_thermo_warnings
.. autofunction:: use_legacy_rate_constants
.. autofunction:: debug_mode_enabled
.. autofunction:: add_module_directory

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
