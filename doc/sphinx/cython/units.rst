.. py:currentmodule:: cantera.with_units

Python Interface With Units
===========================

This interface allows users to specify physical units associated with quantities.
To do so, this interface leverages the `pint <https://pint.readthedocs.io/en/stable/>`__
library to provide consistent unit conversion. Several examples of this interface can
be found in the ``samples/python`` folder in the 
`root of the repository <https://github.com/Cantera/cantera/tree/main/samples/python>`__.
Examples that use this interface are suffixed with `_units`.

The overall goal is to provide a compatible implementation of the `cantera.Solution` and
`cantera.PureFluid` interfaces. Please see those pages for further documentation of the
available properties.

Solution with Units
-------------------

.. autoclass:: Solution
   :no-members:

PureFluid Phases With Units
---------------------------

The following convenience classes are available to create `PureFluid <PureFluid>`
objects with the indicated equation of state:

.. autoclass:: CarbonDioxide
   :no-members:
.. autoclass:: Heptane
   :no-members:
.. autoclass:: Hfc134a
   :no-members:
.. autoclass:: Hydrogen
   :no-members:
.. autoclass:: Methane
   :no-members:
.. autoclass:: Nitrogen
   :no-members:
.. autoclass:: Oxygen
   :no-members:
.. autoclass:: Water
   :no-members:

The full documentation for the `PureFluid <PureFluid>` class and its properties is here:

.. autoclass:: PureFluid
   :no-members:
