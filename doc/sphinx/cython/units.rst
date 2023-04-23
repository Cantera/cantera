.. py:currentmodule:: cantera.with_units

Python Interface With Units
===========================

This interface allows users to specify physical units associated with quantities.
To do so, this interface leverages the `pint <https://pint.readthedocs.io/en/stable/>`__
library to provide consistent unit conversion.

Solution with Units
-------------------

.. autoclass:: Solution

PureFluid Phases With Units
---------------------------

The following convenience classes are available to create `PureFluid <PureFluid>`
objects with the indicated equation of state:

.. autoclass:: CarbonDioxide
.. autoclass:: Heptane
.. autoclass:: Hfc134a
.. autoclass:: Hydrogen
.. autoclass:: Methane
.. autoclass:: Nitrogen
.. autoclass:: Oxygen
.. autoclass:: Water

The full documentation for the `PureFluid <PureFluid>` class and its properties is here:

.. autoclass:: PureFluid
