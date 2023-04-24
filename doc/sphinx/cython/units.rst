.. py:currentmodule:: cantera.with_units

Python Interface With Units
===========================

This interface allows users to specify physical units associated with quantities.
To do so, this interface leverages the `pint <https://pint.readthedocs.io/en/stable/>`__
library to provide consistent unit conversion.

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
