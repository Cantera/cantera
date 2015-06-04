.. py:currentmodule:: cantera

Creating Phase Objects
======================

These classes are composite representations of a substance which has
thermodynamic, chemical kinetic, and (optionally) transport properties.

.. autoclass:: Solution(infile='', phaseid='', source=None, thermo=None, species=(), kinetics=None, reactions=(), **kwargs)

.. autoclass:: Interface(infile='', phaseid='', phases=(), thermo=None, species=(), kinetics=None, reactions=())

.. autoclass:: DustyGas(infile, phaseid='')

Pure Fluid Phases
-----------------

The following convenience functions can be used to create `PureFluid` objects
with the indicated equation of state:

.. autofunction:: CarbonDioxide
.. autofunction:: Heptane
.. autofunction:: Hfc134a
.. autofunction:: Hydrogen
.. autofunction:: Methane
.. autofunction:: Nitrogen
.. autofunction:: Oxygen
.. autofunction:: Water


Utility Functions
-----------------

.. autofunction:: add_directory
