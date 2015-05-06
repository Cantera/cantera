.. py:currentmodule:: cantera

Creating Phase Objects
======================

These classes are composite representations of a substance which has
thermodynamic, chemical kinetic, and (optionally) transport properties.

.. autoclass:: Solution(infile='', phaseid='', source=None, thermo=None, species=(), kinetics=None, reactions=(), **kwargs)

.. autoclass:: Interface(infile='', phaseid='', phases=(), thermo=None, species=(), kinetics=None, reactions=())

.. autoclass:: DustyGas(infile, phaseid='')

Utility Functions
-----------------

.. autofunction:: add_directory
