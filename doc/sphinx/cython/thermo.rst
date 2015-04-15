.. py:currentmodule:: cantera

Thermodynamic Properties
========================

These classes are used to describe the thermodynamic state of a system.

.. autoclass:: ThermoPhase(infile='', phaseid='')
.. autoclass:: InterfacePhase(infile='', phaseid='')
.. autoclass:: PureFluid(infile='', phaseid='')
.. autoclass:: Mixture

Species Thermodynamic Properties
================================

These classes are used to describe the reference-state thermodynamic properties
of a pure species.

.. autoclass:: SpeciesThermo
.. autoclass:: ConstantCp(T_low, T_high, P_ref, coeffs)
.. autoclass:: NasaPoly2(T_low, T_high, P_ref, coeffs)
.. autoclass:: ShomatePoly2(T_low, T_high, P_ref, coeffs)
