.. py:currentmodule:: cantera

Thermodynamic Properties
========================

Phases
------

These classes are used to describe the thermodynamic state of a system.

.. autoclass:: ThermoPhase(infile='', phaseid='')
.. autoclass:: InterfacePhase(infile='', phaseid='')
.. autoclass:: PureFluid(infile='', phaseid='')

Mixture
-------

.. autoclass:: Mixture

Species
-------

.. autoclass:: Species

Species Thermodynamic Properties
--------------------------------

These classes are used to describe the reference-state thermodynamic properties
of a pure species.

.. autoclass:: SpeciesThermo(T_low, T_high, P_ref, coeffs)
.. autoclass:: ConstantCp(T_low, T_high, P_ref, coeffs)
    :no-undoc-members:

.. autoclass:: NasaPoly2(T_low, T_high, P_ref, coeffs)
    :no-undoc-members:

.. autoclass:: ShomatePoly2(T_low, T_high, P_ref, coeffs)
    :no-undoc-members:
