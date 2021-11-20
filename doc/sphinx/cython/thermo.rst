.. py:currentmodule:: cantera


Thermodynamic Properties
========================

.. contents::
   :local:

Phases
------

These classes are used to describe the thermodynamic state of a system.

ThermoPhase
^^^^^^^^^^^
.. autoclass:: ThermoPhase(infile='', phaseid='')

InterfacePhase
^^^^^^^^^^^^^^
.. autoclass:: InterfacePhase(infile='', phaseid='')

PureFluid
^^^^^^^^^
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

SpeciesThermo
^^^^^^^^^^^^^
.. autoclass:: SpeciesThermo(T_low, T_high, P_ref, coeffs)

ConstantCp
^^^^^^^^^^
.. autoclass:: ConstantCp(T_low, T_high, P_ref, coeffs)
    :no-undoc-members:

Mu0Poly
^^^^^^^
.. autoclass:: Mu0Poly(T_low, T_high, P_ref, coeffs)
    :no-undoc-members:

NasaPoly2
^^^^^^^^^
.. autoclass:: NasaPoly2(T_low, T_high, P_ref, coeffs)
    :no-undoc-members:

Nasa9PolyMultiTempRegion
^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: Nasa9PolyMultiTempRegion(T_low, T_high, P_ref, coeffs)
    :no-undoc-members:

ShomatePoly2
^^^^^^^^^^^^
.. autoclass:: ShomatePoly2(T_low, T_high, P_ref, coeffs)
    :no-undoc-members:

Element
-------

.. autoclass:: Element
    :no-undoc-members:

    .. autoattribute:: num_elements_defined
    .. autoattribute:: element_symbols
    .. autoattribute:: element_names
