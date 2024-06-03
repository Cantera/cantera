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
.. autoclass:: ThermoPhase(infile='', name='', *, origin=None, yaml=None, thermo=None, species=(), **kwargs)

InterfacePhase
^^^^^^^^^^^^^^
.. autoclass:: InterfacePhase(infile='', name='', adjacent=(), *, origin=None, yaml=None, thermo=None, species=(), **kwargs)

PureFluid
^^^^^^^^^
.. autoclass:: PureFluid(infile='', name='')

Mixture
-------

.. autoclass:: Mixture(phases)

Species
-------

.. autoclass:: Species(name=None, composition=None, charge=None, size=None, init=True)

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

.. autoclass:: Element(name_or_symbol_or_atomic_number)
    :no-undoc-members:

    .. autoattribute:: num_elements_defined
    .. autoattribute:: element_symbols
    .. autoattribute:: element_names
