.. py:currentmodule:: ctml_writer

.. _sec-species:

********************
Elements and Species
********************

.. _sec-elements:

Elements
========

The element entry defines an element or an isotope of an element. Note that
these entries are not often needed, since the the database file ``elements.xml``
is searched for element definitions when importing phase and interface
definitions.  An explicit element entry is needed only if an isotope not in
``elements.xml`` is required::

    element(symbol='C-13',
            atomic_mass=13.003354826)
    element("O-!8", 17.9991603)

Species
=======

For each species, a :class:`species` entry is required. Species are defined at
the top-level of the input file---their definitions are not embedded in a phase
or interface entry.

Species Name
------------

The name field may contain embedded parentheses, ``+`` or ``-`` signs to
indicate the charge, or just about anything else that is printable and not a
reserved character in XML.  Some example name specifications::

    name = 'CH4'
    name = 'methane'
    name = 'argon_2+'
    name = 'CH2(singlet)'

Elemental Composition
---------------------

The elemental composition is specified in the atoms entry, as follows::

    atoms = "C:1 O:2"             # CO2
    atoms = "C:1, O:2"            # CO2 with optional comma
    atoms = "Y:1 Ba:2 Cu:3 O:6.5" # stoichiometric YBCO
    atoms = ""                    # a surface species representing an empty site
    atoms = "Ar:1 E:-2"           # Ar++

For gaseous species, the elemental composition is well-defined, since the
species represent distinct molecules. For species in solid or liquid solutions,
or on surfaces, there may be several possible ways of defining the species. For
example, an aqueous species might be defined with or without including the water
molecules in the solvation cage surrounding it.

For surface species, it is possible to omit the ``atoms`` field entirely, in
which case it is composed of nothing, and represents an empty surface site. This
can also be done to represent vacancies in solids. A charged vacancy can be
defined to be composed solely of electrons::

    species(name = 'ysz-oxygen-vacancy',
            atoms = 'O:0, E:2',
            ...)

Note that an atom number of zero may be given if desired, but is completely
equivalent to omitting that element.

The number of atoms of an element must be non-negative, except for the special
"element" ``E`` that represents an electron.

Thermodynamic Properties
------------------------

The :class:`phase` and :class:`ideal_interface` entries discussed in the last
chapter implement specific models for the thermodynamic properties appropriate
for the type of phase or interface they represent. Although each one may use
different expressions to compute the properties, they all require thermodynamic
property information for the individual species. For the phase types implemented
at present, the properties needed are:

1. the molar heat capacity at constant pressure :math:`\hat{c}^0_p(T)` for a
   range of temperatures and a reference pressure :math:`P_0`;
2. the molar enthalpy :math:`\hat{h}(T_0, P_0)` at :math:`P_0` and a reference
   temperature :math:`T_0`;
3. the absolute molar entropy :math:`\hat{s}(T_0, P_0)` at :math:`(T_0, P_0)`.

See: :ref:`sec-thermo-models`

.. _sec-species-transport-models:

Species Transport Coefficients
------------------------------

Transport property models in general require coefficients that express the
effect of each species on the transport properties of the phase. The
``transport`` field may be assigned an embedded entry that provides
species-specific coefficients.

Currently, the only entry type is :class:`gas_transport`, which supplies
parameters needed by the ideal-gas transport property models. The field values
and their units of the :class:`gas_transport` entry are compatible with the
transport database parameters described by Kee et al. [1986]. Entries in
transport databases in the format described in their report can be used directly
in the fields of the :class:`gas_transport` entry, without requiring any unit
conversion. The numeric field values should all be entered as pure numbers, with
no attached units string.

.. _sec-thermo-models:

Thermodynamic Property Models
=============================

The entry types described in this section can be used to provide data for the
``thermo`` field of a :class:`species`. Each implements a different
*parameterization* (functional form) for the heat capacity. Note that there is
no requirement that all species in a phase use the same parameterization; each
species can use the one most appropriate to represent how the heat capacity
depends on temperature.

Currently, three entry types are implemented, all of which provide species
properties appropriate for models of ideal gas mixtures, ideal solutions, and
pure compounds. Non-ideal phase models are not yet implemented, but may be in
future releases. When they are, additional entry types may also be added that
provide species-specific coefficients required by specific non-ideal equations
of state.

The NASA Polynomial Parameterization
------------------------------------

The NASA polynomial parameterization is used to compute the species
reference-state thermodynamic properties :math:`\hat{c}^0_p(T)`,
:math:`\hat{h}^0(T)` and :math:`\hat{s}^0(T)`.

The NASA parameterization represents :math:`\hat{c}^0_p(T)` with a fourth-order
polynomial:

.. math::

    \frac{c_p^0(T)}{R} = a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4

    \frac{h^0(T)}{RT} = a_0 + \frac{a1}{2}T + \frac{a_2}{3} T^2 +
                        \frac{a_3}{4} T^3 + \frac{a_4}{5} T^4 + a_5

    \frac{s^0(T)}{R} = a_o \ln T + a_1 T + \frac{a_2}{2} T^2 + \frac{a_3}{3} T^3 +
                       \frac{a_4}{4} T^4 + a_6

Note that this is the "old" NASA polynomial form, used in the original NASA
equilibrium program and in Chemkin. It is not compatible with the form used in
the most recent version of the NASA equilibrium program, which uses 9
coefficients, not 7.

A NASA parameterization is defined by an embedded :class:`NASA` entry. Very
often, two NASA parameterizations are used for two contiguous temperature
ranges. This can be specified by assigning the ``thermo`` field of the
``species`` entry a sequence of two :class:`NASA` entries::

    # use one NASA parameterization for T < 1000 K, and another for T > 1000 K.
    species(name = "O2",
            atoms = " O:2 ",
            thermo = (
                NASA( [ 200.00, 1000.00], [ 3.782456360E+00, -2.996734160E-03,
		        9.847302010E-06, -9.681295090E-09, 3.243728370E-12,
                        -1.063943560E+03, 3.657675730E+00] ),
                NASA( [ 1000.00, 3500.00], [ 3.282537840E+00, 1.483087540E-03,
                        -7.579666690E-07, 2.094705550E-10, -2.167177940E-14,
                        -1.088457720E+03, 5.453231290E+00] ) ) )

The Shomate Parameterization
----------------------------

The Shomate parameterization is:

.. math::

    \hat{c}_p^0(T) = A + Bt + Ct^2 + Dt^3 + \frac{E}{t^2}

    \hat{h}^0(T) = At + \frac{Bt^2}{2} + \frac{Ct^3}{3} + \frac{Dt^4}{4} -
                   \frac{E}{t} + F

    \hat{s}^0(T) = A \ln t + B t + \frac{Ct^2}{2} + \frac{Dt^3}{3} -
                   \frac{E}{2t^2} + G

where :math:`t = T / 1000 K`. It requires 7 coefficients A, B, C, D, E, F, and
G. This parameterization is used to represent reference-state properties in the
`NIST Chemistry WebBook <http://webbook.nist.gov/chemistry>`_. The values of the
coefficients A through G should be entered precisely as shown there, with no
units attached. Unit conversions to SI will be handled internally.

Example usage of the :class:`Shomate` directive::

    # use a single Shomate parameterization.
    species(name = "O2",
            atoms = " O:2 ",
            thermo = Shomate( [298.0, 6000.0],
                              [29.659, 6.137261, -1.186521, 0.09578, -0.219663,
                               -9.861391, 237.948] ) )

Constant Heat Capacity
----------------------

In some cases, species properties may only be required at a single temperature
or over a narrow temperature range. In such cases, the heat capacity can be
approximated as constant, and simpler expressions can be used for the thermodynamic
properties. The :class:`const_cp` parameterization computes the properties as
follows:

.. math::

    \hat{c}_p^0(T) = \hat{c}_p^0(T_0)

    \hat{h}^0(T) = \hat{h}^0(T_0) + \hat{c}_p^0\cdot(T-T_0)

    \hat{s}^0(T) = \hat{s}^0(T_0) + \hat{c}_p^0 \ln (T/T_0)

The parameterization uses four constants: :math:`T_0, \hat{c}_p^0(T_0),
\hat{h}^0(T_0), \hat{s}^0(T)`.

Example::

    thermo = const_cp( t0 = 1200.0,
                       h0 = (-5.0, 'kcal/mol') )

.. See ##REF## for more examples of use of this parameterization.
