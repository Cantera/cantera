.. py:currentmodule:: cantera.ctml_writer

.. _sec-species:

********************
Elements and Species
********************

.. _sec-elements:

Elements
========

The :class:`element` entry defines an element or an isotope of an element. Note that
these entries are not often needed, since the the database file ``elements.xml``
is searched for element definitions when importing phase and interface
definitions.  An explicit element entry is needed only if an isotope not in
``elements.xml`` is required::

    element(symbol='C-13',
            atomic_mass=13.003354826)
    element("O-18", 17.9991603)

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
            # ...,
            )

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
transport database parameters described by Kee et al. [#Kee1986]_. Entries in
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

Currently, several types are implemented which provide species properties
appropriate for models of ideal gas mixtures, ideal solutions, and pure
compounds.

The NASA 7-Coefficient Polynomial Parameterization
--------------------------------------------------

The NASA 7-coefficient polynomial parameterization is used to compute the
species reference-state thermodynamic properties :math:`\hat{c}^0_p(T)`,
:math:`\hat{h}^0(T)` and :math:`\hat{s}^0(T)`.

The NASA parameterization represents :math:`\hat{c}^0_p(T)` with a fourth-order
polynomial:

.. math::

    \frac{c_p^0(T)}{R} = a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4

    \frac{h^0(T)}{RT} = a_0 + \frac{a1}{2}T + \frac{a_2}{3} T^2 +
                        \frac{a_3}{4} T^3 + \frac{a_4}{5} T^4 + \frac{a_5}{T}

    \frac{s^0(T)}{R} = a_0 \ln T + a_1 T + \frac{a_2}{2} T^2 + \frac{a_3}{3} T^3 +
                       \frac{a_4}{4} T^4 + a_6

Note that this is the "old" NASA polynomial form, used in the original NASA
equilibrium program and in Chemkin, which uses 7 coefficients in each of two
temperature regions. It is not compatible with the form used in the most recent
version of the NASA equilibrium program, which uses 9 coefficients for each
temperature region.

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

The NASA 9-Coefficient Polynomial Parameterization
--------------------------------------------------

The NASA 9-coefficient polynomial parameterization [#McBride2002]_ ("NASA9" for
short) is an extension of the NASA 7-coefficient polynomial parameterization
which includes two additional terms in each temperature region, as well as
supporting an arbitrary number of temperature regions.

The NASA9 parameterization represents the species thermodynamic properties with
the following equations:

.. math::

    \frac{C_p^0(T)}{R} = a_0 T^{-2} + a_1 T^{-1} + a_2 + a_3 T
                     + a_4 T^2 + a_5 T^3 + a_6 T^4

    \frac{H^0(T)}{RT} = - a_0 T^{-2} + a_1 \frac{\ln T}{T} + a_2
        + \frac{a_3}{2} T + \frac{a_4}{3} T^2  + \frac{a_5}{4} T^3 +
        \frac{a_6}{5} T^4 + \frac{a_7}{T}

    \frac{s^0(T)}{R} = - \frac{a_0}{2} T^{-2} - a_1 T^{-1} + a_2 \ln T
       + a_3 T + \frac{a_4}{2} T^2 + \frac{a_5}{3} T^3  + \frac{a_6}{4} T^4 + a_8

The following is an example of a species defined using the NASA9
parameterization in three different temperature regions::

    species(name=u'CO2',
            atoms='C:1 O:2',
            thermo=(NASA9([200.00, 1000.00],
                          [ 4.943650540E+04, -6.264116010E+02,  5.301725240E+00,
                            2.503813816E-03, -2.127308728E-07, -7.689988780E-10,
                            2.849677801E-13, -4.528198460E+04, -7.048279440E+00]),
                    NASA9([1000.00, 6000.00],
                          [ 1.176962419E+05, -1.788791477E+03,  8.291523190E+00,
                           -9.223156780E-05,  4.863676880E-09, -1.891053312E-12,
                            6.330036590E-16, -3.908350590E+04, -2.652669281E+01]),
                    NASA9([6000.00, 20000.00],
                          [-1.544423287E+09,  1.016847056E+06, -2.561405230E+02,
                            3.369401080E-02, -2.181184337E-06,  6.991420840E-11,
                           -8.842351500E-16, -8.043214510E+06,  2.254177493E+03])),
            note='Gurvich,1991 pt1 p27 pt2 p24. [g 9/99]')

Thermodynamic data for a range of species can be obtained from the `NASA
ThermoBuild <http://cearun.grc.nasa.gov/cea/index_ds.html>`_ tool. Using the web
interface, an input file can be obtained for a set of species. This input file
should then be modified so that the first line reads "`thermo nasa9`", as in the
following example::

    thermo nasa9
       200.000  1000.000  6000.000 20000.000   9/09/04
    CO                Gurvich,1979 pt1 p25 pt2 p29.
     3 tpis79 C   1.00O   1.00    0.00    0.00    0.00 0   28.0101000    -110535.196
        200.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         8671.104
     1.489045326D+04-2.922285939D+02 5.724527170D+00-8.176235030D-03 1.456903469D-05
    -1.087746302D-08 3.027941827D-12                -1.303131878D+04-7.859241350D+00
       1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         8671.104
     4.619197250D+05-1.944704863D+03 5.916714180D+00-5.664282830D-04 1.398814540D-07
    -1.787680361D-11 9.620935570D-16                -2.466261084D+03-1.387413108D+01
       6000.000  20000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         8671.104
     8.868662960D+08-7.500377840D+05 2.495474979D+02-3.956351100D-02 3.297772080D-06
    -1.318409933D-10 1.998937948D-15                 5.701421130D+06-2.060704786D+03
    CO2               Gurvich,1991 pt1 p27 pt2 p24.
     3 g 9/99 C   1.00O   2.00    0.00    0.00    0.00 0   44.0095000    -393510.000
        200.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         9365.469
     4.943650540D+04-6.264116010D+02 5.301725240D+00 2.503813816D-03-2.127308728D-07
    -7.689988780D-10 2.849677801D-13                -4.528198460D+04-7.048279440D+00
       1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         9365.469
     1.176962419D+05-1.788791477D+03 8.291523190D+00-9.223156780D-05 4.863676880D-09
    -1.891053312D-12 6.330036590D-16                -3.908350590D+04-2.652669281D+01
       6000.000  20000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         9365.469
    -1.544423287D+09 1.016847056D+06-2.561405230D+02 3.369401080D-02-2.181184337D-06
     6.991420840D-11-8.842351500D-16                -8.043214510D+06 2.254177493D+03
    END PRODUCTS
    END REACTANTS

This file (saved for example as `nasathermo.dat`) can then be converted to the
CTI format using the `ck2cti` script::

    ck2cti --thermo=nasathermo.dat

To generate a full phase definition, create an input file defining the phase as
well, saved for example as `nasa.inp`::

    elements
    C O
    end

    species
    CO CO2
    end

The two input files can then be converted together by calling::

    ck2cti --input=nasa.inp --thermo=nasathermo.dat


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
\hat{h}^0(T_0), \hat{s}^0(T)`. The default value of :math:`T_0` is 298.15 K; the
default value for the other parameters is 0.0.

Example::

    thermo = const_cp(h0=(-393.51, 'kJ/mol'),
                      s0=(213.785, 'J/mol/K'),
                      cp0=(37.12, 'J/mol/K'))

Assuming that the :func:`units` function has been used to set the default energy
units to Joules and the default quantity unit to kmol, this may be equivalently
written as::

    thermo = const_cp(h0=-3.9351e8, s0=2.13785e5, cp0=3.712e4)

.. See ##REF## for more examples of use of this parameterization.

.. rubric:: References

.. [#Kee1986] R. J. Kee, G. Dixon-Lewis, J. Warnatz, M. E. Coltrin, and J. A. Miller.
   A FORTRAN Computer Code Package for the Evaluation of Gas-Phase, Multicomponent
   Transport Properties. Technical Report SAND86-8246, Sandia National Laboratories, 1986.

.. [#Mcbride2002] B. J. McBride, M. J. Zehe, S. Gordon. "NASA Glenn Coefficients
   for Calculating Thermodynamic Properties of Individual Species,"
   NASA/TP-2002-211556, Sept. 2002.
