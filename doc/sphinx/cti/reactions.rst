.. py:currentmodule:: cantera.ctml_writer

.. _sec-reactions:

*********
Reactions
*********

Cantera supports a number of different types of reactions, including several
types of homogeneous reactions, surface reactions, and electrochemical
reactions. For each, there is a corresponding entry type. The simplest entry
type is :class:`reaction`, which can be used for any homogeneous reaction that
has a rate expression that obeys the law of mass action, with a rate coefficient
that depends only on temperature.

Common Attributes
=================

All of the entry types that define reactions share some common features. These
are described first, followed by descriptions of the individual reaction types
in the following sections.

The Reaction Equation
---------------------

The reaction equation determines the reactant and product stoichiometry. A
relatively simple parsing strategy is currently used, which assumes that all
coefficient and species symbols on either side of the equation are delimited by
spaces::

    2 CH2 <=> CH + CH3      # OK
    2 CH2<=>CH + CH3        # OK
    2CH2 <=> CH + CH3       # error
    CH2 + CH2 <=> CH + CH3  # OK
    2 CH2 <=> CH+CH3        # error

The incorrect versions here would generate "undeclared species" errors and would
halt processing of the input file. In the first case, the error would be that
the species ``2CH2`` is undeclared, and in the second case it would be species
``CH+CH3``.

Whether the reaction is reversible or not is determined by the form of the
equality sign in the reaction equation. If either ``<=>`` or ``=`` is found,
then the reaction is regarded as reversible, and the reverse rate will be
computed from detailed balance. If, on the other hand, ``=>`` is found, the
reaction will be treated as irreversible.

The rate coefficient is specified with an embedded entry corresponding to the
rate coefficient type. At present, the only implemented type is the modified
Arrhenius function

.. math::

    k_f(T) = A T^b \exp(-E/\hat{R}T)

which is defined with an :class:`Arrhenius` entry::

    rate_coeff = Arrhenius(A=1.0e13, b=0, E=(7.3, 'kcal/mol'))
    rate_coeff = Arrhenius(1.0e13, 0, (7.3, 'kcal/mol'))

Note: the usage of ``n`` as the temperature exponent has been deprecated. It is
still available in version 2.2 but will be removed.

As a shorthand, if the ``rate_coeff`` field is assigned a sequence of three numbers, these are assumed to be :math:`(A, b, E)` in the modified Arrhenius function::

    rate_coeff = [1.0e13, 0, (7.3, 'kcal/mol')] # equivalent to above

The units of the pre-exponential factor *A* can be specified explicitly if
desired. If not specified, they will be constructed using the *quantity*, *length*,
and *time* units specified in the units directive. Since the units of *A* depend on
the reaction order, the units of each reactant concentration (different for bulk
species in solution, surface species, and pure condensed-phase species), and the
units of the rate of progress (different for homogeneous and heterogeneous
reactions), it is usually best not to specify units for *A*, in which case they
will be computed taking all of these factors into account.

Note: if :math:`b \ne 0`, then the term :math:`T^b` should have units of
:math:`K^b`, which would change the units of *A*. This is not done, however, so
the units associated with A are really the units for :math:`k_f` . One way to
formally express this is to replace :math:`T^b` by the non-dimensional quantity
:math:`[T/(1 K)]^b`.

The ID String
-------------

An optional identifying string can be entered in the ``ID`` field, which can
then be used in the ``reactions`` field of a :class:`phase` or interface entry
to identify this reaction. If omitted, the reactions are assigned ID strings as
they are read in, beginning with ``'0001'``, ``'0002'``, etc.

Note that the ID string is only used when selectively importing reactions. If
all reactions in the local file or in an external one are imported into a phase
or interface, then the reaction ``ID`` field is not used.

.. _sec-reaction-options:

Options
-------

Certain conditions are normally flagged as errors by Cantera. In some cases,
they may not be errors, and the options field can be used to specify how they
should be handled.

``skip``
    The ``'skip'`` option can be used to temporarily remove this reaction from
    the phase or interface that imports it, just as if the reaction entry were
    commented out. The advantage of using skip instead of commenting it out is
    that a warning message is printed each time a phase or interface definition
    tries to import it. This serves as a reminder that this reaction is not
    included, which can easily be forgotten when a reaction is "temporarily"
    commented out of an input file.

``duplicate``
    Normally, when a reaction is imported into a phase, it is checked to see
    that it is not a duplicate of another reaction already present in the phase,
    and an error results if a duplicate is found. But in some cases, it may be
    appropriate to include duplicate reactions, for example if a reaction can
    proceed through two distinctly different pathways, each with its own rate
    expression.  Another case where duplicate reactions can be used is if it is
    desired to implement a reaction rate coefficient of the form:

    .. math::

        k_f(T) = \sum_{n=1}^{N} A_n T^{b_n} exp(-E_n/\hat{R}T)

    While Cantera does not provide such a form for reaction rates, it can be
    implemented by defining *N* duplicate reactions, and assigning one rate
    coefficient in the sum to each reaction.  If the ``'duplicate'`` option is
    specified, then the reaction not only *may* have a duplicate, it *must*. Any
    reaction that specifies that it is a duplicate, but cannot be paired with
    another reaction in the phase that qualifies as its duplicate generates an
    error.

``negative_A``
    If some of the terms in the above sum have negative :math:`A_n`, this scheme
    fails, since Cantera normally does not allow negative pre-exponential
    factors. But if there are duplicate reactions such that the total rate is
    positive, then negative *A* parameters are acceptable, as long as the
    ``'negative_A'`` option is specified.

``negative_orders``
    Reaction orders are normally required to be non-negative, since negative
    orders are non-physical and undefined at zero concentration. Cantera allows
    negative orders for a global reaction only if the ``negative_orders``
    override option is specified for the reaction.


Reactions with Pressure-Independent Rate
========================================

The :class:`reaction` entry is used to represent homogeneous reactions with
pressure-independent rate coefficients and mass action kinetics.  Examples of
reaction entries that implement some reactions in the GRI-Mech 3.0 natural gas
combustion mechanism [#Smith1997]_ are shown below::

    units(length = 'cm', quantity = 'mol', act_energy = 'cal/mol')
    ...
    reaction( "O + H2 <=> H + OH", [3.87000E+04, 2.7, 6260])
    reaction( "O + HO2 <=> OH + O2", [2.00000E+13, 0.0, 0])
    reaction( "O + H2O2 <=> OH + HO2", [9.63000E+06, 2.0, 4000])
    reaction( "O + HCCO <=> H + 2 CO", [1.00000E+14, 0.0, 0])
    reaction( "H + O2 + AR <=> HO2 + AR", kf=Arrhenius(A=7.00000E+17, b=-0.8, E=0))
    reaction( equation = "HO2 + C3H7 <=> O2 + C3H8", kf=Arrhenius(2.55000E+10, 0.255, -943))
    reaction( equation = "HO2 + C3H7 => OH + C2H5 + CH2O", kf=[2.41000E+13, 0.0, 0])

Three-Body Reactions
====================

A three-body reaction is a gas-phase reaction of the form:

.. math::

    {\rm A + B + M} \rightleftharpoons {\rm AB + M}

Here *M* is an unspecified collision partner that carries away excess energy to
stabilize the *AB* molecule (forward direction) or supplies energy to break the *AB*
bond (reverse direction).

Different species may be more or less effective in acting as the collision partner. A species that is much lighter than
*A* and *B* may not be able to transfer much of its kinetic energy, and so would be inefficient as a collision partner. On
the other hand, a species with a transition from its ground state that is nearly resonant with one in the *AB** activated
complex may be much more effective at exchanging energy than would otherwise be expected.

These effects can be accounted for by defining a collision efficiency
:math:`\epsilon` for each species, defined such that the forward reaction rate is

.. math::

    k_f(T)[A][B][M]

where

.. math::

    [M] = \sum_k \epsilon_k C_k

where :math:`C_k` is the concentration of species *k*. Since any constant
collision efficiency can be absorbed into the rate coefficient :math:`k_f(T)`, the
default collision efficiency is 1.0.

A three-body reaction may be defined using the :class:`three_body_reaction` entry. The equation string for a three-body
reaction must contain an ``'M'`` or ``'m'`` on both the reactant and product sides of the equation. The collision
efficiencies are specified as a string, with the species name followed by a colon and the efficiency.

Some examples from GRI-Mech 3.0 are shown below::

    three_body_reaction( "2 O + M <=> O2 + M", [1.20000E+17, -1, 0],
                         " AR:0.83 C2H6:3 CH4:2 CO:1.75 CO2:3.6 H2:2.4 H2O:15.4 ")

    three_body_reaction( "O + H + M <=> OH + M", [5.00000E+17, -1, 0],
                         efficiencies = " AR:0.7 C2H6:3 CH4:2 CO:1.5 CO2:2 H2:2 H2O:6 ")

    three_body_reaction(
        equation = "H + OH + M <=> H2O + M",
        rate_coeff = [2.20000E+22, -2, 0],
        efficiencies = " AR:0.38 C2H6:3 CH4:2 H2:0.73 H2O:3.65 "
    )

As always, the field names are optional *if* the field values are entered in the
declaration order.

Falloff Reactions
=================

A *falloff reaction* is one that has a rate that is first-order in [M] at low
pressure, like a three-body reaction, but becomes zero-order in [M] as [M]
increases. Dissociation / association reactions of polyatomic molecules often
exhibit this behavior.

The simplest expression for the rate coefficient for a falloff reaction is the
Lindemann form [#Lindemann1922]_:

.. math::

    k_f(T, [{\rm M}]) = \frac{k_0[{\rm M}]}{1 + \frac{k_0{\rm [M]}}{k_\infty}}

In the low-pressure limit, this approaches :math:`k0{\rm [M]}`, and in the
high-pressure limit it approaches :math:`k_\infty`.

Defining the non-dimensional reduced pressure:

.. math::

    P_r = \frac{k_0 {\rm [M]}}{k_\infty}

The rate constant may be written as

.. math::

    k_f(T, P_r) = k_\infty \left(\frac{P_r}{1 + P_r}\right)

More accurate models for unimolecular processes lead to other, more complex,
forms for the dependence on reduced pressure. These can be accounted for by
multiplying the Lindemann expression by a function :math:`F(T, P_r)`:

.. math::

    k_f(T, P_r) = k_\infty \left(\frac{P_r}{1 + P_r}\right) F(T, P_r)

This expression is used to compute the rate coefficient for falloff
reactions. The function :math:`F(T, P_r)` is the *falloff function*, and is
specified by assigning an embedded entry to the ``falloff`` field.

The Troe Falloff Function
-------------------------

A widely-used falloff function is the one proposed by Gilbert et
al. [#Gilbert1983]_:

.. math::

    \log_{10} F(T, P_r) = \frac{\log_{10} F_{cent}(T)}{1 + f_1^2}

    F_{cent}(T) = (1-A) \exp(-T/T_3) + A \exp (-T/T_1) + \exp(-T_2/T)

    f_1 = (\log_{10} P_r + C) / (N - 0.14 (\log_{10} P_r + C))

    C = -0.4 - 0.67\; \log_{10} F_{cent}

    N = 0.75 - 1.27\; \log_{10} F_{cent}

The :class:`Troe` directive requires specifying the first three parameters
:math:`(A, T_3, T_1)`. The fourth parameter, :math:`T_2`, is optional, defaulting to 0.0.

.. _sec-sri-falloff:

The SRI Falloff Function
------------------------

This falloff function is based on the one originally due to Stewart et
al. [#Stewart1989]_, which required three parameters :math:`(a, b, c)`. Kee et
al. [#Kee1989]_ generalized this function slightly by adding two more parameters
:math:`(d, e)`. (The original form corresponds to :math:`d = 1, e = 0`.) Cantera
supports the extended 5-parameter form, given by:

.. math::

    F(T, P_r) = d \bigl[a \exp(-b/T) + \exp(-T/c)\bigr]^{1/(1+\log_{10}^2 P_r )} T^e

In keeping with the nomenclature of Kee et al. [#Kee1989]_, we will refer to this as
the "SRI" falloff function. It is implemented by the :class:`SRI` directive.

.. :: NOTE: "definingphases.pdf" contains documentation for the Wang-Frenklach falloff
      function, which has a C++ implementation, but doesn't appear to be implemented
      in the CTI or CTML parsers.

Chemically-Activated Reactions
==============================

For these reactions, the rate falls off as the pressure increases, due to
collisional stabilization of a reaction intermediate. Example:

.. math::
     \mathrm{Si + SiH_4 (+M) \leftrightarrow Si_2H_2 + H_2 (+M)}

which competes with:

.. math::
    \mathrm{Si + SiH_4 (+M) \leftrightarrow Si_2H_4 (+M)}

Like falloff reactions, chemically-activated reactions are described by
blending between a "low pressure" and a "high pressure" rate expression. The
difference is that the forward rate constant is written as being proportional
to the *low pressure* rate constant:

.. math::

    k_f(T, P_r) = k_0 \left(\frac{1}{1 + P_r}\right) F(T, P_r)

and the optional blending function *F* may described by any of the
parameterizations allowed for falloff reactions. Chemically-activated
reactions can be defined using the :class:`chemically_activated_reaction`
directive.

An example of a reaction specified with this parameterization::

    chemically_activated_reaction('CH3 + OH (+ M) <=> CH2O + H2 (+ M)',
                                  kLow=[2.823201e+02, 1.46878, (-3270.56495, 'cal/mol')],
                                  kHigh=[5.880000e-14, 6.721, (-3022.227, 'cal/mol')],
                                  falloff=Troe(A=1.671, T3=434.782, T1=2934.21, T2=3919.0))

In this example, the units of :math:`k_0` (`kLow`) are m^3/kmol/s and the
units of :math:`k_\infty` (`kHigh`) are 1/s.

Pressure-Dependent Arrhenius Rate Expressions (P-Log)
=====================================================

The :class:`pdep_arrhenius` class represents pressure-dependent reaction rates
by logarithmically interpolating between Arrhenius rate expressions at various
pressures. Given two rate expressions at two specific pressures:

.. math::

    P_1: k_1(T) = A_1 T^{b_1} e^{E_1 / RT}

    P_2: k_2(T) = A_2 T^{b_2} e^{E_2 / RT}

The rate at an intermediate pressure :math:`P_1 < P < P_2` is computed as

.. math::

    \log k(T,P) = \log k_1(T) + \bigl(\log k_2(T) - \log k_1(T)\bigr)
        \frac{\log P - \log P_1}{\log P_2 - \log P_1}

Multiple rate expressions may be given at the same pressure, in which case the
rate used in the interpolation formula is the sum of all the rates given at that
pressure. For pressures outside the given range, the rate expression at the nearest
pressure is used.

An example of a reaction specified in this format::

    pdep_arrhenius('R1 + R2 <=> P1 + P2',
                   [(0.001315789, 'atm'), 2.440000e+10, 1.04, 3980.0],
                   [(0.039473684, 'atm'), 3.890000e+10, 0.989, 4114.0],
                   [(1.0, 'atm'), 3.460000e+12, 0.442, 5463.0],
                   [(10.0, 'atm'), 1.720000e+14, -0.01, 7134.0],
                   [(100.0, 'atm'), -7.410000e+30, -5.54, 12108.0],
                   [(100.0, 'atm'), 1.900000e+15, -0.29, 8306.0])

The first argument is the reaction equation. Each subsequent argument is a
sequence of four elements specifying a pressure and the Arrhenius parameters at
that pressure.

Chebyshev Reaction Rate Expressions
===================================

Class :class:`chebyshev_reaction` represents a phenomenological rate coefficient
:math:`k(T,P)` in terms of a bivariate Chebyshev polynomial. The rate constant
can be written as:

.. math:: \log k(T,P) = \sum_{t=1}^{N_T} \sum_{p=1}^{N_P} \alpha_{tp}
                            \phi_t(\tilde{T}) \phi_p(\tilde{P})

where :math:`\alpha_{tp}` are the constants defining the rate, :math:`\phi_n(x)`
is the Chebyshev polynomial of the first kind of degree :math:`n` evaluated at
:math:`x`, and

.. math::

    \tilde{T} \equiv \frac{2T^{-1} - T_\mathrm{min}^{-1} - T_\mathrm{max}^{-1}}
                          {T_\mathrm{max}^{-1} - T_\mathrm{min}^{-1}}

    \tilde{P} \equiv \frac{2 \log P - \log P_\mathrm{min} - \log P_\mathrm{max}}
                          {\log P_\mathrm{max} - \log P_\mathrm{min}}

are reduced temperature and reduced pressures which map the ranges
:math:`(T_\mathrm{min}, T_\mathrm{max})` and :math:`(P_\mathrm{min},
P_\mathrm{max})` to :math:`(-1, 1)`.

A Chebyshev rate expression is specified in terms of the coefficient matrix
:math:`\alpha` and the temperature and pressure ranges. An example of a
Chebyshev rate expression where :math:`N_T = 6` and :math:`N_P = 4` is::

    chebyshev_reaction('R1 + R2 <=> P1 + P2',
                       Tmin=290.0, Tmax=3000.0,
                       Pmin=(0.001, 'atm'), Pmax=(100.0, 'atm'),
                       coeffs=[[-1.44280e+01,  2.59970e-01, -2.24320e-02, -2.78700e-03],
                               [ 2.20630e+01,  4.88090e-01, -3.96430e-02, -5.48110e-03],
                               [-2.32940e-01,  4.01900e-01, -2.60730e-02, -5.04860e-03],
                               [-2.93660e-01,  2.85680e-01, -9.33730e-03, -4.01020e-03],
                               [-2.26210e-01,  1.69190e-01,  4.85810e-03, -2.38030e-03],
                               [-1.43220e-01,  7.71110e-02,  1.27080e-02, -6.41540e-04]])

Note that the Chebyshev polynomials are not defined outside the interval
:math:`(-1,1)`, and therefore extrapolation of rates outside the range of
temperatures and pressure for which they are defined is strongly discouraged.

Surface Reactions
=================

Heterogeneous reactions on surfaces are represented by an extended Arrhenius-
like rate expression, which combines the modified Arrhenius rate expression with
further corrections dependent on the fractional surface coverages
:math:`\theta_k` of one or more surface species. The forward rate constant for a
reaction of this type is:

.. math::

    k_f = A T^b \exp \left( - \frac{E_a}{RT} \right)
            \prod_k 10^{a_k \theta_k} \theta_k^{m_k}
            \exp \left( \frac{- E_k \theta_k}{RT} \right)

where :math:`A`, :math:`b`, and :math:`E_a` are the modified Arrhenius
parameters and :math:`a_k`, :math:`m_k`, and :math:`E_k` are the coverage
dependencies from species *k*. A reaction of this form with a single coverage
dependency (on the species ``H(S)``) can be written using class
:class:`surface_reaction` with the ``coverage`` keyword argument supplied to the
class :class:`Arrhenius`::

    surface_reaction("2 H(S) => H2 + 2 PT(S)",
                     Arrhenius(A, b, E_a,
                               coverage=['H(S)', a_1, m_1, E_1]))

For a reaction with multiple coverage dependencies, the following syntax is
used::

    surface_reaction("2 H(S) => H2 + 2 PT(S)",
                     Arrhenius(A, b, E_a,
                               coverage=[['H(S)', a_1, m_1, E_1],
                                         ['PT(S)', a_2, m_2, E_2]]))

Additional Options
==================

Reaction Orders
---------------

Explicit reaction orders different from the stoichiometric coefficients are
sometimes used for non-elementary reactions. For example, consider the global
reaction:

.. math::
    \mathrm{C_8H_{18} + 12.5 O_2 \rightarrow 8 CO_2 + 9 H_2O}

the forward rate constant might be given as [#Westbrook1981]_:

.. math::
    k_f = 4.6 \times 10^{11} [\mathrm{C_8H_{18}}]^{0.25} [\mathrm{O_2}]^{1.5}
          \exp\left(\frac{30.0\,\mathrm{kcal/mol}}{RT}\right)

This reaction could be defined as::

    reaction("C8H18 + 12.5 O2 => 8 CO2 + 9 H2O", [4.6e11, 0.0, 30.0],
             order="C8H18:0.25 O2:1.5")

Special care is required in this case since the units of the pre-exponential
factor depend on the sum of the reaction orders, which may not be an integer.

Normally, reaction orders are required to be positive. However, in some cases
negative reaction orders are found to be better fits for experimental data. In
these cases, the default behavior may be overridden by adding
``negative_orders`` to the reaction options, e.g.::

    reaction("C8H18 + 12.5 O2 => 8 CO2 + 9 H2O", [4.6e11, 0.0, 30.0],
             order="C8H18:-0.25 O2:1.75", options=['negative_orders'])


.. rubric:: References

.. [#Gilbert1983] R. G. Gilbert, K. Luther, and
   J. Troe. *Ber. Bunsenges. Phys. Chem.*, 87:169, 1983.

.. [#Lindemann1922] F. Lindemann. *Trans. Faraday Soc.*, 17:598, 1922.

.. [#Smith1997] Gregory P. Smith, David M. Golden, Michael Frenklach, Nigel
   W. Moriarty, Boris Eiteneer, Mikhail Goldenberg, C. Thomas Bowman, Ronald
   K. Hanson, Soonho Song, William C. Gardiner, Jr., Vitali V. Lissianski, , and
   Zhiwei Qin. GRI-Mech version 3.0, 1997. see
   http://www.me.berkeley.edu/gri_mech.

.. [#Stewart1989] P. H. Stewart, C. W. Larson, and D. Golden.
   *Combustion and Flame*, 75:25, 1989.

.. [#Kee1989] R. J. Kee, F. M. Rupley, and J. A. Miller. Chemkin-II: A Fortran
   chemical kinetics package for the analysis of gas-phase chemical
   kinetics. Technical Report SAND89-8009, Sandia National Laboratories, 1989.

.. [#Westbrook1981] C. K. Westbrook and F. L. Dryer. Simplified reaction
   mechanisms for the oxidation of hydrocarbon fuels in flames. *Combustion
   Science and Technology* **27**, pp. 31--43. 1981.
