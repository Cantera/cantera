.. py:currentmodule:: ctml_writer

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

    k_f(T) = A T^n \exp(-E/\hat{R}T)

whish is defined with an :class:`Arrhenius` entry::

    rate_coeff = Arrhenius(A=1.0e13, n=0, E=(7.3, 'kcal/mol'))
    rate_coeff = Arrhenius(1.0e13, 0, (7.3, 'kcal/mol'))

As a shorthand, if the ``rate_coeff`` field is assigned a sequence of three numbers, these are assumed to be :math:`(A, n, E)` in the modified Arrhenius function::

    rate_coeff = [1.0e13, 0, (7.3, 'kcal/mol')] # equivalent to above

The units of the pre-exponential factor *A* can be specified explicitly if
desired. If not specified, they will be constructed using the *quantity*, length,
and time units specified in the units directive. Since the units of *A* depend on
the reaction order, the units of each reactant concentration (different for bulk
species in solution, surface species, and pure condensed-phase species), and the
units of the rate of progress (different for homogeneous and heterogeneous
reactions), it is usually best not to specify units for *A*, in which case they
will be computed taking all of these factors into account.

Note: if :math:`n \ne 0`, then the term :math:`T^n` should have units of
:math:`K^n`, which would change the units of *A*. This is not done, however, so
the units associated with A are really the units for :math:`k_f` . One way to
formally express this is to replace :math:`T^n` by the non-dimensional quantity
:math:`[T/(1 K)]^n`.

The ID String
-------------

An optional identifying string can be entered in the ``ID`` field, which can
then be used in the ``reactions`` field of a :class:`phase` or interface entry
to identify this reaction. If omitted, the reactions are assigned ID strings as
they are read in, beginning with ``'0001'``, ``'0002'``, etc.

Note that the ID string is only used when selectively importing reactions. If
all reactions in the local file or in an external one are imported into a phase
or interface, then the reaction ``ID`` field is not used.

Options
-------

Certain conditions are normally flagged as errors by Cantera. In some cases,
theey may not be errors, and the options field can be used to specify how they
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

Reactions with Pressure-Independent Rate
========================================

The :class:`reaction` entry is used to represent homogeneous reactions with
pressure-independent rate coefficients and mass action kinetics.  Examples of
reaction entries that implement some reactions in the GRI-Mech 3.0 natural gas
combustion mechanism [Smith et al., 1997] are shown below::

    units(length = 'cm', quantity = 'mol', act_energy = 'cal/mol')
    ...
    reaction( "O + H2 <=> H + OH", [3.87000E+04, 2.7, 6260])
    reaction( "O + HO2 <=> OH + O2", [2.00000E+13, 0.0, 0])
    reaction( "O + H2O2 <=> OH + HO2", [9.63000E+06, 2.0, 4000])
    reaction( "O + HCCO <=> H + 2 CO", [1.00000E+14, 0.0, 0])
    reaction( "H + O2 + AR <=> HO2 + AR", [7.00000E+17, -0.8, 0])
    reaction( "HO2 + C3H7 <=> O2 + C3H8", [2.55000E+10, 0.255, -943])
    reaction( "HO2 + C3H7 => OH + C2H5 + CH2O", [2.41000E+13, 0.0, 0])

Three-Body Reactions
====================

A three-body reaction is a gas-phase reaction of the form:

.. math::

    {\rm A + B} \rightleftharpoons {\rm AB + M}

Here M is an unspecified collision partner that carries away excess energy to
stabilize the AB molecule (forward direction) or supplies energy to break the AB
bond (reverse direction).

Different species may be more or less effective in acting as the collision partner. A species that is much lighter than
A and B may not be able to transfer much of its kinetic energy, and so would be inefficient as a collision partner. On
the other hand, a species with a transition from its ground state that is nearly resonant with one in the AB* activated
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
reaction must contain an ``'M'`` or ``'m'`` on both the reactant and product sides of the equation.

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
Lindemann form [Lindemann, 1922]:

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

A widely-used falloff function is the one proposed by Gilbert et al. [1983]:

.. math::

    \log_{10} F(T, P_r) = \frac{\log_{10} F_{cent}(T)}{1 + f_1^2}

    F_{cent}(T) = (1-A) \exp(-T/T_3) + A \exp (-T/T_1) + \exp(-T_2/T)

    f_1 = (\log_{10} P_r + C) / (N - 0.14 (\log_{10} P_r + C))

    C = -0.4 - 0.67\; \log_{10} F_{cent}

    N = 0.75 - 1.27\; \log_{10} F_{cent}

The :class:`Troe` directive requires specifying the first three parameters
:math:`(A, T_3, T_1)`. The fourth paramteter, :math:`T_2`, is optional, defaulting to 0.0.

The SRI Falloff Function
------------------------

This falloff function is based on the one originally due to Stewart et
al. [1989], which required three parameters :math:`(a, b, c)`. Kee et al. [1989]
generalized this function slightly by adding two more parameters :math:`(d,
e)`. (The original form corresponds to :math:`d = 1, e = 0`.) Cantera supports
the extended 5-parameter form, given by:

.. math::

    F(T, P_r) = d \bigl[a \exp(-b/T) + \exp(-T/c)\bigr]^{1/(1+\log_{10}^2 P_r )} T^e

In keeping with the nomenclature of [Kee et al., 1989], we will refer to this as
the "SRI" falloff function. It is implemented by the :class:`SRI` directive.

.. :: NOTE: "definingphases.pdf" contains documentation for the Wang-Frenklach falloff
      function, which has a C++ implementation, but doesn't appear to be implemented
      in the CTI or CTML parsers.

