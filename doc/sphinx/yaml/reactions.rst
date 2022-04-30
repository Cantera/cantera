.. highlight:: yaml

.. _sec-yaml-reactions:

*********
Reactions
*********

The fields common to all ``reaction`` entries are:

``equation``
    The stoichiometric equation for the reaction. Each term (that is,
    stoichiometric coefficient, species name, ``+`` or ``<=>``) in the equation
    must be separated by a space.

    Reversible reactions may be written using ``<=>`` or ``=`` to separate
    reactants and products. Irreversible reactions are written using ``=>``.

``type``
    A string specifying the type of reaction or rate coefficient
    parameterization. The default is ``elementary``. Reaction types are:

    - :ref:`elementary <sec-yaml-elementary>`
    - :ref:`three-body <sec-yaml-three-body>`
    - :ref:`Blowers-Masel <sec-yaml-Blowers-Masel>`
    - :ref:`two-temperature-plasma <sec-yaml-two-temperature-plasma>`
    - :ref:`falloff <sec-yaml-falloff>`
    - :ref:`chemically-activated <sec-yaml-chemically-activated>`
    - :ref:`pressure-dependent-Arrhenius <sec-yaml-pressure-dependent-Arrhenius>`
    - :ref:`Chebyshev <sec-yaml-Chebyshev>`

    Reactions without a specified ``type`` on surfaces or edges are
    automatically treated as :ref:`interface-Arrhenius <sec-yaml-interface-Arrhenius>`
    reactions, unless a ``sticking-coefficient`` implies a
    :ref:`sticking-Arrhenius <sec-yaml-sticking-Arrhenius>` reaction. Interface
    reactions that involve charge transfer between phases are automatically treated as
    :ref:`electrochemical <sec-yaml-electrochemical-reaction>` reactions.

    Reactions on surfaces or edges specifying ``type`` as ``Blowers-Masel`` are treated
    as :ref:`interface-Blowers-Masel <sec-yaml-interface-Blowers-Masel>` or
    :ref:`sticking-Blowers-Masel <sec-yaml-sticking-Blowers-Masel>`.

``duplicate``
    Boolean indicating whether the reaction is a known duplicate of another
    reaction. The default is ``false``.

``orders``
    An optional mapping of species to explicit reaction orders to use. Reaction
    orders for reactant species not explicitly mentioned are taken to be their
    respective stoichiometric coefficients. See
    `Reaction orders <https://cantera.org/science/kinetics.html#reaction-orders>`__
    for additional information.

``negative-orders``
    Boolean indicating whether negative reaction orders are allowed. The
    default is ``false``.

``nonreactant-orders``
    Boolean indicating whether orders for non-reactant species are allowed.
    The default is ``false``.

Depending on the reaction ``type``, other fields may be necessary to specify
the rate of the reaction.


Reaction rate expressions
=========================

.. _sec-yaml-Arrhenius-rate:

Arrhenius
---------

Arrhenius rate expressions are specified as a mapping with fields:

``A``
    The pre-exponential factor :math:`A`
``b``
    The temperature exponent :math:`b`
``Ea``
    The activation energy :math:`E_a`

or a corresponding three-element list. The following are equivalent::

    {A: -2.70000E+13 cm^3/mol/s, b: 0, Ea: 355 cal/mol}
    [-2.70000E+13 cm^3/mol/s, 0, 355 cal/mol]


.. _sec-yaml-Blowers-Masel-rate:

Blowers-Masel
-------------

Blowers-Masel rate expressions calculate the rate constant based on the Blowers Masel
approximation as
`described here <https://cantera.org/science/kinetics.html#sec-blowers-masel>`__.
The rate parameters are specified as a mapping with fields:

``A``
    The pre-exponential factor :math:`A`
``b``
    The temperature exponent :math:`b`
``Ea0``
    The intrinsic activation energy :math:`E_{a0}`
``w``
    The average of the bond dissociation energy of the bond breaking and that being
    formed in the reaction :math:`w`

or a corresponding four-element list. The following are equivalent::

    {A: 3.87e+04 cm^3/mol/s, b: 2.7, Ea0: 6260.0 cal/mol, w: 1e9 cal/mol}
    [3.87e+04 cm^3/mol/s, 2.7, 6260.0 cal/mol, 1e9 cal/mol]


.. _sec-yaml-two-temperature-plasma-rate:

Two-Temperature Plasma
----------------------

Two-temperature plasma reactions involve an electron as one of the reactants, where the
electron temperature may differ from the gas temperature as
`described here <https://cantera.org/science/kinetics.html#two-temperature-plasma-reactions>`__.
The rate parameters are specified as a mapping with fields:

``A``
    The pre-exponential factor
``b``
    The temperature exponent, which is applied to the electron temperature
``Ea-gas``
    The activation energy term :math:`E_{a,g}` that is related to the gas temperature
``Ea-electron``
    The activation energy term :math:`E_{a,e}` that is related to the electron
    temperature

or a corresponding four-element list. The following are equivalent::

    {A: 17283, b: -3.1, Ea-gas: -5820 J/mol, Ea-electron: 1081 J/mol}
    [17283, -3.1, -5820 J/mol, 1081 J/mol]


.. _sec-yaml-efficiencies:

Efficiencies
============

Some reaction types include parameters for the "efficiency" of different species
as third-body colliders. For these reactions, the following additional fields
are supported:

``efficiencies``
    A mapping of species names to efficiency values

``default-efficiency``
    The efficiency for use for species not included in the ``efficiencies``
    mapping. Defaults to 1.0.


Reaction types
==============

.. _sec-yaml-elementary:

``elementary``
--------------

A homogeneous reaction with a pressure-independent rate coefficient and mass
action kinetics, as
`described here <https://cantera.org/science/kinetics.html#reactions-with-a-pressure-independent-rate>`__.

Additional fields are:

``rate-constant``
    An :ref:`Arrhenius-type <sec-yaml-Arrhenius-rate>` list or mapping.

``negative-A``
    A boolean indicating whether a negative value for the pre-exponential factor
    is allowed. The default is ``false``.

Example::

    equation: N + NO <=> N2 + O
    rate-constant: {A: -2.70000E+13 cm^3/mol/s, b: 0, Ea: 355 cal/mol}
    negative-A: true


.. _sec-yaml-three-body:

``three-body``
--------------

A three body reaction as
`described here <https://cantera.org/science/kinetics.html#three-body-reactions>`__.

The reaction equation should include the third body collision partner ``M``.

Includes the fields of an ``elementary`` reaction, plus the fields for
specifying :ref:`efficiencies <sec-yaml-efficiencies>`.

Example::

    equation: 2 O + M = O2 + M
    type: three-body
    rate-constant: [1.20000E+17 cm^6/mol^2/s, -1, 0]
    efficiencies: {AR: 0.83, H2O: 5}


.. _sec-yaml-Blowers-Masel:

``Blowers-Masel``
-----------------

Includes the fields of an :ref:`elemntary <sec-yaml-elementary>` reaction, except that
the ``rate-constant`` field is a
:ref:`Blowers-Masel-type <sec-yaml-Blowers-Masel-rate>` list or mapping.

Example::

    equation: O + H2 <=> H + OH
    type: Blowers-Masel
    rate-constant: {A: 3.87e+04 cm^2/mol/s, b: 2.7, Ea0: 6260.0 cal/mol, w: 1e9 cal/mol}


.. _sec-yaml-two-temperature-plasma:

``two-temperature-plasma``
--------------------------

Includes the fields of an :ref:`elementary <sec-yaml-elementary>` reaction, except that
the ``rate-constant`` field is a
:ref:`Two-temperature-plasma-type <sec-yaml-two-temperature-plasma-rate>` list or
mapping.

Example::

    equation: O + H => O + H
    type: two-temperature-plasma
    rate-constant: {A: 17283, b: -3.1, Ea-gas: -5820 J/mol, Ea-electron: 1081 J/mol}


.. _sec-yaml-falloff:

``falloff``
-----------

A falloff reaction as
`described here <https://cantera.org/science/kinetics.html#falloff-reactions>`__.

The reaction equation should include the pressure-dependent third body collision
partner ``(+M)`` or ``(+name)`` where ``name`` is the name of a species. The
latter case is equivalent to setting the efficiency for ``name`` to 1 and the
efficiency for all other species to 0.

Includes field for specifying :ref:`efficiencies <sec-yaml-efficiencies>` as well
as:

``high-P-rate-constant``
    An :ref:`sec-yaml-Arrhenius-rate` expression for the high-pressure limit

``low-P-rate-constant``
    An :ref:`sec-yaml-Arrhenius-rate` expression for the low-pressure limit

``Troe``
    Parameters for the
    `Troe <https://cantera.org/science/kinetics.html#the-troe-falloff-function>`__
    falloff function. A mapping containing the keys ``A``, ``T3``, ``T1`` and
    optionally ``T2``. The default value for ``T2`` is 0.

``SRI``
    Parameters for the
    `SRI <https://cantera.org/science/kinetics.html#the-sri-falloff-function>`__
    falloff function. A mapping containing the keys ``A``, ``B``, ``C``, and
    optionally ``D`` and ``E``. The default values for ``D`` and ``E`` are 1.0
    and 0.0, respectively.

``Tsang``
    Parameters for the
    `Tsang <https://cantera.org/science/kinetics.html#tsang-s-approximation-to-f-cent>`__
    falloff function. A mapping containing the keys ``A`` and ``B``. The default value
    for ``B`` is 0.0.

Example::

    equation: H + CH2 (+ N2) <=> CH3 (+N2)
    type: falloff
    high-P-rate-constant: [6.00000E+14 cm^3/mol/s, 0, 0]
    low-P-rate-constant: {A: 1.04000E+26 cm^6/mol^2/s, b: -2.76, Ea: 1600}
    Troe: {A: 0.562, T3: 91, T1: 5836}


.. _sec-yaml-chemically-activated:

``chemically-activated``
------------------------

A chemically activated reaction as
`described here <https://cantera.org/science/kinetics.html#chemically-activated-reactions>`__.

The parameters are the same as for :ref:`sec-yaml-falloff` reactions.

Example::

    equation: CH3 + OH (+M) <=> CH2O + H2 (+M)
    type: chemically-activated
    high-P-rate-constant: [5.88E-14, 6.721, -3022.227]
    low-P-rate-constant: [282320.078, 1.46878, -3270.56495]


.. _sec-yaml-pressure-dependent-Arrhenius:

``pressure-dependent-Arrhenius``
--------------------------------

A pressure-dependent reaction using multiple Arrhenius expressions as
`described here <https://cantera.org/science/kinetics.html#pressure-dependent-arrhenius-rate-expressions-p-log>`__.

The only additional field in this reaction type is:

``rate-constants``
    A list of mappings, where each mapping is the mapping form of an
    :ref:`sec-yaml-Arrhenius-rate` expression with the addition of a pressure ``P``.

Example::

    equation: H + CH4 <=> H2 + CH3
    type: pressure-dependent-Arrhenius
    rate-constants:
    - {P: 0.039474 atm, A: 2.720000e+09 cm^3/mol/s, b: 1.2, Ea: 6834.0}
    - {P: 1.0 atm, A: 1.260000e+20, b: -1.83, Ea: 15003.0}
    - {P: 1.0 atm, A: 1.230000e+04, b: 2.68, Ea: 6335.0}
    - {P: 1.01325 MPa, A: 1.680000e+16, b: -0.6, Ea: 14754.0}


.. _sec-yaml-Chebyshev:

``Chebyshev``
-------------

A reaction parameterized as a bivariate Chebyshev polynomial as
`described here <https://cantera.org/science/kinetics.html#chebyshev-reaction-rate-expressions>`__.

Additional fields are:

``temperature-range``
    A list of two values specifying the minimum and maximum temperatures at
    which the rate constant is valid

``pressure-range``
    A list of two values specifying the minimum and maximum pressures at
    which the rate constant is valid

``data``
    A list of lists containing the Chebyshev coefficients

Example::

    equation: CH4 <=> CH3 + H
    type: Chebyshev
    temperature-range: [290, 3000]
    pressure-range: [0.0098692326671601278 atm, 98.692326671601279 atm]
    data: [[-1.44280e+01,  2.59970e-01, -2.24320e-02, -2.78700e-03],
           [ 2.20630e+01,  4.88090e-01, -3.96430e-02, -5.48110e-03],
           [-2.32940e-01,  4.01900e-01, -2.60730e-02, -5.04860e-03],
           [-2.93660e-01,  2.85680e-01, -9.33730e-03, -4.01020e-03],
           [-2.26210e-01,  1.69190e-01,  4.85810e-03, -2.38030e-03],
           [-1.43220e-01,  7.71110e-02,  1.27080e-02, -6.41540e-04]]


.. _sec-yaml-interface-Arrhenius:

``interface-Arrhenius``
-----------------------

A reaction occurring on a surface between two bulk phases, or along an edge
at the intersection of two surfaces, as
`described here <https://cantera.org/science/kinetics.html#sec-surface>`__.

Includes the fields of an :ref:`sec-yaml-elementary` reaction plus:

``coverage-dependencies``
    A mapping of species names to coverage dependence parameters, where these
    parameters are contained in either a mapping with the fields:

    ``a``
        Coefficient for exponential dependence on the coverage

    ``m``
        Power-law exponent of coverage dependence

    ``E``
        Activation energy dependence on coverage, which uses the same sign convention
        as the leading-order activation energy term

    or a list containing the three elements above, in the given order.

    Note that parameters ``a``, ``m`` and ``E`` correspond to parameters
    :math:`\eta_{ki}`, :math:`\mu_{ki}` and :math:`\epsilon_{ki}` in Eq 11.113 of
    [Kee, R. J., Coltrin, M. E., & Glarborg, P.(2003). Chemically reacting flow:
    theory and practice. John Wiley & Sons], respectively.

Examples::

    - equation: 2 H(s) => H2 + 2 Pt(s)
      rate-constant: {A: 3.7e21 cm^2/mol/s, b: 0, Ea: 67400 J/mol}
      coverage-dependencies: {H(s): {a: 0, m: 0, E: -6000 J/mol}}

    - equation: CH4 + PT(S) + O(S) => CH3(S) + OH(S)
      rate-constant: {A: 5.0e+18, b: 0.7, Ea: 4.2e+04}
      coverage-dependencies:
        O(S): [0, 0, 8000]
        PT(S): [0, -1.0, 0]


.. _sec-yaml-interface-Blowers-Masel:

``interface-Blowers-Masel``
---------------------------

Includes the same fields as :ref:`interface-Arrhenius <sec-yaml-interface-Arrhenius>`,
while using the :ref:`Blowers-Masel <sec-yaml-Blowers-Masel-rate>` parameterization
for the rate constant.

Example::

    equation: 2 H(s) => H2 + 2 Pt(s)
    type: Blowers-Masel
    rate-constant: {A: 3.7e21 cm^2/mol/s, b: 0, Ea0: 67400 J/mol, w: 1000000 J/mol}
    coverage-dependencies: {H(s): {a: 0, m: 0, E: -6000 J/mol}}


.. _sec-yaml-sticking-Arrhenius:

``sticking-Arrhenius``
----------------------

A sticking reaction occurring on a surface adjacent to a bulk phase, as
`described here <https://cantera.org/science/kinetics.html#sec-sticking>`__.

Includes the fields of an :ref:`sec-yaml-interface-Arrhenius` reaction plus:

``sticking-coefficient``
    An :ref:`Arrhenius-type <sec-yaml-Arrhenius-rate>` expression for the sticking
    coefficient

``Motz-Wise``
    A boolean indicating whether to use the Motz-Wise correction factor for sticking
    coefficients near unity. Defaults to ``false``.

``sticking-species``
    The name of the sticking species. Required if the reaction includes multiple
    non-surface species.

Example::

    equation: OH + PT(S) => OH(S)
    sticking-coefficient: {A: 1.0, b: 0, Ea: 0}


.. _sec-yaml-sticking-Blowers-Masel:

``sticking-Blowers-Masel``
--------------------------

Includes the same fields as :ref:`sticking-Arrhenius <sec-yaml-sticking-Arrhenius>`,
while using the :ref:`Blowers-Masel <sec-yaml-Blowers-Masel-rate>` parameterization
for the sticking coefficient.

Example::

    equation: OH + PT(S) => OH(S)
    type: Blowers-Masel
    sticking-coefficient: {A: 1.0, b: 0, Ea0: 0, w: 100000}
    Motz-Wise: true


.. _sec-yaml-electrochemical-reaction:

``electrochemical``
-------------------

Interface reactions involving charge transfer between phases,
as `described here <https://cantera.org/documentation/dev/doxygen/html/d6/ddd/classCantera_1_1ElectrochemicalReaction.html#details>`__.

Includes the fields of an :ref:`sec-yaml-interface-Arrhenius` reaction, plus:

``beta``
    The symmetry factor for the reaction. Default is 0.5.

``exchange-current-density-formulation``
    Set to ``true`` if the rate constant parameterizes the exchange current
    density. Default is ``false``.

Example::

    equation: LiC6 <=> Li+(e) + C6
    rate-constant: [5.74, 0.0, 0.0]
    beta: 0.4
