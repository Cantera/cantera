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

    - ``elementary`` (:ref:`details <sec-yaml-elementary>`)
    - ``three-body`` (:ref:`details <sec-yaml-three-body>`)
    - ``Blowers-Masel`` (:ref:`details <sec-yaml-Blowers-Masel>`)
    - ``two-temperature-plasma`` (:ref:`details <sec-yaml-two-temperature-plasma>`)
    - ``falloff`` (:ref:`details <sec-yaml-falloff>`)
    - ``chemically-activated`` (:ref:`details <sec-yaml-chemically-activated>`)
    - ``pressure-dependent-Arrhenius`` (:ref:`details <sec-yaml-pressure-dependent-Arrhenius>`)
    - ``Chebyshev`` (:ref:`details <sec-yaml-Chebyshev>`)

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
    respective stoichiometric coefficients. See :ref:`sec-reaction-orders` for
    additional information.

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
approximation as :ref:`described here <sec-blowers-masel>`. The rate parameters are
specified as a mapping with fields:

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
:ref:`described here <sec-two-temperature-plasma-rate>`.
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

.. _sec-yaml-rate-types:

Reaction types
==============

.. _sec-yaml-elementary:

``elementary``
--------------

A homogeneous reaction with a pressure-independent rate coefficient and mass action
kinetics, as :ref:`described here <sec-arrhenius-rate>`.

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

A three body reaction as :ref:`described here <sec-three-body-reaction>`.

The reaction equation must include a third body collision partner, which may be either
a specific species or the generic third body ``M``.

Includes the fields of an ``elementary`` reaction, plus the fields for
specifying :ref:`efficiencies <sec-yaml-efficiencies>`.

Example::

    equation: 2 O + M = O2 + M
    type: three-body
    rate-constant: [1.20000E+17 cm^6/mol^2/s, -1, 0]
    efficiencies: {AR: 0.83, H2O: 5}

The ``type`` field of the YAML entry may be omitted. Reactions containing the generic
third body M are automatically identified as three-body reactions. Reactions are also
identified as three-body reactions if all of the following conditions are met:

- There is exactly one species appearing as both a reactant and product
- All reactants and products have integral stoichiometric coefficients
- The sum of the stoichiometric coefficients for either the reactants or products is 3.

Examples::

    - equation: H + 2 O2 <=> HO2 + O2  # Reaction 34
      rate-constant: {A: 2.08e+19, b: -1.24, Ea: 0.0}
    - equation: H + O2 + N2 <=> HO2 + N2  # Reaction 36
      rate-constant: {A: 2.6e+19, b: -1.24, Ea: 0.0}

.. caution::
    If a corresponding reaction with the generic third body M also appears in the
    mechanism, such as::

        - equation: H + O2 + M <=> HO2 + M  # Reaction 33
          rate-constant: {A: 2.8e+18, b: -0.86, Ea: 0.0}
          efficiencies: {O2: 0.0, H2O: 0.0, CO: 0.75, CO2: 1.5, C2H6: 1.5, N2: 0.0, AR: 0.0}

    then the third body efficiency for any third bodies that are given in the explicit
    form of Reaction 34 or Reaction 35 above must be set to zero, as shown here for
    O2 and N2, or the reactions must be marked as duplicate.

.. versionchanged:: 3.0

    Three body reactions are detected automatically and the the ``type`` field may be
    omitted. Reactions with explicit third bodies are required to be marked as
    duplicates of reactions with the generic third body if the corresponding efficiency
    is not zero.

.. versionadded:: 3.1

    Reactions with explicit third bodies and the corresponding reaction with "M" issue
    warnings instead of raising errors by default. The
    :ref:`sec-yaml-phase-explicit-third-body-duplicates` field of the phase entry can be
    used to control how these reactions are handled.


.. _sec-yaml-Blowers-Masel:

``Blowers-Masel``
-----------------

Includes the fields of an :ref:`elementary <sec-yaml-elementary>` reaction, except that
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
:ref:`Two-temperature-plasma <sec-yaml-two-temperature-plasma-rate>` list or
mapping.

Example::

    equation: O + H => O + H
    type: two-temperature-plasma
    rate-constant: {A: 17283, b: -3.1, Ea-gas: -5820 J/mol, Ea-electron: 1081 J/mol}


.. _sec-yaml-falloff:

``falloff``
-----------

A falloff reaction as :ref:`described here <sec-falloff-rate>`.

The reaction equation should include the pressure-dependent third body collision
partner ``(+M)`` or ``(+name)`` where ``name`` is the name of a species. The
latter case is equivalent to setting the efficiency for ``name`` to 1 and the
efficiency for all other species to 0.

Includes field for specifying :ref:`efficiencies <sec-yaml-efficiencies>` as well as:

``high-P-rate-constant``
    An :ref:`sec-yaml-Arrhenius-rate` expression for the high-pressure limit

``low-P-rate-constant``
    An :ref:`sec-yaml-Arrhenius-rate` expression for the low-pressure limit

``Troe``
    Parameters for the :ref:`Troe <sec-troe-falloff>` falloff function. A mapping
    containing the keys ``A``, ``T3``, ``T1`` and optionally ``T2``. The default value
    for ``T2`` is 0.

``SRI``
    Parameters for the :ref:`SRI <sec-sri-falloff>` falloff function. A mapping
    containing the keys ``A``, ``B``, ``C``, and optionally ``D`` and ``E``. The default
    values for ``D`` and ``E`` are 1.0 and 0.0, respectively.

``Tsang``
    Parameters for the :ref:`Tsang <sec-tsang-falloff>` falloff function. A mapping
    containing the keys ``A`` and ``B``. The default value for ``B`` is 0.0.

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
:ref:`described here <sec-chemically-activated-rate>`.

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
:ref:`described here <sec-plog-rate>`.

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
:ref:`described here <sec-chebyshev-rate>`.

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


.. _sec-yaml-linear-burke:

``linear-burke``
-------------

A complex-forming reaction (one that depends on both P and X) parameterized
according to the reduced-pressure linear mixture rule as
:ref:`described here <sec-linear-burke>`.

The pressure-dependent aspect of the rate constant can be parameterized in the user's choice of
:ref:`Troe <sec-yaml-falloff>`,
:ref:`pressure-dependent-arrhenius <sec-yaml-pressure-dependent-Arrhenius>`, or
:ref:`Chebyshev <sec-yaml-Chebyshev>` representations. The same parameters used for a standalone
Troe, PLOG, or Chebyshev reaction are then inserted directly below ``eps`` or ``eig0`` for a given collider
(note: Troe cannot be given its own ``efficiencies`` key). At minimum, this treatment must be applied to ``"M"``.
However, additional colliders may also be given their own Troe, PLOG, or Chebyshev
parameterization if so desired. Mixing and matching of types within the same reaction is allowed (e.g., a PLOG
table for ``"M"``, Troe parameters for ``"H2"``, and Chebyshev data for ``"NH3"``).

Additional fields are:

``collider-list``
    A list of dictionaries, where each entry contains parameters corresponding
    to individual colliders (species in the bath gas).

``collider``
    The name of the collider species, which must be entered inside quotations (e.g.,
    ``"H2O"``). The first collider defined must be ``"M"``, which represents the generic
    reference collider (often ``Ar`` or ``N2``) that represents all species lacking their
    own explicit parameterization.

``eps`` or ``eig0``
    The fractional contribution of each bath gas component (collider) to the reduced
    pressure. ``eps`` represents the third-body efficiency of the collider relative
    to that of the reference collider ``"M"`` (``eps: {A:1, b:0, Ea: 0}`` must be provided for
    ``"M"`` by necessity). ``eig0`` represents the absolute value of the least negative chemically
    significant eigenvalue of the master equation, evaluated for a collider at its low-pressure
    limit. All explicitly defined colliders must include either ``eps`` or ``eig0``, but the choice
    must remain consistent throughout a single reaction (either all colliders are defined with ``eps``,
    or all are defined with ``eig0``). In both cases, the parameters are entered in Arrhenius format to
    enable representation of their temperature-dependence.

A mathematical description of this YAML implementation can be found in Eq. 8 of
:cite:t:`singal2024`. [CITATION NOT YET ADDED]

Examples::

    equation: H + OH <=> H2O
    type: linear-burke
    collider-list:
    - collider: 'M' # N2 is reference collider (Troe format)
      eps: {A: 1, b: 0, Ea: 0}
      low-P-rate-constant: {A: 4.530000e+21, b: -1.820309e+00, Ea: 4.987000e+02}
      high-P-rate-constant: {A: 2.510000e+13, b: 2.329303e-01, Ea: -1.142000e+02}
      Troe: {A: 9.995044e-01, T3: 1.0e-30, T1: 1.0e+30}
    - collider: 'AR'
      eps: {A: 2.20621e-02, b: 4.74036e-01, Ea: -1.13148e+02}
    - collider: 'H2O'
      eps: {A: 1.04529e-01, b: 5.50787e-01, Ea: -2.32675e+02}

    equation: H + O2 (+M) <=> HO2 (+M)  # Including "(+M)" is optional
    type: linear-burke
    collider-list: 
    - collider: "M" # Argon is reference collider (PLOG format)
      eps: {A: 1, b: 0, Ea: 0}
      rate-constants:
      - {P: 1.316e-02 atm, A: 9.39968e+14, b: -2.14348e+00, Ea: 7.72730e+01}
      - {P: 1.316e-01 atm, A: 1.07254e+16, b: -2.15999e+00, Ea: 1.30239e+02}
      - {P: 3.947e-01 atm, A: 3.17830e+16, b: -2.15813e+00, Ea: 1.66994e+02}
      - {P: 1.000e+00 atm, A: 7.72584e+16, b: -2.15195e+00, Ea: 2.13473e+02}
      - {P: 3.000e+00 atm, A: 2.11688e+17, b: -2.14062e+00, Ea: 2.79031e+02}
      - {P: 1.000e+01 atm, A: 6.53093e+17, b: -2.13213e+00, Ea: 3.87493e+02}
      - {P: 3.000e+01 atm, A: 1.49784e+18, b: -2.10026e+00, Ea: 4.87579e+02}
      - {P: 1.000e+02 atm, A: 3.82218e+18, b: -2.07057e+00, Ea: 6.65984e+02}
    - collider: "HE"
      eps: {A: 3.37601e-01, b: 1.82568e-01, Ea: 3.62408e+01}
    - collider: "N2"
      eps: {A: 1.24932e+02, b: -5.93263e-01, Ea: 5.40921e+02}
    - collider: "H2"
      eps: {A: 3.13717e+04, b: -1.25419e+00, Ea: 1.12924e+03}
    - collider: "CO2"
      eps: {A: 1.62413e+08, b: -2.27622e+00, Ea: 1.97023e+03}
    - collider: "NH3"
      eps: {A: 4.97750e+00, b: 1.64855e-01, Ea: -2.80351e+02}
    - collider: "H2O"
      eps: {A: 3.69146e+01, b: -7.12902e-02, Ea: 3.19087e+01}

    equation: H2O2 <=> 2 OH
    type: linear-burke
    collider-list:
    - collider: 'M' # Argon is reference collider (Chebyshev format)
      eps: {A: 1, b: 0, Ea: 0}
      temperature-range: [200.0, 2000.0]
      pressure-range: [1.000e-01 atm, 1.000e+02 atm]
      data:
      - [-1.5843e+01, 8.7088e-01, -9.4364e-02, -2.8099e-03, -4.4803e-04, 1.5809e-03, -2.5088e-04]
      - [2.3154e+01, 5.2739e-01, 2.8862e-02, -5.4601e-03, 7.0783e-04, -3.0282e-03, 7.8121e-04]
      - [-3.8008e-01, 8.6349e-02, 4.0292e-02, -7.2269e-03, 5.7570e-04, 2.7944e-03, -1.4912e-03]
      - [-1.4800e-01, -7.1798e-03, 2.2052e-02, 6.2269e-03, -5.9801e-03, -8.2205e-06, 1.9243e-03]
      - [-6.0604e-02, -1.4203e-02, 1.3414e-03, 9.6228e-03, 1.7002e-03, -3.6506e-03, -4.3168e-04]
      - [-2.4557e-02, -9.7102e-03, -5.8753e-03, 3.0456e-03, 5.8666e-03, 1.5037e-03, -2.0073e-03]
      - [-1.5400e-02, -5.2427e-03, -6.9148e-03, -5.9440e-03, -1.2183e-03, 2.1694e-03, 1.5925e-03]
    - collider: 'N2'
      eps: {A: 1.14813e+00, b: 4.60090e-02, Ea: -2.92413e+00}
    - collider: 'CO2'
      eps: {A: 8.98839e+01, b: -4.27974e-01, Ea: 2.41392e+02}
    - collider: 'H2O2'
      eps: {A: 6.45295e-01, b: 4.26266e-01, Ea: 4.28932e+01}
    - collider: 'H2O'
      eps: {A: 1.36377e+00, b: 3.06592e-01, Ea: 2.10079e+02}


.. _sec-yaml-interface-Arrhenius:

``interface-Arrhenius``
-----------------------

A reaction occurring on a surface between two bulk phases, or along an edge at the
intersection of two surfaces, as :ref:`described here <sec-surface-rate>`.

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
        as the leading-order activation energy term. This can be a scalar value for
        the linear dependency or a list of four values for the polynomial dependency
        given in the order of 1st, 2nd, 3rd, and 4th-order coefficients

    or a list containing the three elements above, in the given order.

    Note that parameters ``a``, ``m`` and ``E`` correspond to parameters
    :math:`\eta_{ki}`, :math:`\mu_{ki}` and :math:`\epsilon_{ki}` in Eq 11.113 of
    :cite:t:`kee2003`, respectively.

Examples::

    - equation: 2 H(s) => H2 + 2 Pt(s)
      rate-constant: {A: 3.7e21 cm^2/mol/s, b: 0, Ea: 67400 J/mol}
      coverage-dependencies: {H(s): {a: 0, m: 0, E: -6000 J/mol}}

    - equation: 2 O(S) => O2 + 2 Pt(S)
      rate-constant: {A: 3.7e+21, b: 0, Ea: 213200 J/mol}
      coverage-dependencies: {O(S): {a: 0.0, m: 0.0,
        E: [1.0e3 J/mol, 3.0e3 J/mol , -7.0e4 J/mol , 5.0e3 J/mol]}

    - equation: CH4 + PT(S) + O(S) => CH3(S) + OH(S)
      rate-constant: {A: 5.0e+18, b: 0.7, Ea: 4.2e+04}
      coverage-dependencies:
        O(S): [0, 0, 8000]
        PT(S): [0, -1.0, 0]

    - equation: 2 O(S) => O2 + 2 Pt(S)
      rate-constant: {A: 3.7e+21, b: 0, Ea: 213200 J/mol}
      coverage-dependencies:
        O(S): [0, 0, [1.0e6, 3.0e6, -7.0e7, 5.0e6]]

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
:ref:`described here <sec-sticking-rate>`.

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

Interface reactions involving :ref:`charge transfer <sec-electrochemical-reactions>`
between phases.

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
