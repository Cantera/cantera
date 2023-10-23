.. highlight:: yaml

.. _sec-yaml-species:

*******
Species
*******

The fields of a ``species`` entry are:

``name``
    String identifier used for the species. Required.

``composition``
    Mapping that specifies the elemental composition of the species,
    for example, ``{C: 1, H: 4}``. Required.

``thermo``
    Mapping containing the reference state thermodynamic model specification
    and parameters. See :ref:`sec-yaml-species-thermo`.

``equation-of-state``
    A mapping or list of mappings. Each mapping contains an equation of state
    model specification for the species, any parameters for that model, and any
    parameters for interactions with other species. See
    :ref:`sec-yaml-species-eos`.

``critical-parameters``
    Mapping containing parameters related to the critical state of a species. Used in
    models that incorporate "real gas" effects, such as
    :ref:`sec-yaml-eos-redlich-kwong`.
    See :ref:`sec-yaml-species-crit-props`.

``transport``
    Mapping containing the species transport model specification and
    parameters. See :ref:`sec-yaml-species-transport`.

``coverage-dependencies``
    Mapping where keys are the names of other species whose coverages affects
    thermodynamic properties of this species. The map values are the dependency entries
    including ``model`` and model-specific parameters as :ref:`described here
    <sec-yaml-species-coverage>`.

``sites``
    The number of sites occupied by a surface or edge species. Default is 1.

``Debye-Huckel``
    Additional model parameters used in the Debye-Hückel model. See
    :ref:`sec-yaml-Debye-Huckel` for more information.


.. _sec-yaml-species-thermo:

Species thermo models
=====================

Fields of a species ``thermo`` entry used by all models are:

``model``
    String specifying the model to be used. Required. Supported model strings
    are:

    - ``NASA7`` (:ref:`details <sec-yaml-nasa7>`)
    - ``NASA9`` (:ref:`details <sec-yaml-nasa9>`)
    - ``Shomate`` (:ref:`details <sec-yaml-shomate>`)
    - ``constant-cp`` (:ref:`details <sec-yaml-constcp>`)
    - ``piecewise-Gibbs`` (:ref:`details <sec-yaml-piecewise-gibbs>`)

``reference-pressure``
    The reference pressure at which the given thermodynamic properties apply.
    Defaults to 1 atm.


.. _sec-yaml-nasa7:

NASA 7-coefficient polynomials
------------------------------

The polynomial form :ref:`described here <sec-thermo-nasa7>`,
given for one or two temperature regions. Additional fields of a ``NASA7``
thermo entry are:

``temperature-ranges``
    A list giving the temperature intervals on which the polynomials are valid.
    For one temperature region, this list contains the minimum and maximum
    temperatures for the polynomial. For two temperature regions, this list
    contains the minimum, intermediate, and maximum temperatures.

``data``
    A list with one item per temperature region, where that item is a 7 item
    list of polynomial coefficients. The temperature regions are arranged in
    ascending order. Note that this is different from the standard CHEMKIN
    formulation that uses two temperature regions listed in descending order.

Example::

    thermo:
      model: NASA7
      temperature-ranges: [300.0, 1000.0, 5000.0]
      data:
      - [3.298677, 0.0014082404, -3.963222e-06, 5.641515e-09,
        -2.444854e-12, -1020.8999, 3.950372]
      - [2.92664, 0.0014879768, -5.68476e-07, 1.0097038e-10,
        -6.753351e-15, -922.7977, 5.980528]


.. _sec-yaml-nasa9:

NASA 9-coefficient polynomials
------------------------------

The polynomial form :ref:`described here <sec-thermo-nasa9>`,
given for any number of temperature regions. Additional fields of a ``NASA9``
thermo entry are:

``temperature-ranges``
    A list giving the temperature intervals on which the polynomials are valid.
    This list contains the minimum temperature, the intermediate temperatures
    between each set pair of regions, and the maximum temperature.

``data``
    A list with one item per temperature region, where that item is a 9 item
    list of polynomial coefficients. The temperature regions are arranged in
    ascending order.

Example::

    thermo:
      model: NASA9
      temperature-ranges: [200.00, 1000.00, 6000.0, 20000]
      reference-pressure: 1 bar
      data:
      - [2.210371497E+04, -3.818461820E+02, 6.082738360E+00, -8.530914410E-03,
         1.384646189E-05, -9.625793620E-09, 2.519705809E-12, 7.108460860E+02,
         -1.076003744E+01]
      - [5.877124060E+05, -2.239249073E+03, 6.066949220E+00, -6.139685500E-04,
         1.491806679E-07,  -1.923105485E-11, 1.061954386E-15, 1.283210415E+04,
         -1.586640027E+01]
      - [8.310139160E+08, -6.420733540E+05, 2.020264635E+02, -3.065092046E-02,
         2.486903333E-06, -9.705954110E-11, 1.437538881E-15, 4.938707040E+06,
         -1.672099740E+03]


.. _sec-yaml-shomate:

Shomate polynomials
-------------------

The polynomial form :ref:`described here <sec-thermo-shomate>`,
given for one or two temperature regions. Additional fields of a ``Shomate``
thermo entry are:

``temperature-ranges``
    A list giving the temperature intervals on which the polynomials are valid.
    For one temperature region, this list contains the minimum and maximum
    temperatures for the polynomial. For two temperature regions, this list
    contains the minimum, intermediate, and maximum temperatures.

``data``
    A list with one item per temperature region, where that item is a 7 item
    list of polynomial coefficients. The temperature regions are arranged in
    ascending order.

Example::

    thermo:
      model: Shomate
      temperature-ranges: [298, 1300, 6000]
      reference-pressure: 1 bar
      data:
      - [25.56759, 6.096130, 4.054656, -2.671301, 0.131021,
        -118.0089, 227.3665]
      - [35.15070, 1.300095, -0.205921, 0.013550, -3.282780,
        -127.8375, 231.7120]


.. _sec-yaml-constcp:

Constant heat capacity
----------------------

The constant heat capacity model :ref:`described here <sec-thermo-const-cp>`.
Additional fields of a ``constant-cp`` thermo entry are:

``T0``
    The reference temperature. Defaults to 298.15 K.

``h0``
    The molar enthalpy at the reference temperature. Defaults to 0.0.

``s0``
    The molar entropy at the reference temperature. Defaults to 0.0.

``cp0``
    The heat capacity at constant pressure. Defaults to 0.0.

``T-min``
    The minimum temperature at which this thermo data should be used.
    Defaults to 0.0.

``T-max``
    The maximum temperature at which this thermo data should be used.
    Defaults to infinity.

Example::

    thermo:
      model: constant-cp
      T0: 1000 K
      h0: 9.22 kcal/mol
      s0: -3.02 cal/mol/K
      cp0: 5.95 cal/mol/K


.. _sec-yaml-piecewise-gibbs:

Piecewise Gibbs
---------------

A model based on piecewise interpolation of the Gibbs free energy as
:ct:`described here <Mu0Poly>`. Additional fields of a ``piecewise-Gibbs`` entry are:

``h0``
    The molar enthalpy at the reference temperature of 298.15 K. Defaults to
    0.0.

``dimensionless``
    A boolean flag indicating whether the values of the Gibbs free energy are
    given in a dimensionless form, that is, divided by :math:`RT`. Defaults to
    ``false``.

``data``
    A mapping of temperatures to values of the Gibbs free energy. The Gibbs free
    energy can be either in molar units (if ``dimensionless`` is ``false``) or
    nondimensionalized by the corresponding temperature (if ``dimensionless`` is
    ``true``). A value must be provided at :math:`T^\circ = 298.15` K.

``T-min``
    The minimum temperature at which this thermo data should be used.
    Defaults to 0.0.

``T-max``
    The maximum temperature at which this thermo data should be used.
    Defaults to infinity.

Example::

    thermo:
      model: piecewise-Gibbs
      h0: -230.015 kJ/mol
      dimensionless: true
      data: {298.15: -91.50963, 333.15: -85.0}


.. _sec-yaml-species-crit-props:

Species critical state parameters
=================================

``critical-temperature``
    The critical temperature of the species [K]

``critical-pressure``
    The critical pressure of the species [Pa]

``acentric-factor``
    Pitzer's acentric factor :math:`omega` [-]


.. _sec-yaml-species-eos:

Species equation of state models
================================

``model``
    String specifying the model to be used. Required. Supported model strings
    are:

    - ``constant-volume`` (:ref:`details <sec-yaml-eos-constant-volume>`)
    - ``density-temperature-polynomial`` (:ref:`details <sec-yaml-eos-density-temperature-polynomial>`)
    - ``HKFT`` (:ref:`details <sec-yaml-eos-hkft>`)
    - ``liquid-water-IAPWS95`` (:ref:`details <sec-yaml-eos-liquid-water-iapws95>`)
    - ``molar-volume-temperature-polynomial`` (:ref:`details <sec-yaml-eos-molar-volume-temperature-polynomial>`)
    - ``Peng-Robinson`` (:ref:`details <sec-yaml-eos-peng-robinson>`)
    - ``Redlich-Kwong`` (:ref:`details <sec-yaml-eos-redlich-kwong>`)


.. _sec-yaml-eos-constant-volume:

Constant volume
---------------

A constant volume model as :ct:`described here <PDSS_ConstVol>`.

Any one of the following may be specified:

``molar-volume``
    The molar volume of the species.

``molar-density``
    The molar density of the species.

``density``
    The mass density of the species.

Example::

    equation-of-state:
      model: constant-volume
      molar-volume: 1.3 cm^3/mol


.. _sec-yaml-eos-density-temperature-polynomial:

Density temperature polynomial
------------------------------

A model in which the density varies with temperature as
:ct:`described here <PDSS_SSVol>`.

Additional fields:

``data``
    Vector of 4 coefficients for a cubic polynomial in temperature

Example::

    equation-of-state:
      model: density-temperature-polynomial
      units: {mass: g, length: cm}
      data: [0.536504, -1.04279e-4, 3.84825e-9, -5.2853e-12]


.. _sec-yaml-eos-hkft:

HKFT
----

The Helgeson-Kirkham-Flowers-Tanger model for aqueous species as
:ct:`described here <PDSS_HKFT>`.

Additional fields:

``h0``
    Enthalpy of formation at the reference temperature and pressure

``s0``
    Entropy of formation at the reference temperature and pressure

``a``
    4-element vector containing the coefficients :math:`a_1, \ldots , a_4`

``c``
    2-element vector containing the coefficients :math:`c_1` and :math:`c_2`

``omega``
    The :math:`\omega` parameter at the reference temperature and pressure

Example::

    equation-of-state:
      model: HKFT
      h0: -57433. cal/gmol
      s0: 13.96 cal/gmol/K
      a: [0.1839 cal/gmol/bar, -228.5 cal/gmol,
         3.256 cal*K/gmol/bar, -27260. cal*K/gmol]
      c: [18.18 cal/gmol/K, -29810. cal*K/gmol]
      omega: 33060 cal/gmol


.. _sec-yaml-eos-liquid-water-iapws95:

Liquid Water IAPWS95
--------------------

A detailed equation of state for liquid water as :ct:`described here <PDSS_Water>`.


.. _sec-yaml-eos-molar-volume-temperature-polynomial:

Molar volume temperature polynomial
-----------------------------------

A model in which the molar volume varies with temperature as
:ct:`described here <PDSS_SSVol>`.

Additional fields:

``data``
    Vector of 4 coefficients for a cubic polynomial in temperature


.. _sec-yaml-eos-peng-robinson:

Peng-Robinson
-------------

A model where species follow the Peng-Robinson real gas equation of state as
:ct:`described here <PengRobinson>`.

Additional fields:

``a``
    Pure-species ``a`` coefficient [Pa*m^6/kmol^2]

``b``
    Pure-species ``b`` coefficient [m^3/kmol]

``acentric-factor``
    Pitzer's acentric factor [-]

``binary-a``
    Optional mapping where the keys are species names and the values are the ``a``
    coefficients for binary interactions between the two species.

Example::

    equation-of-state:
      model: Peng-Robinson
      units: {length: cm, quantity: mol}
      a: 5.998873E+11
      b: 18.9714
      acentric-factor: 0.344
      binary-a:
        H2: 4 bar*cm^6/mol^2
        CO2: 7.897e7 bar*cm^6/mol^2


.. _sec-yaml-eos-redlich-kwong:

Redlich-Kwong
-------------

A model where species follow the Redlich-Kwong real gas equation of state as
:ct:`described here <RedlichKwongMFTP>`.

Additional fields:

``a``
    Pure-species ``a`` coefficient. Scalar or list of two values for a
    temperature-dependent expression.

``b``
    Pure-species ``b`` coefficient.

``binary-a``
    Mapping where the keys are species and the values are the ``a``
    coefficients for binary interactions between the two species.


.. _sec-yaml-species-transport:

Species transport models
========================

``model``
    String specifying the model type. The only model that is specifically
    handled is ``gas``.

Gas transport
-------------

Species transport properties are a rare exception to Cantera's use of SI units,
and use the units in which these properties are customarily reported. No
conversions are supported.

The additional fields of a ``gas`` transport entry are:

``geometry``
    A string specifying the geometry of the molecule. One of ``atom``,
    ``linear``, or ``nonlinear``.

``diameter``
    The Lennard-Jones collision diameter [Å]

``well-depth``
    The Lennard-Jones well depth [K]

``dipole``
    The permanent dipole moment [Debye]. Default 0.0.

``polarizability``
    The dipole polarizability [Å^3]. Default 0.0.

``rotational-relaxation``
    The rotational relaxation collision number at 298 K [-]. Default 0.0.

``acentric-factor``
    Pitzer's acentric factor [-]. Default 0.0. This value may also be specified as part
    of the :ref:`critical-parameters <sec-yaml-species-crit-props>` field, in which case
    the value provided there supersedes this one.

``dispersion-coefficient``
    The dispersion coefficient, normalized by :math:`e^2` [Å^5]. Default 0.0.

``quadrupole-polarizability``
    The quadrupole polarizability [Å^5]. Default 0.0.

Example::

    transport:
      model: gas
      geometry: linear
      well-depth: 107.4
      diameter: 3.458
      polarizability: 1.6
      rotational-relaxation: 3.8


.. _sec-yaml-species-coverage:

Species coverage dependencies
=============================

The coverage-dependent surface species formulation calculates coverage-dependent
correction factors to the ideal surface phase properties. Used in conjunction with the
:ref:`coverage-dependent-surface <sec-yaml-coverage-dependent-surface>` phase model.
Full details are :ct:`described here <CoverageDependentSurfPhase>`.

Fields of a species ``coverage-dependencies`` map entry used by all models are:

``model``
    String specifying the model to be used. Required. Supported model strings are:

    - ``linear`` (:ref:`details <sec-yaml-species-coverage-linear>`)
    - ``polynomial`` (:ref:`details <sec-yaml-species-coverage-polynomial>`)
    - ``piecewise-linear`` (:ref:`details <sec-yaml-species-coverage-piecewise-linear>`)
    - ``interpolative`` (:ref:`details <sec-yaml-species-coverage-interpolative>`)


.. _sec-yaml-species-coverage-linear:

Linear dependency model
-----------------------

``enthalpy``
    Slope of the coverage-dependent enthalpy.

``entropy``
    Slope of the coverage-dependent entropy.

Example::

    coverage-dependencies:
      O_Pt:
        model: linear
        units: {energy: eV, quantity: molec}
        enthalpy: 0.48
        entropy: -0.031
      # + other entries (optional)


.. _sec-yaml-species-coverage-polynomial:

Polynomial dependency model
---------------------------

``enthalpy-coefficients``
    Array of polynomial coefficients in order of 1st, 2nd, 3rd, and
    4th-order used in coverage-dependent enthalpy calculation.

``entropy-coefficients``
    Array of polynomial coefficients in order of 1st, 2nd, 3rd, and
    4th-order used in coverage-dependent entropy calculation.

Example::

    coverage-dependencies:
      OC_Pt:
        model: polynomial
        units: {energy: J, quantity: mol}
        enthalpy-coefficients: [0.0, -3.86e4, 0.0, 4.2e5]
        entropy-coefficients: [0.8e3, 0.0, -1.26e4, 0.0]
      # + other entries (optional)


.. _sec-yaml-species-coverage-piecewise-linear:

Piecewise-linear dependency model
---------------------------------

``enthalpy-low``
    Slope of the coverage-dependent enthalpy for the lower coverage region.

``entropy-low``
    Slope of the coverage-dependent entropy for the lower coverage region.

``enthalpy-high``
    Slope of the coverage-dependent enthalpy for the higher coverage region.

``entropy-high``
    Slope of the coverage-dependent entropy for the higher coverage region.

``enthalpy-change``
    Coverage that separates the lower and higher coverage regions of the
    coverage-dependent enthalpy.

``entropy-change``
    Coverage that separates the lower and higher coverage regions of the
    coverage-dependent entropy.

``heat-capacity-a``
    Coefficient :math:`c^{(a)}` used in the
    :ct:`coverage-dependent heat capacity <HeatCapacityDependency>` model.

``heat-capacity-b``
    Coefficient :math:`c^{(b)}` used in the
    :ct:`coverage-dependent heat capacity <HeatCapacityDependency>` model.

Example::

    coverage-dependencies:
      CO2_Pt:
        model: piecewise-linear
        units: {energy: kJ, quantity: mol}
        enthalpy-low: 0.5e2
        enthalpy-high: 1.0e2
        enthalpy-change: 0.4
        entropy-low: 0.1e2
        entropy-high: -0.2e2
        entropy-change: 0.4
        heat-capacity-a: 0.02e-1
        heat-capacity-b: -0.156e-1
      # + other entries (optional)


.. _sec-yaml-species-coverage-interpolative:

Interpolative dependency model
------------------------------

``enthalpy-coverages``
    Array of discrete coverage values used in coverage-dependent enthalpy.

``entropy-coverages``
    Array of discrete coverage values used in coverage-dependent entropy.

``enthalpies``
    Array of discrete enthalpy values corresponding to the
    coverages in ``enthalpy-coverages``.

``entropies``
    Array of discrete entropy values corresponding to the
    coverages in ``entropy-coverages``.

Example::

    coverage-dependencies:
      C_Pt:
        model: interpolative
        units: {energy: kcal, quantity: mol}
        enthalpy-coverages: [0.0, 0.2, 0.4, 0.7, 0.9, 1.0]
        entropy-coverages: [0.0, 0.5, 1.0]
        enthalpies: [0.0, 0.5, 1.0, 2.7, 3.5, 4.0]
        entropies: [0.0, -0.7, -2.0]
      # + other entries (optional)
