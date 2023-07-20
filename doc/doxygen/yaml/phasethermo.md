# Phase Thermodynamic Models {#sec-yaml-phase-thermo-models}

[TOC]

Thermodynamic models are specified using the `thermo` field of @ref sec-yaml-phases.


## binary-solution-tabulated {#sec-yaml-binary-solution-tabulated}

A phase implementing tabulated standard state thermodynamics for one
species in a binary solution, see @ref Cantera.BinarySolutionTabulatedThermo.

Includes the fields of <tt>@ref sec-yaml-ideal-condensed</tt>, plus:

<b>`tabulated-species`</b>

The name of the species to which the tabulated enthalpy and entropy
is added.

<b>`tabulated-thermo`</b>

A mapping containing three (optionally four) lists of equal lengths:

-   `mole-fractions`: A list of mole fraction values for the tabulated species.
-   `enthalpy`: The extra molar enthalpy to be added to the tabulated species at
    these mole fractions.
-   `entropy`: The extra molar entropy to be added to the tabulated species at
    these mole fractions.
-   `molar-volume`: The molar volume of the phase at these mole fractions. This
    input is optional.

## compound-lattice {#sec-yaml-compound-lattice}

A phase that is comprised of a fixed additive combination of other
lattice phases, see @ref Cantera.LatticeSolidPhase.

Additional fields:

<b>`composition`</b>

A mapping of component phase names to their relative stoichiometries.

**Example:**

``` yaml
thermo: compound-lattice
composition: {Li7Si3(s): 1.0, Li7Si3-interstitial: 1.0}
```

## coverage-dependent-surface {#sec-yaml-coverage-dependent-surface}

A coverage-dependent surface phase. That is, a surface phase where the
enthalpy, entropy, and heat capacity of each species may depend on its
coverage and the coverage of other species in the phase. For full details,
see @ref Cantera.CoverageDependentSurfPhase.
The majority of coverage dependency parameters are provided in the species
entry, see @ref sec-yaml-coverage-dependent-surface.

Additional fields:

<b>`site-density`</b>

The molar density of surface sites.

<b>`reference-state-coverage`</b>

The reference state coverage denoting the low-coverage limit
(ideal-surface) thermodynamic properties.

**Example:**

``` yaml
- name: covdep
  thermo: coverage-dependent-surface
  species: [Pt, OC_Pt, CO2_Pt, C_Pt, O_Pt]
  state:
    T: 500.0
    P: 1.01325e+05
    coverages: {Pt: 0.5, OC_Pt: 0.5, CO2_Pt: 0.0, C_Pt: 0.0, O_Pt: 0.0}
  site-density: 2.72e-09
  reference-state-coverage: 0.22
```

## Debye-Huckel {#sec-yaml-Debye-Huckel}

The Debye-Hückel model, see @ref Cantera.DebyeHuckel.

Additional parameters for this model are contained in the `activity-data` field:

<b>`activity-data`</b>

The activity data field contains the following fields:

-   `model`: One of `dilute-limit`, `B-dot-with-variable-a`,
    `B-dot-with-common-a`, `beta_ij`, or `Pitzer-with-beta_ij`
-   `A_Debye`: The value of the Debye \"A\" parameter, or the string `variable`
    to use a calculation based on the water equation of state.
    Defaults to the constant value of 1.172576 kg^0.5/gmol^0.5, a
    nominal value for water at 298 K and 1 atm.
-   `B_Debye`: The Debye \"B\" parameter. Defaults to 3.2864e+09
    kg^0.5/gmol^0.5/m, a nominal value for water.
-   `max-ionic-strength`: The maximum ionic strength
-   `use-Helgeson-fixed-form`: Boolean, `true` or `false`
-   `default-ionic-radius`: Ionic radius to use for species where the ionic radius has
    not been specified.
-   `B-dot`: The value of B-dot.
-   `beta`: List of mappings providing values of \f$\beta_{ij}\f$ for different
    species pairs. Each mapping contains a `species` key that
    contains a list of two species names, and a `beta` key that
    contains the corresponding value of \f$\beta_{ij}\f$.

**Example:**

``` yaml
thermo: Debye-Huckel
activity-data:
  model: beta_ij
  max-ionic-strength: 3.0
  use-Helgeson-fixed-form: true
  default-ionic-radius: 3.042843 angstrom
  beta:
  - species: [H+, Cl-]
    beta: 0.27
  - species: [Na+, Cl-]
    beta: 0.15
  - species: [Na+, OH-]
    beta: 0.06
```

In addition, the Debye-Hückel model uses several species-specific
properties which may be defined in the `Debye-Huckel` field of the
*species* entry. These properties are:

<b>`ionic-radius`</b>

Size of the species.

<b>`electrolyte-species-type`</b>

One of `solvent`, `charged-species`, `weak-acid-associated`,
`strong-acid-associated`, `polar-neutral`, or `nonpolar-neutral`.
The type `solvent` is the default for the first species in the
phase. The type `charged-species` is the default for species with a
net charge. Otherwise, the default is and `nonpolar-neutral`.

<b>`weak-acid-charge`</b>

Charge to use for species that can break apart into charged species.

**Example:**

``` yaml
name: NaCl(aq)
composition: {Na: 1, Cl: 1}
thermo:
  model: piecewise-Gibbs
  h0: -96.03E3 cal/mol
  dimensionless: true
  data: {298.15: -174.5057463, 333.15: -174.5057463}
equation-of-state:
  model: constant-volume
  molar-volume: 1.3
Debye-Huckel:
  ionic-radius: 4 angstrom
  electrolyte-species-type: weak-acid-associated
  weak-acid-charge: -1.0
```

## edge {#sec-yaml-edge}

A one-dimensional edge between two surfaces, see @ref Cantera.EdgePhase

Additional fields:

<b>`site-density`</b>

The molar density of sites per unit length along the edge

**Example:**

``` yaml
thermo: edge
site-density: 5.0e-17 mol/cm
```

## electron-cloud {#sec-yaml-electron-cloud}

A phase representing an electron cloud, such as conduction electrons in
a metal, see @ref Cantera.MetalPhase.

Additional fields:

<b>`density`</b>

The density of the bulk metal

## fixed-stoichiometry {#sec-yaml-fixed-stoichiometry}

A phase with fixed composition, see @ref Cantera.StoichSubstance.

## HMW-electrolyte {#sec-yaml-HMW-electrolyte}

A dilute or concentrated liquid electrolyte phase that obeys the Pitzer
formulation for nonideality, see @ref Cantera.HMWSoln.

Additional parameters for this model are contained in the `activity-data` field:

<b>`activity-data`</b>

The activity data field contains the following fields:

-   `temperature-model`: The form of the Pitzer temperature model. One of `constant`,
    `linear` or `complex`. The default is `constant`.

-   `A_Debye`: The value of the Debye \"A\" parameter, or the string `variable`
    to use a calculation based on the water equation of state. The
    default is 1.172576 kg^0.5/gmol^0.5, a nominal value for water
    at 298 K and 1 atm.

-   `max-ionic-strength`: The maximum ionic strength
-   `interactions`: A list of mappings, where each mapping describes a binary or
    ternary interaction among species. Fields of this mapping
    include:

    -   `species`: A list of one to three species names
    -   `beta0`: The \f$\beta^{(0)}\f$ parameters for an cation/anion
        interaction. 1, 2, or 5 values depending on the value of
        `temperature-model`.
    -   `beta1`: The \f$\beta^{(1)}\f$ parameters for an cation/anion
        interaction. 1, 2, or 5 values depending on the value of
        `temperature-model`.
    -   `beta2`: The \f$\beta^{(2)}\f$ parameters for an cation/anion
        interaction. 1, 2, or 5 values depending on the value of
        `temperature-model`.
    -   `Cphi`: The \f$C^\phi\f$ parameters for an cation/anion interaction. 1,
        2, or 5 values depending on the value of
        `temperature-model`.
    -   `alpha1`: The \f$\alpha^{(1)}\f$ parameter for an cation/anion interaction.
    -   `alpha2`: The \f$\alpha^{(2)}\f$ parameter for an cation/anion interaction.
    -   `theta`: The \f$\theta\f$ parameters for a like-charged binary
        interaction. 1, 2, or 5 values depending on the value of
        `temperature-model`.
    -   `lambda`: The \f$\lambda\f$ parameters for binary interactions involving
        at least one neutral species. 1, 2, or 5 values depending on
        the value of `temperature-model`.
    -   `psi`: The \f$\Psi\f$ parameters for ternary interactions involving
        three charged species. 1, 2, or 5 values depending on the
        value of `temperature-model`.
    -   `zeta`: The \f$\zeta\f$ parameters for ternary interactions involving
        one neutral species. 1, 2, or 5 values depending on the
        value of `temperature-model`.
    -   `mu`: The \f$\mu\f$ parameters for a neutral species self-interaction.
        1, 2, or 5 values depending on the value of
        `temperature-model`.
-   `cropping-coefficients`:
    -   `ln_gamma_k_min`: Default -5.0.
    -   `ln_gamma_k_max`: Default 15.0.
    -   `ln_gamma_o_min`: Default -6.0.
    -   `ln_gamma_o_max`: Default 3.0.

**Example:**

``` yaml
thermo: HMW-electrolyte
activity-data:
  temperature-model: complex
  A_Debye: 1.175930 kg^0.5/gmol^0.5
  interactions:
  - species: [Na+, Cl-]
    beta0: [0.0765, 0.008946, -3.3158E-6, -777.03, -4.4706]
    beta1: [0.2664, 6.1608E-5, 1.0715E-6, 0.0, 0.0]
    beta2: [0.0, 0.0, 0.0, 0.0, 0.0]
    Cphi: [0.00127, -4.655E-5, 0.0, 33.317, 0.09421]
    alpha1: 2.0
  - species: [H+, Cl-]
    beta0: [0.1775]
    beta1: [0.2945]
    beta2: [0.0]
    Cphi: [0.0008]
    alpha1: 2.0
  - species: [Na+, OH-]
    beta0: 0.0864
    beta1: 0.253
    beta2: 0.0
    Cphi: 0.0044
    alpha1: 2.0
    alpha2: 0.0
  - {species: [Cl-, OH-], theta: -0.05}
  - {species: [Na+, Cl-, OH-], psi: -0.006}
  - {species: [Na+, H+], theta: 0.036}
  - {species: [Cl-, Na+, H+], psi: [-0.004]}
```

## ideal-condensed {#sec-yaml-ideal-condensed}

A condensed phase ideal solution, see @ref Cantera.IdealSolidSolnPhase.

Additional fields:

<b>`standard-concentration-basis`</b>

A string specifying the basis for the standard concentration. One of
`unity`, `species-molar-volume`, or `solvent-molar-volume`.

## ideal-gas {#sec-yaml-ideal-gas}

The ideal gas model, see @ref Cantera.IdealGasPhase.

**Example:**

``` yaml
- name: ohmech
  thermo: ideal-gas
  elements: [O, H, Ar, N]
  species: [H2, H, O, O2, OH, H2O, HO2, H2O2, AR, N2]
  kinetics: gas
  transport: mixture-averaged
  state: {T: 300.0, P: 1 atm}
```

## ideal-molal-solution {#sec-yaml-ideal-molal-solution}

A phase based on the mixing-rule assumption that all molality-based
activity coefficients are equal to one, see @ref Cantera.IdealMolalSoln.

Additional fields:

<b>`standard-concentration-basis`</b>

A string specifying the basis for the standard concentration. One of
`unity`, `species-molar-volume`, or `solvent-molar-volume`.

<b>`cutoff`</b>

Parameters for cutoff treatments of activity coefficients

-   `model`: `poly` or `polyExp`
-   `gamma_o`: gamma_o value for the cutoff process at the zero solvent point
-   `gamma_k`: gamma_k minimum for the cutoff process at the zero solvent point
-   `X_o`: value of the solute mole fraction that centers the cutoff
    polynomials for the cutoff = 1 process
-   `c_0`: Parameter in the polyExp cutoff treatment having to do with rate
    of exponential decay
-   `slope_f`: Parameter in the `polyExp` cutoff treatment
-   `slope_g`: Parameter in the `polyExp` cutoff treatment

**Example:**

``` yaml
thermo: ideal-molal-solution
standard-concentration-basis: solvent-molar-volume
cutoff:
  model: polyexp
  gamma_o: 0.0001
  gamma_k: 10.0
  X_o: 0.2
  c_0: 0.05
  slope_f: 0.6
  slope_g: 0.0
```

## ideal-solution-VPSS {#sec-yaml-ideal-solution-VPSS}

An ideal solution model using variable pressure standard state methods, see
@ref Cantera.IdealSolnGasVPSS.

Additional fields:

<b>`standard-concentration-basis`</b>

A string specifying the basis for the standard concentration. One of
`unity`, `species-molar-volume`, or `solvent-molar-volume`.

## ideal-surface {#sec-yaml-ideal-surface}

An ideal surface phase, see @ref Cantera.SurfPhase.

Additional fields:

<b>`site-density`</b>

The molar density of surface sites

**Example:**

``` yaml
- name: Pt_surf
  thermo: ideal-surface
  adjacent-phases: [gas]
  elements: [Pt, H, O, C]
  species: [PT(S), H(S), H2O(S), OH(S), CO(S), CO2(S), CH3(S), CH2(S)s,
    CH(S), C(S), O(S)]
  kinetics: surface
  reactions: all
  state:
    T: 900.0
    coverages: {O(S): 0.0, PT(S): 0.5, H(S): 0.5}
  site-density: 2.7063e-09
```

## ions-from-neutral-molecule {#sec-yaml-ions-from-neutral-molecule}

A model that handles the specification of the chemical potentials for
ionic species, given a specification of the chemical potentials for the
same phase expressed in terms of combinations of the ionic species that
represent neutral molecules, see @ref Cantera.IonsFromNeutralVPSSTP.

@deprecated This phase model is deprecated and will be removed after %Cantera 3.0.

Additional fields:

<b>`neutral-phase`</b>

The `name` of the phase definition for the phase containing the
neutral molecules.

**Example:**

``` yaml
- name: KCl-ions
  thermo: ions-from-neutral-molecule
  neutral-phase: KCl-neutral
  species: [K+, Cl-]
- name: KCl-neutral
  species: [KCl(l)]
  thermo: Margules
```

## lattice {#sec-yaml-lattice}

A simple thermodynamic model for a bulk phase, assuming a lattice of
solid atoms, see @ref Cantera.LatticePhase.

Additional fields:

<b>`site-density`</b>

The molar density of lattice sites

## liquid-water-IAPWS95 {#sec-yaml-liquid-water-IAPWS95}

An equation of state for liquid water, see @ref Cantera.WaterSSTP.

## Margules {#sec-yaml-Margules}

A phase employing the Margules approximation for the excess Gibbs free
energy, see @ref Cantera.MargulesVPSSTP.

Additional fields:

<b>`interactions`</b>

A list of mappings, where each mapping has the following fields:

-   `species`: A list of two species names
-   `excess-enthalpy`: A list of two values specifying the first and second excess
    enthalpy coefficients for the interaction of the specified
    species. Defaults to \[0, 0\].
-   `excess-entropy`: A list of two values specifying the first and second excess
    entropy coefficients for the interaction of the specified
    species. Defaults to \[0, 0\].
-   `excess-volume-enthalpy`: A list of two values specifying the first and second
    enthalpy coefficients for the excess volume interaction of the specified
    species. Defaults to \[0, 0\].
-   `excess-volume-entropy`: A list of two values specifying the first and second
    entropy coefficients for the excess volume interaction of the specified
    species. Defaults to \[0, 0\].

**Example:**

``` yaml
thermo: Margules
interactions:
- species: [KCl(l), LiCl(l)]
  excess-enthalpy: [-17570, -377]
  excess-entropy: [-7.627, 4.958]
```

## Maskell-solid-solution {#sec-yaml-Maskell-solid-solution}

A condensed phase non-ideal solution with two species, see
@ref Cantera.MaskellSolidSolnPhase.

@deprecated This phase model is deprecated and will be removed after %Cantera 3.0.

Additional fields:

<b>`excess-enthalpy`</b>

The molar excess enthalpy

<b>`product-species`</b>

String specifying the \"product\" species

**Example:**

``` yaml
thermo: Maskell-solid-solution
excess-enthalpy: 5 J/mol
product-species: H(s)
```

## Peng-Robinson {#sec-yaml-Peng-Robinson}

A multi-species Peng-Robinson phase, see @ref Cantera.PengRobinson.

The parameters for each species are contained in the corresponding
species entries. See @ref sec-yaml-eos-peng-robinson.

## plasma {#sec-yaml-plasma}

A phase for plasma. This phase handles plasma properties such as the
electron energy distribution and electron temperature with different
models, see @ref Cantera.PlasmaPhase.

Additional fields:

<b>`electron-energy-distribution`</b>

A mapping with the following fields:

-   `type`: String specifying the type of the electron energy distribution
    to be used. Supported model strings are:
    -   `isotropic`
    -   `discretized`
-   `shape-factor`: A constant in the isotropic distribution, which is shown as x in
    the detailed description of this class. The value needs to be a
    positive number. This field is only used with `isotropic`. Defaults to 2.0.
-   `mean-electron-energy`: Mean electron energy of the isotropic distribution. The
    default sets the electron temperature equal gas temperature and uses the
    corresponding electron energy as mean electron energy. This
    field is only used with `isotropic`.
-   `energy-levels`: A list of values specifying the electron energy levels. The
    default uses 1001 equal spaced points from 0 to 1 eV.
-   `distribution`: A list of values specifying the discretized electron energy
    distribution. This field is only used with `discretized`.
-   `normalize`: A flag specifying whether normalizing the discretized electron energy
    distribution or not. This field is only used with `discretized`. Defaults to `true`.

**Example:**

``` yaml
- name: isotropic-electron-energy-plasma
  thermo: plasma
  kinetics: gas
  transport: ionized-gas
  electron-energy-distribution:
    type: isotropic
    shape-factor: 2.0
    mean-electron-energy: 1.0 eV
    energy-levels: [0.0, 0.1, 1.0, 10.0]
- name: discretized-electron-energy-plasma
  thermo: plasma
  kinetics: gas
  transport: ionized-gas
  electron-energy-distribution:
    type: discretized
    energy-levels: [0.0, 0.1, 1.0, 10.0]
    distribution: [0.0, 0.2, 0.7, 0.01]
    normalize: False
```

## pure-fluid {#sec-yaml-pure-fluid}

A phase representing a pure fluid equation of state for one of several
species, see @ref Cantera.PureFluidPhase.

Additional fields:

<b>`pure-fluid-name`</b>

Name of the pure fluid model to use: - `carbon-dioxide` -
`heptane` - `HFC-134a` - `hydrogen` - `methane` - `nitrogen` -
`oxygen` - `water`

## Redlich-Kister {#sec-yaml-Redlich-Kister}

A phase employing the Redlich-Kister approximation for the excess Gibbs
free energy, see @ref Cantera.RedlichKisterVPSSTP.

Additional fields:

<b>`interactions`</b>

A list of mappings, where each mapping has the following fields:

-   `species`: A list of two species names
-   `excess-enthalpy`: A list of polynomial coefficients for the excess enthalpy of the
    specified binary interaction
-   `excess-entropy`: A list of polynomial coefficients for the excess entropy of the
    specified binary interaction

**Example:**

``` yaml
thermo: Redlich-Kister
interactions:
- species: [Li(C6), V(C6)]
  excess-enthalpy: [-3.268e+06, 3.955e+06, -4.573e+06, 6.147e+06, -3.339e+06,
                    1.117e+07, 2.997e+05, -4.866e+07, 1.362e+05, 1.373e+08,
                    -2.129e+07, -1.722e+08, 3.956e+07, 9.302e+07, -3.280e+07]
  excess-entropy: [0.0]
```

## Redlich-Kwong {#sec-yaml-Redlich-Kwong}

A multi-species Redlich-Kwong phase, see @ref Cantera.RedlichKwongMFTP.

The parameters for each species are contained in the corresponding
species entries. See @ref sec-yaml-eos-redlich-kwong.
