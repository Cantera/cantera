# Phase Thermodynamic Models

On this page, we list the phase thermodynamic models implemented in Cantera, with
links to the documentation for their YAML input parameters and the documentation for
the C++ classes which implement these models. This API documentation may also provide
references or a mathematical description of the model.


## Models for Gaseous Mixtures

(sec-ideal-gas-phase)=
Ideal Gas Mixture
: A mixture which follows the ideal gas law. Defined in the YAML format by specifying
  [`ideal-gas`](sec-yaml-ideal-gas) in the `thermo` field of the phase definition.
  Implemented by class {ct}`IdealGasPhase`.

(sec-Redlich-Kwong-phase)=
Redlich-Kwong Real Gas Mixture
: A multi-species real gas following the Redlich-Kwong equation of state. Defined in the
  YAML format by specifying [`Redlich-Kwong`](sec-yaml-Redlich-Kwong) in the `thermo`
  field of the phase definition. Implemented by class {ct}`RedlichKwongMFTP`.

(sec-Peng-Robinson-phase)=
Peng-Robinson Real Gas Mixture
: A multi-species real gas following the Peng-Robinson equation of state. Defined in the
  YAML format by specifying [`Peng-Robinson`](sec-yaml-Peng-Robinson) in the `thermo`
  field of the phase definition. Implemented by class {ct}`PengRobinson`.

(sec-plasma-phase)=
Plasma
: A phase that extends the ideal gas model to handle plasma properties such as the
  electron energy distribution and electron temperature with different models. Defined
  in the YAML format by specifying [`plasma`](sec-yaml-plasma) in the `thermo` field of
  the phase definition. Implemented by class {ct}`PlasmaPhase`.


## Models for Surfaces and Interfaces

(sec-ideal-surface-phase)=
Ideal Surface
: An interface between two bulk phases where the species behave as an ideal solution and
  the composition is described by the coverage of each species on the surface. Defined
  in the YAML format by specifying [`ideal-surface`](sec-yaml-ideal-surface) in the
  `thermo` field of the phase definition. Implemented by class {ct}`SurfPhase`.

(sec-coverage-dependent-surface-phase)=
Surface Phase with Coverage-Dependent Thermo
: A coverage-dependent surface phase. That is, a surface phase where the enthalpy,
  entropy, and heat capacity of each species may depend on its coverage and the coverage
  of other species in the phase. Defined in the YAML format by specifying
  [`coverage-dependent-surface`](sec-yaml-coverage-dependent-surface) in the `thermo`
  field of the phase definition. Implemented by class {ct}`CoverageDependentSurfPhase`.

(sec-edge-phase)=
Edge
: A one-dimensional edge between two surfaces, which defines a triple phase boundary.
  Defined in the YAML format by specifying [`edge`](sec-yaml-edge) in the `thermo` field
  of the phase definition. Implemented by class {ct}`EdgePhase`.


## Single-species Phase Models

(sec-fixed-stoichiometry-phase)=
Stoichiometric Substance
: A *stoichiometric substance* is one that is modeled as having a precise, fixed
  composition, given by the composition of the one species present. Defined in the YAML
  format by specifying [`fixed-stoichiometry`](sec-yaml-fixed-stoichiometry) in the
  `thermo` field of the phase definition. Implemented by class {ct}`StoichSubstance`.

(sec-electron-cloud-phase)=
Electron Cloud
: A phase representing an electron cloud, such as conduction electrons in a metal.
  Defined in the YAML format by specifying [`electron-cloud`](sec-yaml-electron-cloud)
  in the `thermo` field of the phase definition. Implemented by class {ct}`MetalPhase`.

(sec-pure-fluid-phase)=
Pure Fluid Phases
: A phase representing a pure fluid equation of state for one of several pure substances
  including liquid, vapor, two-phase, and supercritical regions. Defined in the YAML
  format by specifying [`pure-fluid`](sec-yaml-pure-fluid) in the `thermo` field of the
  phase definition. Implemented by class {ct}`PureFluidPhase`.

(sec-liquid-water-IAPWS95-phase)=
Liquid Water using the IAPWS95 Equation of State
: An implementation of the IAPWS95 equation of state for water {cite:p}`wagner2002`, for
  the liquid region only. Defined in the YAML format by specifying
  [`liquid-water-IAPWS95`](sec-yaml-liquid-water-IAPWS95) in the `thermo` field of the
  phase definition. Implemented by class {ct}`WaterSSTP`.


## Ideal Solid and Liquid Solutions

(sec-ideal-molal-solution-phase)=
Ideal Molal Solution
: An ideal solution based on the mixing-rule assumption that all molality-based activity
  coefficients are equal to one. Defined in the YAML format by specifying
  [`ideal-molal-solution`](sec-yaml-ideal-molal-solution) in the `thermo` field of the
  phase definition. Implemented by class {ct}`IdealMolalSoln`.

(sec-ideal-condensed-phase)=
Ideal Condensed Phase
: An ideal liquid or solid solution based on the mixing-rule assumption that all molar
  concentration-based activity coefficients are equal to one. Defined in the YAML format
  by specifying [`ideal-condensed`](sec-yaml-ideal-condensed) in the `thermo` field of
  the phase definition. Implemented by class {ct}`IdealSolidSolnPhase`.

(sec-ideal-solution-VPSS-phase)=
Ideal Condensed Phase with VPSS Species
: An ideal solution model using variable pressure standard state methods. This allows
  the standard state molar volume of species to be specified as a function of
  temperature. Defined in the YAML format by specifying
  [`ideal-solution-VPSS`](sec-yaml-ideal-solution-VPSS) in the `thermo` field of the
  phase definition. Implemented by class {ct}`IdealSolnGasVPSS`.

(sec-lattice-phase)=
Lattice Phase
: A simple thermodynamic model for a bulk phase, assuming an incompressible lattice of
  solid atoms. Defined in the YAML format by specifying [`lattice`](sec-yaml-lattice) in
  the `thermo` field of the phase definition. Implemented by class {ct}`LatticePhase`.

(sec-compound-lattice-phase)=
Compound Lattice Phase
: A phase that is comprised of a fixed additive combination of other lattice phases.
  Defined in the YAML format by specifying [`compound-lattice`](sec-yaml-compound-lattice)
  in the `thermo` field of the phase definition. Implemented by class
  {ct}`LatticeSolidPhase`.


## Non-ideal Solid and Liquid Solutions

(sec-binary-solution-tabulated-phase)=
Binary Solution with Tabulated Enthalpy and Entropy
: A phase representing a non-ideal binary solution where the excess enthalpy and entropy
  are interpolated between tabulated values as a function of mole fraction. Defined in
  the YAML format by specifying
  [`binary-solution-tabulated`](sec-yaml-binary-solution-tabulated) in the `thermo`
  field of the phase definition. Implemented by class
  {ct}`BinarySolutionTabulatedThermo`.

(sec-Debye-Huckel-phase)=
Debye-Huckel Solution
: A dilute liquid electrolyte which obeys the Debye-Hückel formulation for nonideality.
  Defined in the YAML format by specifying [`Debye-Huckel`](sec-yaml-Debye-Huckel) in
  the `thermo` field of the phase definition. Implemented by class {ct}`DebyeHuckel`.

(sec-HMW-electrolyte-phase)=
Harvie--Møller--Weare electrolyte
: A dilute or concentrated liquid electrolyte phase that obeys the Pitzer formulation
  for nonideality. Defined in the YAML format by specifying
  [`HMW-electrolyte`](sec-yaml-HMW-electrolyte) in the `thermo` field of the phase
  definition. Implemented by class {ct}`HMWSoln`.

(sec-Margules-phase)=
Margules Solution
: A condensed phase employing the Margules approximation for the excess Gibbs free
  energy. Defined in the YAML format by specifying [`Margules`](sec-yaml-Margules) in
  the `thermo` field of the phase definition. Implemented by class {ct}`MargulesVPSSTP`.

(sec-Redlich-Kister-phase)=
Redlich-Kister Solution
: A phase employing the Redlich-Kister approximation for the excess Gibbs free energy.
  Defined in the YAML format by specifying [`Redlich-Kister`](sec-yaml-Redlich-Kister)
  in the `thermo` field of the phase definition. Implemented by class
  {ct}`RedlichKisterVPSSTP`.
