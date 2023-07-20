# Phase Definitions {#sec-yaml-phases}

[TOC]

## Overview

For each phase that appears in a problem, a corresponding entry should
be present in the input file(s). The phase entry specifies the elements
and species present in that phase, and the models to be used for
computing thermodynamic, kinetic, and transport properties.

### Naming the Phase

The `name` entry is a string that identifies the phase. It must be
unique within the file among all phase definitions of any type. Phases
are referenced by name when importing them. The `name` is also used to
identify the phase within multiphase mixtures or at phase boundaries.

### Setting the Thermodynamic Model {#sec-yaml-phase-setting-thermo}

The thermodynamic model used to represent a phase is specified in the
`thermo` field. Supported models are:

-   <tt>@ref sec-yaml-binary-solution-tabulated</tt>:
    A binary mixture where the excess enthalpy and entropy are interpolated between
    tabulated values as a function of mole fraction
-   <tt>@ref sec-yaml-compound-lattice</tt>:
    A phase that is comprised of a fixed additive combination of other lattice phases
-   <tt>@ref sec-yaml-coverage-dependent-surface</tt>:
    A coverage-dependent surface phase
-   <tt>@ref sec-yaml-Debye-Huckel</tt>:
    A dilute liquid electrolyte which obeys the Debye-Hückel formulation for nonideality
-   <tt>@ref sec-yaml-edge</tt>:
    A one-dimensional edge between two surfaces
-   <tt>@ref sec-yaml-electron-cloud</tt>:
    A phase representing free electrons in a metal
-   <tt>@ref sec-yaml-fixed-stoichiometry</tt>:
    An incompressible, single-species phase
-   <tt>@ref sec-yaml-HMW-electrolyte</tt>:
    A dilute or concentrated liquid electrolyte which obeys the Pitzer formulation
    for nonideality
-   <tt>@ref sec-yaml-ideal-condensed</tt>:
    A condensed phase ideal solution
-   <tt>@ref sec-yaml-ideal-gas</tt>:
    A mixture which obeys the ideal gas law
-   <tt>@ref sec-yaml-ideal-molal-solution</tt>:
    An ideal solution based on the mixing-rule assumption that all molality-based
    activity coefficients are equal to one
-   <tt>@ref sec-yaml-ideal-solution-VPSS</tt>:
    An ideal solution; Uses \"variable pressure standard state\" methods for
    calculating thermodynamic properties
-   <tt>@ref sec-yaml-ideal-surface</tt>:
    A surface between two bulk phases
-   <tt>@ref sec-yaml-ions-from-neutral-molecule</tt>:
    A phase for representing ionic species based on another phase where those ions
    are components of neutral molecules
-   <tt>@ref sec-yaml-lattice</tt>:
    A simple model for an incompressible lattice of solid atoms
-   <tt>@ref sec-yaml-liquid-water-IAPWS95</tt>:
    An implementation of the IAPWS95 equation of state for water, for the
    liquid region only
-   <tt>@ref sec-yaml-Margules</tt>:
    A model that employs the Margules approximation for the excess Gibbs free energy
-   <tt>@ref sec-yaml-Maskell-solid-solution</tt>:
    A condensed, binary, non-ideal solution
-   <tt>@ref sec-yaml-Peng-Robinson</tt>:
    A multi-species mixture obeying the Peng-Robinson equation of state
-   <tt>@ref sec-yaml-plasma</tt>:
    A plasma phase handling plasma properties such as the electron energy distribution
    and electron temperature
-   <tt>@ref sec-yaml-pure-fluid</tt>:
    A phase representing one of several pure substances including liquid, vapor,
    two-phase, and supercritical regions
-   <tt>@ref sec-yaml-Redlich-Kister</tt>:
    A model that employs the Redlich-Kister approximation for the excess Gibbs
    free energy
-   <tt>@ref sec-yaml-Redlich-Kwong</tt>:
    A multi-species mixture obeying the Redlich-Kwong equation of state.

Some thermodynamic models use additional fields in the `phase` entry,
which are described in the documentation linked above
(see also @subpage sec-yaml-phase-thermo-models).

### Declaring the Elements {#sec-yaml-phase-declaring-elements}

In most cases, it is not necessary to specify the elements present in a
phase. If no `elements` field is present, elements will be added
automatically using the definitions of the standard chemical elements
based on the composition of the species present in the phase.

If non-standard elements such as isotopes need to be represented, or the
ordering of elements within the phase is important, the elements in the
phase may be declared in the optional `elements` entry.

If all of the elements to be added are either standard chemical elements or defined in
the <tt>[elements](@ref sec-yaml-elements)</tt> section of the current input file, the
elements can be specified as a list of element symbols. For example:

``` yaml
elements: [H, C, O, Ar]
```

To add elements from other top-level sections, from a different file, or
from multiple such sources, a list of single-key mappings can be used
where the key of each mapping specifies the source and the value is a
list of element names. The keys can be:

-   The name of a section within the current file.
-   The name of an input file and a section in that file, separated by a
    slash, for example `myfile.yaml/my_elements`. If a relative path is
    specified, the directory containing the current file is searched
    first, followed by the %Cantera data path.
-   The name `default` to reference the standard chemical elements.

**Example:**

``` yaml
elements:
- default: [C, H, Ar]
- isotopes: [O18]
- myelements.yaml/uranium: [U235, U238]
```

The order of the elements specified in the input file determines the
order of the elements in the phase when it is imported by %Cantera.

### Declaring the Species {#sec-yaml-phase-declaring-species}

If the species present in the phase corresponds to those species defined
in the `species` section of the input file, the `species` field may be
omitted, and those species will be added to the phase automatically. As
a more explicit alternative, the `species` field may be set to the
string `all`.

To include specific species from the `species` section of the input
file, the `species` entry can be a list of species names from that
section. For example:

``` yaml
species: [H2, O2, H2O]
```

If species are defined in multiple input file sections, the `species`
entry can be a list of single-key mappings, where the key of each
mapping specifies the source and the value is either the string `all` or
a list of species names. Each key can be either the name of a section
within the current input file or the name of a different file and a
section in that file, separated by a slash. If a relative path is
specified, the directory containing the current file is searched first,
followed by the %Cantera data path. Example:

``` yaml
species:
- species: [O2, N2]
- more_species: all
- subdir/myfile.yaml/species: [NO2, N2O]
```

The order of species specified in the input file determines the order of
the species in the phase when it is imported by %Cantera.

Species containing elements that are not declared within the phase may
be skipped by setting the `skip-undeclared-elements` field to `true`.
For example, to add all species from the `species` section that contain
only hydrogen or oxygen, the phase definition could contain:

``` yaml
phases:
- name: hydrogen-and-oxygen
  elements: [H, O]
  species: all
  skip-undeclared-elements: true
```

### Setting the Kinetics Model {#sec-yaml-phase-setting-kinetics}

The kinetics model to be used, if any, is specified in the `kinetics`
field. Supported model strings are:

-   `none`
-   <tt>[gas](@ref Cantera.GasKinetics)</tt>
-   <tt>[surface](@ref Cantera.InterfaceKinetics)</tt>
-   <tt>[edge](@ref Cantera.EdgeKinetics)</tt>

If omitted, no kinetics model will be used.

### Declaring the Reactions {#sec-yaml-phase-declaring-reactions}

If a kinetics model has been specified, reactions may be added to the
phase. By default, all reactions from the `reactions` section of the
input file will be added. Equivalently, the `reactions` entry may be
specified as the string `all`.

To disable automatic addition of reactions from the `reactions` section,
the `reactions` entry may be set to `none`. This may be useful if
reactions will be added programmatically after the phase is constructed.
The `reactions` field must be set to `none` if a kinetics model has been
specified but there is no `reactions` section in the input file.

To include only those reactions from the `reactions` section where all
of the species involved are declared as being in the phase, the
`reactions` entry can be set to the string `declared-species`.

To include reactions from multiple sections or other files, the
`reactions` entry can be given as a list of section names, for example:

``` yaml
reactions:
- OH_submechanism
- otherfile.yaml/C1-reactions
- otherfile.yaml/C2-reactions
```

To include reactions from multiple sections or other files while only
including reactions involving declared species, a list of single-key
mappings can be used, where the key is the section name (or file and
section name) and the value is either the string `all` or the string
`declared-species`. For example:

``` yaml
reactions:
- OH_submechanism: all
- otherfile.yaml/C1-reactions: all
- otherfile.yaml/C2-reactions: declared-species
```

To permit reactions containing third-body efficiencies for species not
present in the phase, the additional field
`skip-undeclared-third-bodies` may be added to the phase entry with the
value `true`.

### Setting the Transport Model {#sec-yaml-phase-setting-transport}

To enable transport property calculation, the transport model to be used
can be specified in the `transport` field. Supported models are:

-   `none`:
    No transport model
-   <tt>[high-pressure](@ref Cantera.HighPressureGasTransport)</tt>:
    A model for high-pressure gas transport properties based on a method
    of corresponding states
-   <tt>[ionized-gas](@ref Cantera.IonGasTransport)</tt>:
    A model implementing the Stockmayer-(n,6,4) model for transport of
    ions in a gas
-   <tt>[mixture-averaged](@ref Cantera.MixTransport)</tt>:
    The mixture-averaged transport model for ideal gases
-   <tt>[mixture-averaged-CK](@ref Cantera.MixTransport)</tt>:
    The mixture-averaged transport model for ideal gases, using
    polynomial fits corresponding to Chemkin-II
-   <tt>[multicomponent](@ref Cantera.MultiTransport)</tt>:
    The multicomponent transport model for ideal gases
-   <tt>[multicomponent-CK](@ref Cantera.MultiTransport)</tt>:
    The multicomponent transport model for ideal gases, using polynomial
    fits corresponding to Chemkin-II
-   <tt>[unity-Lewis-number](@ref Cantera.UnityLewisTransport)</tt>:
    A transport model for ideal gases, where diffusion coefficients for
    all species are set so that the Lewis number is 1
-   <tt>[water](@ref Cantera.WaterTransport)</tt>:
    A transport model for pure water applicable in both liquid and vapor phases

### Declaring Adjacent Phases {#sec-yaml-phase-declaring-adjacent}

For interface phases (surfaces and edges), the names of phases adjacent
to the interface can be specified, in which case these additional phases
can be automatically constructed when creating the interface phase. This
behavior is useful when the interface has reactions that include species
from one of these adjacent phases, since those phases must be known when
adding such reactions to the interface.

If the definitions of the adjacent phases are contained in the
`phases` section of the same input file as the interface,
they can be specified as a list of names:

``` yaml
adjacent: [gas, bulk]
```

Alternatively, if the adjacent phase definitions are in other sections
or other input files, they can be specified as a list of single-key
mappings where the key is the section name (or file and section name)
and the value is the phase name:

``` yaml
adjacent:
- {sectionname: gas} # a phase defined in a different section of the same YAML file
- {path/to/other-file.yaml/phases: bulk} # a phase defined in the 'phases' section
                                         # of a different YAML file
```

Since an interface kinetics mechanism is defined for the
lowest-dimensional phase involved in the mechanism, only
higher-dimensional adjacent phases should be specified. For example,
when defining a surface, adjacent bulk phases may be specified, but
adjacent edges must not.

### Setting the Initial State {#sec-yaml-phase-setting-state}

The state of a phase can be set using two properties to set the
thermodynamic state, plus the composition. This state is specified as a
mapping in the `state` field of the `phase` entry.

The thermodynamic state can be set in terms of two of the following
properties, with the valid property pairs depending on the phase model:

-   `temperature` or `T`
-   `pressure` or `P`
-   `enthalpy` or `H`
-   `entropy` or `S`
-   `int-energy`, `internal-energy` or `U`
-   `specific-volume` or `V`
-   `density` or `D`
-   `vapor-fraction` or `Q`

The thermodynamic state can be set using the following property pairs,
with some exceptions for phases where setting that property pair is not
implemented. All properties are on a per unit mass basis where relevant:

-   `T` and `P`
-   `T` and `D`
-   `T` and `V`
-   `H` and `P`
-   `U` and `V`
-   `S` and `V`
-   `S` and `P`
-   `S` and `T`
-   `P` and `V`
-   `U` and `P`
-   `V` and `H`
-   `T` and `H`
-   `S` and `H`
-   `D` and `P`

The composition can be set using one of the following fields, depending
on the phase type. The composition is specified as a mapping of species
names to values. Where necessary, the values will be automatically
normalized.

-   `mass-fractions` or `Y`
-   `mole-fractions` or `X`
-   `coverages`
-   `molalities` or `M`

**Examples:**

``` yaml
state:
  T: 300 K
  P: 101325 Pa
  X: {O2: 1.0, N2: 3.76}

state:
  density: 100 kg/m^3
  T: 298
  Y:
    CH4: 0.2
    C3H8: 0.1
    CO2: 0.7
```

For pure fluid phases, the temperature, pressure, and vapor fraction may
all be specified if and only if they define a consistent state.

### Examples

The following input file defines two equivalent gas phases including all
reactions and species defined in the input file. The species and
reaction data is not shown for clarity. In the second case, the phase
definition is simplified by having the elements added based on the
species definitions, taking the species definitions from the default
`species` section, and reactions from the default `reactions` section.

``` yaml
phases:
- name: gas1
  thermo: ideal-gas
  elements: [O, H, N, Ar]
  species: [H2, H, O, O2, OH, H2O, HO2, H2O2, N2, AR]
  kinetics: gas
  reactions: all
  transport: mixture-averaged
  state:
    T: 300.0
    P: 1.01325e+05
- name: gas2
  thermo: ideal-gas
  kinetics: gas
  transport: mixture-averaged
  state: {T: 300.0, 1 atm}

species:
- H2: ...
- H: ...
...
- AR: ...

reactions:
...
```

An input file defining an interface and its adjacent bulk phases, with
full species data not shown for clarity:

``` yaml
phases:
- name: graphite
  thermo: lattice
  species:
  - graphite-species: all
  state: {T: 300, P: 101325, X: {C6: 1.0, LiC6: 1e-5}}
  density: 2.26 g/cm^3

- name: electrolyte
  thermo: lattice
  species: [{electrolyte-species: all}]
  density: 1208.2 kg/m^3
  state:
    T: 300
    P: 101325
    X: {Li+(e): 0.08, PF6-(e): 0.08, EC(e): 0.28, EMC(e): 0.56}

- name: anode-surface
  thermo: ideal-surface
  adjacent: [graphite, electrolyte]
  kinetics: surface
  reactions: [graphite-anode-reactions]
  species: [{anode-species: all}]
  site-density: 1.0 mol/cm^2
  state: {T: 300, P: 101325}

graphite-species:
- name: C6
  ...
- name: LiC6
  ...

electrolyte-species:
- name: Li+(e)
  ...
- name: PF6-(e)
  ...
- name: EC(e)
  ...
- name: EMC(e)
  ...

anode-species:
- name: (int)
  ...

graphite-anode-reactions:
- equation: LiC6 <=> Li+(e) + C6
  rate-constant: [5.74, 0.0, 0.0]
  beta: 0.4
```

## Phase API Reference

A `phase` is a mapping that contains definitions for the elements,
species, and optionally reactions that can take place in that phase.

The fields of a `phase` entry are:

<b>`name`</b>

String identifier used for the phase. Required.

<b>`elements`</b>

Specification for the elements present in the phase. This can be:

-   Omitted, in which case the standard elements will be added as
    needed by the species included in the phase.
-   A list of element symbols, which can be either defined in the
    `elements` section of the file or taken from the standard
    elements.
-   A list of single-key mappings of section names to lists of
    element symbols. These sections can be in the same file as the
    phase definition, or from another file if written as
    `file-path/sectionname`. If a relative path is specified, the
    directory containing the current file is searched first,
    followed by the %Cantera data path. Standard elements can be
    included by referencing the fictitious section `default`.

See also @ref sec-yaml-phase-declaring-elements.

<b>`species`</b>

Specification for the species present in the phase. This can be:

-   a list of species that appear in the `species` section of the
    file.
-   The string `all`, to indicate that all species in the `species`
    section should be included. This is the default if no `species`
    entry is present.
-   A list of single-key mappings of section names to either the
    string `all` or a list of species names. These sections can be
    in the same file as the phase definition, or from another file
    if written as `file-path/sectionname`. If a relative path is
    specified, the directory containing the current file is searched
    first, followed by the %Cantera data path.

See also @ref sec-yaml-phase-declaring-species. Species may be skipped
depending on the setting of the `skip-undeclared-elements` option.

<b>`skip-undeclared-elements`</b>

If set to `true`, do not add species that contain elements that are
not explicitly included in the phase. The default is `false`, where
the presence of such species is considered an error.

<b>`skip-undeclared-third-bodies`</b>

If set to `true`, ignore third body efficiencies for species that
are not defined in the phase. The default is `false`, where the
presence of such third body specifications is considered an error.

<b>`state`</b>

A mapping specifying the thermodynamic state. See @ref sec-yaml-phase-setting-state.

<b>`adjacent-phases`</b>

For interface phases, specification of adjacent phases that
participate in reactions on the interface. This can be:

-   a list of phase names that appear in the `phases` section of the
    file.
-   A list of single-key mappings of section names to a list of
    phase names. These sections can be in the same file as the
    current phase definition, or from another file if written as
    `file-path/section-name`. If a relative path is specified, the
    directory containing the current file is searched first,
    followed by the %Cantera data path.

See also @ref sec-yaml-phase-declaring-adjacent.

<b>`thermo`</b>

String specifying the phase thermodynamic model to be used.
For a list of supported model strings, see @ref sec-yaml-phase-setting-thermo.

<b>`kinetics`</b>

String specifying the kinetics model to be used. For supported model
strings, see sec-yaml-phase-setting-kinetics.

<b>`reactions`</b>

Source of reactions to include in the phase, if a kinetics model has
been specified. This can be:

-   The string `all`, which indicates that all reactions from the
    `reactions` section of the file should be included. This is the
    default if no `reactions` entry is present.
-   The string `declared-species`, which indicates that all
    reactions from the `reactions` section involving only species
    present in the phase should be included.
-   The string `none`, which indicates that no reactions should be
    added. This can be used if reactions will be added
    programmatically after the phase is constructed.
-   A list of sections from which to include reactions. These
    sections can be in the same file as the phase definition, or
    from another file if written as `file-path/sectionname`. If a
    relative path is specified, the directory containing the current
    file is searched first, followed by the %Cantera data path.
-   A list of single-key mappings of section names to rules for
    adding reactions, where for each section name, that rule is
    either `all` or `declared-species` and is applied as described
    above.

See also @ref sec-yaml-phase-declaring-reactions.

<b>`Motz-Wise`</b>

Boolean indicating whether the Motz-Wise correction should be
applied to sticking reactions. Applicable only to interface phases.
The default is `false`. The value set at the phase level may be
overridden on individual reactions.

<b>`transport`</b>

String specifying the transport model to be used. For a list of supported
model strings, see @ref sec-yaml-phase-setting-transport.
