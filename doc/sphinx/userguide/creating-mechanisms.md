# Creating YAML Mechanism Files from Scratch

Virtually every Cantera simulation involves one or more phases of matter. Depending on
the calculation being performed, it may be necessary to evaluate thermodynamic
properties, transport properties, and/or reaction rates for the phase(s) present. Before
the properties can be evaluated, each phase must be defined, meaning that the models to
use to compute its properties and reaction rates must be specified, along with any
parameters the models require.

Because the amount of data required can be quite large, this data is imported from a
YAML file that can be read by the application, so that a given phase model can be
re-used for other simulations. This guide describes how to write such files to define
phases and interfaces for use in Cantera simulations.

```{seealso}
See the [](input-tutorial) for an introduction to the YAML syntax used by Cantera and
a description of how dimensional values are handled.
```

## Phases

For each phase that appears in a problem, a corresponding entry should be present in the
`phases` section of the input file. The phase entry specifies the elements and species
present in that phase, and the models to be used for computing thermodynamic, kinetic,
and transport properties.

### Naming the Phase

The `name` entry is a string that identifies the phase. It must be unique within the
file among all phase definitions of any type. Phases are referenced by name when
importing them. The `name` is also used to identify the phase within multiphase mixtures
or at phase boundaries.

### Setting the Thermodynamic Model

The thermodynamic model used to represent a phase is specified in the `thermo` field. A
[complete list of supported models](sec-yaml-phase-thermo-models) can be found in the
YAML Input File Reference. Some thermodynamic models use additional fields in the phase
entry, which are described in the linked documentation.

### Declaring the Elements

In most cases, it is not necessary to specify the elements present in a phase. If no
`elements` field is present, elements will be added automatically using the definitions
of the standard chemical elements based on the composition of the species present in the
phase.

If non-standard elements such as isotopes need to be represented, or the ordering of
elements within the phase is important, the elements in the phase may be declared in the
optional `elements` entry.

If all of the elements to be added are either standard chemical elements or defined in
the [`elements`](sec-yaml-guide-elements) section of the current input file, the
elements can be specified as a list of element symbols. For example:

``` yaml
phases:
- name: my-mechanism
  elements: [H, C, O, Ar]
  ...
```

To add elements from other top-level sections, from a different file, or from multiple
such sources, a list of single-key mappings can be used where the key of each mapping
specifies the source and the value is a list of element names. The keys can be:

- The name of a section within the current file.
- The name of an input file and a section in that file, separated by a slash, for
  example `myfile.yaml/my_elements`. If a relative path is specified, the directory
  containing the current file is searched first, followed by the Cantera data path.
- The name `default` to reference the standard chemical elements.

Example:

```yaml

my-isotopes:
- name: O18
  atomic-weight: 17.9991603
phases:
- name: my-phase
  elements:
  - default: [C, H, Ar]
  - my-isotopes: [O18]
  - myelements.yaml/uranium: [U235, U238]
  species: ...
  ...
```

The order of the elements specified in the input file determines the order of the
elements in the phase when it is imported by Cantera.

### Declaring the Species

If the species present in the phase corresponds to those species defined in the
[`species`](sec-yaml-guide-species) section of the input file, the `species` field may
be omitted, and those species will be added to the phase automatically. As a more
explicit alternative, the `species` field may be set to the string `all`.

To include specific species from the `species` section of the input file, the `species`
entry can be a list of species names from that section. For example:

```yaml
phases:
- name: my-phase
  species: [H2, O2, H2O]
  ...
```

If species are defined in multiple input file sections, the `species` entry can be a
list of single-key mappings, where the key of each mapping specifies the source and the
value is either the string `all` or a list of species names. Each key can be either the
name of a section within the current input file or the name of a different file and a
section in that file, separated by a slash. If a relative path is specified, the
directory containing the current file is searched first, followed by the Cantera data
path. Example:

```yaml
phases:
- name: my-phase
  species:
  - species: [O2, N2]
  - more_species: all
  - subdir/myfile.yaml/species: [NO2, N2O]
  ...
```

The order of species specified in the input file determines the order of the species in
the phase when it is imported by Cantera.

Species containing elements that are not declared within the phase may be skipped by
setting the `skip-undeclared-elements` field to `true`. For example, to add all species
from the `species` section that contain only hydrogen or oxygen, the phase definition
could contain:

```yaml
phases:
- name: hydrogen-and-oxygen
  elements: [H, O]
  species: all
  skip-undeclared-elements: true
  ...
```

### Setting the Kinetics Model

The kinetics model to be used, if any, is specified in the `kinetics` field. Supported
kinetics models are `bulk`, `surface`, and `edge`, depending on the dimensionality of
the phase. If omitted, no kinetics model will be used. For additional details, see the
[list of supported models](sec-yaml-phase-kinetics) in the YAML Input File Reference.

### Declaring the Reactions

If a kinetics model has been specified, reactions may be added to the phase. By default,
all reactions from the `reactions` section of the input file will be added.
Equivalently, the `reactions` entry may be specified as the string `all`.

To disable automatic addition of reactions from the `reactions` section, the `reactions`
entry may be set to `none`. This may be useful if reactions will be added
programmatically after the phase is constructed. The `reactions` field must be set to
`none` if a kinetics model has been specified but there is no `reactions` section in the
input file.

To include only those reactions from the `reactions` section where all of the species
involved are declared as being in the phase, the `reactions` entry can be set to the
string `declared-species`.

To include reactions from multiple sections or other files, the `reactions` entry can be
given as a list of section names, for example:

```yaml
phases:
- name: my-phase
  ...
  reactions:
  - OH_submechanism
  - otherfile.yaml/C1-reactions
  - otherfile.yaml/C2-reactions
  ...
```

To include reactions from multiple sections or other files while only including
reactions involving declared species, a list of single-key mappings can be used, where
the key is the section name (or file and section name) and the value is either the
string `all` or the string `declared-species`. For example:

```yaml
phases:
- name: my-phase
  ...
  reactions:
  - OH_submechanism: all
  - otherfile.yaml/C1-reactions: all
  - otherfile.yaml/C2-reactions: declared-species
  ...
```

To permit reactions containing third-body efficiencies for species not present in the
phase, the additional field `skip-undeclared-third-bodies` may be added to the phase
entry with the value `true`.

### Setting the Transport Model

To enable transport property calculation, the transport model to be used can be
specified in the `transport` field. A [complete list of supported models](sec-yaml-phase-transport)
can be found in the YAML Input File Reference. For most transport models, additional
parameters need to be specified within each species definition.

### Declaring Adjacent Phases

For interface phases (surfaces and edges), the names of phases adjacent to the interface
can be specified, in which case these additional phases can be automatically constructed
when creating the interface phase. This behavior is useful when the interface has
reactions that include species from one of these adjacent phases, since those phases
must be known when adding such reactions to the interface.

If the definitions of the adjacent phases are contained in the `phases` section of the
same input file as the interface, they can be specified as a list of names:

```yaml
phases:
- name: my-surface-phase
  ...
  adjacent: [gas, bulk]
  ...
- name: gas
  ...
- name: bulk
  ...
```

Alternatively, if the adjacent phase definitions are in other sections or other input
files, they can be specified as a list of single-key mappings where the key is the
section name (or file and section name) and the value is the phase name:

```yaml
phases:
- name: my-surface-phase
  ...
  adjacent:
  - {my-other-phases: gas}
  # a phase defined in the 'phases' section of a different YAML file
  - {path/to/other-file.yaml/phases: bulk}
  ...
my-other-phases:
- name: gas
  ...
```

Since an interface kinetics mechanism is defined for the lowest-dimensional phase
involved in the mechanism, only higher-dimensional adjacent phases should be specified.
For example, when defining a surface, adjacent bulk phases may be specified, but
adjacent edges must not.

### Setting the Initial State

The initial state of a phase can be set using two properties to set the thermodynamic
state, plus the composition. This state is specified as a mapping in the `state` field
of `phase` entry.

The thermodynamic state can be set by specifying two properties, such as temperature
(`T`) and pressure (`P`) or internal energy (`U`) and density (`D`). The full list of
[property names and valid combinations](sec-yaml-setting-state) can be found in the YAML
Input File Reference. In addition the composition can be set by providing a mapping that
gives the mass fractions (`X`), mole fractions (`Y`), `coverages` (for surface phases),
or molalities (`M`, for certain solution models). Where necessary, the values will be
automatically normalized.

Some examples of setting the state:

```yaml
phases:
- name: my-phase
  ...
  state:
    T: 300 K
    P: 101325 Pa
    X: {O2: 1.0, N2: 3.76}
- name: my-other-phase
  ...
  state:
    density: 100 kg/m^3
    T: 298
    Y:
      CH4: 0.2
      C3H8: 0.1
      CO2: 0.7
```

For pure fluid phases, the temperature, pressure, and vapor fraction may all be
specified if and only if they define a consistent state.

### Examples

The following input file defines two equivalent gas phases including all reactions and
species defined in the input file. The species and reaction data is not shown for
clarity. In the second case, the phase definition is simplified by having the elements
added based on the species definitions, taking the species definitions from the default
`species` section, and reactions from the default `reactions` section.

```yaml
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

An input file defining an interface and its adjacent bulk phases, with full species data
not shown for clarity:

```yaml
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

(sec-yaml-guide-species)=
## Species

For each species declared as part of a phase description, a species definition is
required to describe the composition, thermodynamic, and transport parameters of that
species.

The default location for species entries is in the `species` section of the input file.
Species defined in this section will automatically be considered for addition to phases
defined in the same file. Species can be defined in other sections of the input file or
in other input files, and these species definitions can be used in phase definitions by
explicitly referencing the section name.

### Species Name

The name of a species is given in the `name` field of a `species` entry. Names may
include almost all printable characters, with the exception of spaces. The use of some
characters such as `[`, `]`, and `,` may require that species names be enclosed in
quotes when written in YAML. Some valid species names given in a YAML list include:

```yaml
[CH4, methane, argon_2+, "C[CH2]", CH2(singlet), "H2O,l"]
```

### Elemental Composition

The elemental composition of a species is specified as a mapping in the `composition`
entry.

For gaseous species, the elemental composition is well-defined, since the species
represent distinct molecules. For species in solid or liquid solutions, or on surfaces,
there may be several possible ways of defining the species. For example, an aqueous
species might be defined with or without including the water molecules in the solvation
cage surrounding it.

For surface species, it is possible for the `composition` mapping to be empty, in which
case the species is composed of nothing, and represents an empty surface site. This can
also be done to represent vacancies in solids. A charged vacancy can be defined to be
composed solely of electrons.

The number of atoms of an element must be non-negative, except for the special "element"
`E` that represents an electron.

Examples:

```yaml
composition: {C: 1, O: 2}  # carbon dioxide
composition: {Ar: 1, E: -2}  # Ar++
composition: {Y: 1, Ba: 2, Cu: 3, O: 6.5}  # stoichiometric YBCO
composition: {}  # A surface species representing an empty site
```

### Thermodynamic Properties

In addition to the thermodynamic model used at the phase level for computing properties,
parameterizations are usually required for the enthalpy, entropy, and specific heat
capacities of individual species under standard conditions. These parameterizations are
provided in the `thermo` field of each `species` entry, with the type of
parameterization specified by the `model` field. A
[complete list of parameterizations](sec-yaml-species-thermo) and their model-specific
fields can be found in the YAML Input File Reference.

An example `thermo` field using the 7-coefficient NASA polynomials in two temperature
regions:

```yaml
thermo:
  model: NASA7
  temperature-ranges: [200.0, 1000.0, 3500.0]
  data:
  - [5.14987613, -0.0136709788, 4.91800599e-05, -4.84743026e-08, 1.66693956e-11,
    -1.02466476e+04, -4.64130376]
  - [0.074851495, 0.0133909467, -5.73285809e-06, 1.22292535e-09, -1.0181523e-13,
    -9468.34459, 18.437318]
```

### Species Equation of State

For some phase thermodynamic models, additional equation of state parameterizations are
needed for each species. This information is provided in the `equation-of-state` field
of each `species` entry, with the type of parameterization used specified by the `model`
field of the `equation-of-state` field. A
[complete list of equation of state parameterizations](sec-yaml-species-eos) and their
model-specific fields can be found in the YAML Input File Reference.

(sec-yaml-guide-species-transport)=
### Species Transport Coefficients

Transport-related parameters for each species are needed in order to calculate transport
properties of a phase. These parameters are provided in the `transport` field of each
`species` entry, with the type of the parameterization used specified by the `model`
field of the `transport` field. The only model type specifically handled is `gas`. The
parameters used depend on the transport model specified at the phase level. The full set
of possible parameters is described in the {ref}`API documentation
<sec-yaml-species-transport>`.

An example of a `transport` entry for a gas-phase species:

```yaml
transport:
  model: gas
  geometry: linear
  well-depth: 107.4
  diameter: 3.458
  polarizability: 1.6
  rotational-relaxation: 3.8
```

(sec-yaml-guide-reactions)=
## Reactions

Cantera supports a number of different types of reactions, including several types of
homogeneous reactions, surface reactions, and electrochemical reactions. The reaction
entries for all reaction types some common features. These general fields of a reaction
entry are described first, followed by fields used for specific reaction types.


### The Reaction Equation

The reaction equation, specified in the `equation` field of the reaction entry,
determines the reactant and product stoichiometry. All tokens (species names,
stoichiometric coefficients, `+`, and `<=>`) in the reaction equation must be separated
with spaces. Some examples of correctly and incorrectly formatted reaction equations are
shown below:

```yaml
- equation: 2 CH2 <=> CH + CH3  # OK
- equation: 2 CH2<=>CH + CH3  # error - spaces required around '<=>''
- equation: 2CH2 <=> CH + CH3  # error - space required between '2' and 'CH2'
- equation: CH2 + CH2 <=> CH + CH3  # OK
- equation: 2 CH2 <=> CH+CH3  # error - spaces required around '+'
```

Whether the reaction is reversible or not is determined by the form of the equality sign
in the reaction equation. If either `<=>` or `=` is found, then the reaction is regarded
as reversible, and the reverse rate will be computed based on the equilibrium constant.
If, on the other hand, `=>` is found, the reaction will be treated as irreversible.

### Reaction type

The type of the rate coefficient parameterization may be specified in the `type` field
of the `reaction` entry. The [fields for specific reaction types](sec-yaml-reactions)
and additional parameters defining the rate constant for each of these reaction types
are described in the YAML Input File Reference.

The default parameterization is `elementary`. Reactions involving surface species are
automatically identified as [`interface`](sec-yaml-interface-Arrhenius) reactions,
reactions involving surface species with the `type` specified as `Blowers-Masel` are
treated as [`interface-Blowers-Masel`](sec-yaml-interface-Blowers-Masel), and reactions
involving charge transfer are automatically identified as
[`electrochemical`](sec-yaml-electrochemical-reaction) reactions.

### Arrhenius Expressions

Most reaction types in Cantera are parameterized by one or more modified Arrhenius
expressions, such as

$$  A T^b e^{-E_a / RT}  $$

where $A$ is the pre-exponential factor, $T$ is the temperature, $b$ is the temperature
exponent, $E_a$ is the activation energy, and $R$ is the gas constant. Rates in this
form can be written as YAML mappings. For example:

```yaml
{A: 1.0e13, b: 0, E: 7.3 kcal/mol}
```

The units of $A$ can be specified explicitly if desired. If not specified, they will be
determined based on the `quantity`, `length`, and `time` units specified in the
governing `units` fields. Since the units of $A$ depend on the reaction order, the units
of each reactant concentration (dependent on phase type and dimensionality), and the
units of the rate of progress (different for homogeneous and heterogeneous reactions),
it is usually best not to specify units for $A$, in which case they will be computed
taking all of these factors into account.

```{note}
If $b \ne 0$, then the term $T^b$ should have units of $\mathrm{K}^b$, which would
change the units of $A$. This is not done, however, so the units associated with $A$
are really the units for $k_f$. One way to formally express this is to replace $T^b$
by the non-dimensional quantity $[T/(1\;\mathrm{K})]^b$.
```

The key `E` is used to specify $E_a$.

(sec-yaml-reaction-options)=
### Duplicate Reactions

When a reaction is imported into a phase, it is checked to see that it is not a
duplicate of another reaction already present in the phase, and normally an error
results if a duplicate is found. But in some cases, it may be appropriate to include
duplicate reactions, for example if a reaction can proceed through two distinctly
different pathways, each with its own rate expression. Another case where duplicate
reactions can be used is if it is desired to implement a reaction rate coefficient of
the form:

$$  k_f(T) = \sum_{n=1}^{N} A_n T^{b_n} \exp(-E_n/RT)  $$

While Cantera does not provide such a form for reaction rates, it can be implemented by
defining $N$ duplicate reactions, and assigning one rate coefficient in the sum to each
reaction. By adding the field:

```yaml
duplicate: true
```

to a reaction entry, then the reaction not only *may* have a duplicate, it *must*. Any
reaction that specifies that it is a duplicate, but cannot be paired with another
reaction in the phase that qualifies as its duplicate generates an error.

### Negative Pre-exponential Factors

If some of the terms in the above sum have negative $A_n$, this scheme fails, since
Cantera normally does not allow negative pre-exponential factors. But if there are
duplicate reactions such that the total rate is positive, then the fact that negative
$A$ parameters are acceptable can be indicated by adding the field:

```yaml
negative-A: true
```

### Reaction Orders

Explicit reaction orders different from the stoichiometric coefficients are sometimes
used for non-elementary reactions. For example, consider the global reaction:

$$  \mathrm{C_8H_{18} + 12.5 O_2 \rightarrow 8 CO_2 + 9 H_2O}  $$

the forward rate constant might be given as {cite:p}`westbrook1981`:

$$
k_f = 4.6 \times 10^{11} [\mathrm{C_8H_{18}}]^{0.25} [\mathrm{O_2}]^{1.5}
      \exp\left(\frac{30.0\,\mathrm{kcal/mol}}{RT}\right)
$$

This reaction could be defined as:

```yaml
- equation: C8H18 + 12.5 O2 => 8 CO2 + 9 H2O
  rate-constant: {A: 4.6e11, b: 0.0, Ea: 30.0 kcal/mol}
  orders: {C8H18: 0.25, O2: 1.5}
```

Special care is required in this case since the units of the pre-exponential factor
depend on the sum of the reaction orders, which may not be an integer.

Note that you can change reaction orders only for irreversible reactions.

#### Negative Reaction Orders

Normally, reaction orders are required to be positive. However, in some cases negative
reaction orders provide better fits for experimental data. In these cases, the default
behavior may be overridden by adding the `negative-orders` field to the reaction entry.
For example:

```yaml
- equation: C8H18 + 12.5 O2 => 8 CO2 + 9 H2O
  rate-constant: {A: 4.6e11, b: 0.0, Ea: 30.0 kcal/mol}
  orders: {C8H18: -0.25, O2: 1.75}
  negative-orders: true
```

#### Non-reactant Orders

Some global reactions could have reactions orders for non-reactant species. In this
case, the `nonreactant-orders` field must be added to the reaction entry:

```yaml
- equation: C8H18 + 12.5 O2 => 8 CO2 + 9 H2O
  rate-constant: {A: 4.6e11, b: 0.0, Ea: 30.0 kcal/mol}
  orders: {C8H18: -0.25, CO: 0.15}
  negative-orders: true
  nonreactant-orders: true
```

(sec-yaml-guide-elements)=
## Elements

Cantera provides built-in definitions for the chemical elements, including values for
their atomic weights taken from IUPAC / CIAAW. These elements can be used by specifying
the corresponding atomic symbols when specifying the composition of species.

In order to give a name to a particular isotope or a virtual element representing a
surface site, a custom `element` entry can be used. The default location for `element`
entries is the `elements` section of the input file. Elements defined in this section
will automatically be considered for addition to phases defined in the same file.
Elements can be defined in other sections of the input file if those sections are named
explicitly in the `elements` field of the phase definition.

An element entry has the following fields:

- `symbol`: The symbol to be used for the element, for example when specifying the
  composition of a species.
- `atomic-weight`: The atomic weight of the element, in unified atomic mass units
  (dalton)
- `atomic-number`: The atomic number of the element. Optional.
- `entropy298`: The standard molar entropy of the element at 298.15 K. Optional.

An example `elements` section:

```yaml
elements:
- symbol: C13
  atomic-weight: 13.003354826
  atomic-number: 6
- symbol: O-18
  atomic-weight: 17.9991603
```
