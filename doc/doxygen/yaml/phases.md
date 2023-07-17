# Phase Definitions {#sec-yaml-phases}

[TOC]

## Fields

A `phase` is a mapping that contains definitions for the elements,
species, and optionally reactions that can take place in that phase.

The fields of a `phase` entry are:

`name`

String identifier used for the phase. Required.

`elements`

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

`species`

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

Species may be skipped depending on the setting of the
`skip-undeclared-elements` option.

`skip-undeclared-elements`

If set to `true`, do not add species that contain elements that are
not explicitly included in the phase. The default is `false`, where
the presence of such species is considered an error.

`skip-undeclared-third-bodies`

If set to `true`, ignore third body efficiencies for species that
are not defined in the phase. The default is `false`, where the
presence of such third body specifications is considered an error.

`state`

A mapping specifying the thermodynamic state. See
`sec-yaml-setting-state`{.interpreted-text role="ref"}.

`adjacent-phases`

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

`thermo`

String specifying the phase thermodynamic model to be used.
Supported model strings are:

-   <tt>@ref sec-yaml-binary-solution-tabulated</tt>
-   <tt>@ref sec-yaml-compound-lattice</tt>
-   <tt>@ref sec-yaml-coverage-dependent-surface</tt>
-   <tt>@ref sec-yaml-Debye-Huckel</tt>
-   <tt>@ref sec-yaml-edge</tt>
-   <tt>@ref sec-yaml-electron-cloud</tt>
-   <tt>@ref sec-yaml-fixed-stoichiometry</tt>
-   <tt>@ref sec-yaml-HMW-electrolyte</tt>
-   <tt>@ref sec-yaml-ideal-condensed</tt>
-   <tt>@ref sec-yaml-ideal-gas</tt>
-   <tt>@ref sec-yaml-ideal-molal-solution</tt>
-   <tt>@ref sec-yaml-ideal-solution-VPSS</tt>
-   <tt>@ref sec-yaml-ideal-surface</tt>
-   <tt>@ref sec-yaml-ions-from-neutral-molecule</tt>
-   <tt>@ref sec-yaml-lattice</tt>
-   <tt>@ref sec-yaml-liquid-water-IAPWS95</tt>
-   <tt>@ref sec-yaml-Margules</tt>
-   <tt>@ref sec-yaml-Maskell-solid-solution</tt>
-   <tt>@ref sec-yaml-Peng-Robinson</tt>
-   <tt>@ref sec-yaml-plasma</tt>
-   <tt>@ref sec-yaml-pure-fluid</tt>
-   <tt>@ref sec-yaml-Redlich-Kister</tt>
-   <tt>@ref sec-yaml-Redlich-Kwong</tt>

`kinetics`

String specifying the kinetics model to be used. Supported model
strings are:

-   `none`
-   `gas`: see @ref Cantera.GasKinetics
-   `surface`: see @ref Cantera.InterfaceKinetics
-   `edge`: see @ref Cantera.EdgeKinetics

`reactions`

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

`Motz-Wise`

Boolean indicating whether the Motz-Wise correction should be
applied to sticking reactions. Applicable only to interface phases.
The default is `false`. The value set at the phase level may be
overridden on individual reactions.

`transport`

String specifying the transport model to be used. Supported model
strings are:

-   `none`
-   `high-pressure`: see @ref Cantera.HighPressureGasTransport
-   `ionized-gas`: see @ref Cantera.IonGasTransport
-   `mixture-averaged`: see @ref Cantera.MixTransport
-   `mixture-averaged-CK`: see @ref Cantera.MixTransport
-   `multicomponent`: see @ref Cantera.MultiTransport
-   `multicomponent-CK`: see @ref Cantera.MultiTransport
-   `unity-Lewis-number`: see @ref Cantera.UnityLewisTransport
-   `water`: see @ref Cantera.WaterTransport

## Setting the State {#sec-yaml-setting-state}

The state of a `phase` can be set using two properties to set the
thermodynamic state, plus the composition.

The composition can be set using one of the following fields, depending
on the phase type. The composition is specified as a mapping of species
names to values. Where necessary, the values will be automatically
normalized.

-   `mass-fractions` or `Y`
-   `mole-fractions` or `X`
-   `coverages`
-   `molalities` or `M`

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

The following synonyms are also implemented for use in any of the pairs:

-   `temperature`, `T`
-   `pressure`, `P`
-   `enthalpy`, `H`
-   `entropy`, `S`
-   `int-energy`, `internal-energy`, `U`
-   `specific-volume`, `V`
-   `density`, `D`
