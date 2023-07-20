# %Cantera YAML Documentation {#sec-yaml-documentation}

This short guide describes %Cantera YAML input files that define phases and
interfaces for use in %Cantera simulations. Each link below represents an entry
point to the @subpage sec-yaml-api; while you certainly can read them in order,
you can also jump to whichever section addresses your current needs. If you need
tips on troubleshooting the YAML file syntax rules, please look at the
@subpage sec-yaml-format-tutorial.

-   <b>@ref sec-yaml-general</b> -
    Structure of a %Cantera YAML input file.
-   <b>@ref sec-yaml-phases</b> -
    For each phase or interface that appears in a problem, a corresponding
    entry should be present in the input file(s). Phase entries specify lists
    of elements, species and reactions, as well as other relevant parameters.
-   <b>@ref sec-yaml-species</b> -
    For each species declared as part of a phase description, thermodynamic
    properties and other data need to be defined.
-   <b>@ref sec-yaml-reactions</b> -
    %Cantera supports a number of different types of reactions, including
    several types of homogeneous reactions, surface reactions, and
    electrochemical reactions. For each, there is a corresponding entry
    type.

> **Note:** %Cantera provides @ref sec-yaml-conversion to generate YAML files
> from Chemkin input data as well as legacy %Cantera formats.
