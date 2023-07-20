# %Cantera YAML Documentation {#sec-yaml-guide}

This short guide describes %Cantera YAML input files that define phases and
interfaces for use in %Cantera simulations. Each link below represents a
standalone module; while you certainly can read them in order, you can
also jump to whichever section addresses your current needs. If you need
tips on troubleshooting the YAML file syntax rules, please look at the
[YAML Format Tutorial](\ref sec-yaml-format-tutorial).

-   **[Phases and Interfaces:](\ref sec-yaml-phases)**
    For each phase or interface that appears in a problem, a corresponding
    entry should be present in the input file(s). Phase entries specify lists
    of elements, species and reactions, as well as other relevant parameters.
-   **[Species Definitions:](\ref sec-yaml-species)**
    For each species declared as part of a phase description, thermodynamic
    properties and other data need to be defined.
-   **[Reaction Definitions:](\ref sec-yaml-reactions)**
    %Cantera supports a number of different types of reactions, including
    several types of homogeneous reactions, surface reactions, and
    electrochemical reactions. For each, there is a corresponding entry
    type.

**Additional Information**

-   <b>@subpage sec-yaml-format-tutorial</b> -
    This module describes the basics of the YAML format as used by %Cantera,
    how dimensional values are represented, and how to understand error
    messages that occur while reading input files.
-   <b>@subpage sec-yaml-api</b> -
    Comprehensive documentation of the %Cantera YAML format, containing the
    specification for each of the entry types discussed previously.
