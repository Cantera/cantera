# YAML Format Reference {#sec-yaml-api}

All calculations in %Cantera require an input file to describe the properties of the
relevant phase(s) of matter. %Cantera uses the
[YAML data language](https://yaml.org/spec/1.2/spec.html#Introduction) to specify
information on thermodynamics, kinetics, and transport in a single file.

-   @subpage sec-yaml-general
-   @subpage sec-yaml-units
-   @subpage sec-yaml-phases
    -   @ref sec-yaml-phase-thermo-models
-   @subpage sec-yaml-elements
-   @subpage sec-yaml-species
    -   @ref sec-yaml-species-thermo-models
    -   @ref sec-yaml-species-eos-models
    -   @ref sec-yaml-species-transport-models
    -   @ref sec-yaml-species-coverage-models
-   @subpage sec-yaml-reactions
    -   @ref sec-yaml-reaction-rate-blocks
    -   @ref sec-yaml-reaction-models

For an introductory tutorial of the YAML data language, refer to the
@ref sec-yaml-format-tutorial.

## Conversion Utilities {#sec-yaml-conversion}

%Cantera provides the following conversion utilities to obtain YAML input from other
input file formats:

-   <tt>[ck2yaml](https://cantera.org/tutorials/ck2yaml-tutorial.html)</tt> Conversion
    from Chemkin-format files to YAML
-   <tt>[cti2yaml](https://cantera.org/tutorials/legacy2yaml.html#cti2yaml)</tt>
    Conversion of legacy %Cantera CTI files to YAML
-   <tt>[ctml2yaml](https://cantera.org/tutorials/legacy2yaml.html#cti2yaml)</tt>
    Conversion of legacy %Cantera XML (CTML) files to YAML
