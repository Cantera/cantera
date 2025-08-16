# YAML Input File Reference

## [General Structure](./general)

This section of the documentation describes the [sections](sec-yaml-sections) of a YAML
input file and how dimensional quantities with different [units](sec-yaml-units) can
be written.

## [Phase Definitions](./phases)

This section describes how to define a [phase](sec-yaml-phases), which includes
specifying the [species names](sec-yaml-phase-species), the [thermodynamic
model](sec-yaml-phase-thermo), the [kinetics model](sec-yaml-phase-kinetics), the
[transport model](sec-yaml-phase-transport), and the [initial
state](sec-yaml-setting-state), as well as other, optional properties.

## [Elements](./elements)

YAML element definitions are needed only when defining custom elements that are not
standard chemical elements, or defining specific isotopes.

## [Species](./species)

Species definitions specify the name and composition of a species, and include entries
defining parameters needed for [species thermo](sec-yaml-species-thermo), [equation of
state](sec-yaml-species-eos), [transport property](sec-yaml-species-transport), and
[coverage dependency](sec-yaml-species-coverage) models.

## [Reactions](./reactions)

YAML reaction definitions include specification of common elements such as the reaction
equation and [efficiencies](sec-yaml-efficiencies), as well as parameters specific to
the type of [rate parameterization](sec-yaml-rate-types).
See also the [electron collision data format](reactions.html#sec-yaml-electron-collisions)
used in plasma-phase simulations.

## Mechanism Conversion

Cantera provides scripts for converting mechanisms from the Chemkin
([`ck2yaml`](./ck2yaml)), LXCat ([`lxcat2yaml`](./lxcat2yaml)), CTI
([`cti2yaml`](./cti2yaml)), and CTML/XML ([`ctml2yaml`](/yaml/ctml2yaml)) formats to
YAML. Cantera also provides the [`yaml2ck`](/yaml/yaml2ck) script for converting YAML
input files to the Chemkin format.


```{toctree}
:hidden:

general
phases
elements
species
reactions
ck2yaml
cti2yaml
ctml2yaml
yaml2ck
lxcat2yaml
```
