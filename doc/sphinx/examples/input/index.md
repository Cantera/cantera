# Example Input Files

The following example input files, demonstrate how some of the features of Cantera's
YAML format can be used to define phases, species, and reactions. These files are a
subset of the input files included with Cantera.

[`h2o2.yaml`](h2o2.md)
: A basic gas phase reaction mechanism. Demonstrates multiple phase definitions using
  the same set of species, specification of transport properties, and parameters for
  the Redlich-Kwong equation of state.

[`diamond.yaml`](diamond.md)
: A heterogenous reaction mechanism. Demonstrates defining surfaces and adjacent bulk
  phases and reading species definitions from other input files. This input file is used
  in the example [`diamond_cvd.py`](/examples/python/kinetics/diamond_cvd).

[`lithium_ion_battery.yaml`](lithium_ion_battery.md)
: An input file for an LCO/graphite lithium-ion battery. Demonstrates multiple surface
  phase definitions, input for the
  [`binary-solution-tabulated`](sec-yaml-binary-solution-tabulated) thermo model, and
  electrochemical reactions. This input file is used in the example
  [`lithium_ion_battery.py`](/examples/python/kinetics/lithium_ion_battery).

[`example_data/covdepsurf.yaml`](covdepsurf.md)
: An input file for a surface phase where the non-ideal interactions between species are
  accounted for in the calculation of the enthalpy and entropy. Demonstrates the
  [`coverage-dependent-surface`](sec-yaml-coverage-dependent-surface) model and the
  corresponding [`coverage-dependencies`](sec-yaml-species-coverage) field of the
  `species` entry. This input file is used in the example
  [`coverage_dependent_surf.py`](/examples/python/thermo/coverage_dependent_surf).

[`example_data/oxygen-plasma-itikawa.yaml`](oxygen-plasma-itikawa.md)
: An input file representing some plasma interactions. Demonstrates the setup of a
  `plasma` phase, reactions parameterized using the
  [`two-temperature-plasma`](sec-yaml-two-temperature-plasma) and
  [`electron-collision-plasma`](sec-yaml-electron-collision-plasma) rate types, and the
  use of multiple `reactions` sections.

```{toctree}
:maxdepth: 2
:hidden:

h2o2
diamond
lithium_ion_battery
covdepsurf
oxygen-plasma-itikawa
```
