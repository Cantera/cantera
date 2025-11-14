# Radiation Models and RadLib Integration

Cantera's one-dimensional flow domains support multiple radiation property models.
The built-in **TabularPlanckMean** model only requires tabulated absorption data,
while the optional [RadLib](https://github.com/BYUignite/radlib) integration adds
Planck-mean, WSGG, and RCSLW property calculators. This page explains how to enable
RadLib when building Cantera, how to configure the runtime inputs, and how to tune
RadLib-specific options from Python or C++.

## Enabling RadLib at build time

RadLib is an optional dependency. When running `scons`, choose one of the following
options:

```
# Use a system installation:
scons build system_radlib=y
# or supply include / library directories explicitly:
scons build system_radlib=default radlib_include=/opt/radlib/include \
    radlib_libdir=/opt/radlib/lib

# Use the bundled submodule (builds RadLib locally under build/ext):
git submodule update --init ext/radlib
scons build system_radlib=n
```

The configuration script searches `radlib_include` / `radlib_libdir` when provided,
otherwise the standard compiler paths are used. Successful detection sets the
`CT_ENABLE_RADLIB` flag so that the RadLib adapters are compiled and linked into
Cantera. If RadLib cannot be found (or the option is left at the default) the build
falls back to the TabularPlanckMean model and radiation calls to `RadLib.*` will raise
an informative `CanteraError`.

## Runtime configuration via `radiation-parameters.yaml`

At runtime, the `Flow1D` constructor looks for a file named
`radiation-parameters.yaml`. The search order is the same as every other Cantera
input: the directory containing the parent file, the working directory, and then the
data directories on the Cantera search path (set via `CANTERA_DATA`). A warning is
issued the first time the file cannot be found or parsed, and the defaults are used.

The file can specify global radiation settings under a top-level `Radiation` key:

```yaml
Radiation:
  property: RadLib.RCSLW     # TabularPlanckMean, RadLib.PlanckMean, RadLib.WSGG, ...
  solver: OpticallyThin
  emissivity:
    left: 0.0
    right: 0.0
  radlib:
    fvsoot: 1.0e-6
    nGray: 25
    Tref: 1500.0
    Pref: 101325.0

PMAC:
  CO2:
    fit-type: polynomial
    data: [18.741, -121.31, 273.5, -194.05, 56.31, -5.8169]
```

The `PMAC` section is optional data used by the TabularPlanckMean calculator; only
species present in the current mechanism and listed in `PMAC` are extracted from the
solution vector.

The `Radiation.radlib` block provides default RadLib settings (soot fraction, number
of gray gases, and reference state) that are cached inside the domain. These defaults
are preserved whenever a RadLib model is selected at construction time or via
`set_radiation_models`.

## Adjusting RadLib options programmatically

In C++, the `Flow1D` interface exposes

```cpp
flow->setRadLibOptions(fvsoot, nGray, Tref, Pref);
double fv = flow->radlibFvSoot();
int nGray = flow->radlibNGray();
```

and similar getters for the reference temperature and pressure. In Python, the same
functionality is available through `FlowBase.set_radlib_options` and the
`FlowBase.radlib_options` property:

```python
gas = ct.Solution("h2o2.yaml")
flow = ct.FreeFlow(gas)
flow.set_radlib_options(fvsoot=5e-6, ngray=30, Tref=1500.0, Pref=2 * ct.one_atm)
flow.set_radiation_models("RadLib.RCSLW")
print(flow.radlib_options)
# {'fvsoot': 5e-06, 'nGray': 30, 'Tref': 1500.0, 'Pref': 202650.0}
```

These settings persist until they are changed again.

## Behavior when RadLib is not enabled

If Cantera was built without RadLib support, attempting to select a RadLib property
(`RadLib.PlanckMean`, `RadLib.WSGG`, or `RadLib.RCSLW`) or calling
`setRadLibOptions()` raises a `CanteraError` explaining that RadLib was not
enabled and how to rebuild with the required dependency. Existing Tabular and
optically-thin functionality continues to work unchanged.
