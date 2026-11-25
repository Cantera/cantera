(sec-lxcat2yaml)=

# LXCat to YAML conversion

The `lxcat2yaml` utility converts electron-collision cross-section data from
LXCat XML files to Cantera's YAML format. Both the legacy LXCat XML format and
LXCat XML version 1.1 are supported.

The converter writes collision data in the root-level `electron-collisions`
section. Each entry includes the collision kind and tabulated cross sections:

- A **collision-only** entry specifies `target` and `kind`. It is used by the
  electron energy distribution function (EEDF) solver but does not create a
  chemical reaction.
- A **reactive** entry specifies `equation`, `type: electron-collision-plasma`,
  and `kind`. Its cross-section data are included directly in the same entry;
  no separate collision reference is needed.
- Elastic and effective collisions are always written as collision-only
  entries, even when inserting data into a mechanism.

Every entry also includes `product`, preserving the original LXCat heavy-product
label independently of any species mapping. When no heavy product is specified,
`product: none` is written as a literal string (not YAML `null`). If a process
contains several heavy products, their LXCat labels are joined with ` + `.
Cross-section-specific XML comments are copied to the optional `note` field.

:::{versionchanged} 4.0
Converted cross sections are written directly in the root-level
`electron-collisions` section. Reactive entries include an `equation`
and can participate in both kinetic reaction calculations and the
electron Boltzmann solver. Collision-only entries do not create
kinetic reactions.
:::

## Basic conversion

To convert an LXCat XML file without inserting the results into an existing
mechanism, run:

```bash
lxcat2yaml --input=my-cross-sections.xml --output=my-cross-sections.yaml
```

If the `lxcat2yaml` executable is not available on the system path, the
converter can instead be invoked as a Python module:

```bash
python -m cantera.lxcat2yaml \
    --input=my-cross-sections.xml \
    --output=my-cross-sections.yaml
```

If `--output` is omitted, the input filename extension is replaced with
`.yaml`.

**Without `--insert`, all converted processes are collision-only entries.**
No `reactions` section or reaction equations are generated, even if the LXCat
file provides product species. For example:

```yaml
electron-collisions:
- target: Ar
  kind: effective
  product: Ar
  note: EFFECTIVE Momentum transfer CROSS SECTION.
  energy-levels: [0.0, 0.001, 0.002]
  cross-sections: [7.5e-20, 7.5e-20, 7.1e-20]
- target: Ar
  kind: ionization
  product: Ar^+
  threshold: 15.8
  energy-levels: [15.8, 16.0, 17.0]
  cross-sections: [0.0, 2.02e-22, 1.34e-21]
```

The generated `electron-collisions` section can be used in a plasma
mechanism independently of chemical reactions.

For example, a plasma phase with `reactions: none` can still use
these cross sections to calculate the EEDF. To include reactive
entries in chemical kinetics, the phase must use a kinetics model
and select the corresponding reaction entries.

## Inserting data into a mechanism

Use `--insert` to add the converted collision data to an existing Cantera
mechanism:

```bash
lxcat2yaml \
    --input=my-cross-sections.xml \
    --database=Phelps \
    --mech=plasma-mechanism.yaml \
    --phase=plasma \
    --insert \
    --output=plasma-mechanism-with-cross-sections.yaml
```

The `--mech` option is required with `--insert`. When the input mechanism has
multiple phases, `--phase` selects the phase whose species are used for
filtering and reaction generation. The `--database` option selects the `id`
attribute of an LXCat `Database` element; if omitted, all databases are
considered.

When a mechanism is supplied:

- Processes whose target species is absent from the selected phase are skipped.
- Processes with an available target are retained as collision data, even if
  their product species are absent from the phase.
- With `--insert`, an inelastic process becomes a kinetic reaction only if
  all its products can be matched to species in the selected phase and the
  electron products are consistent with its collision kind.
- Elastic and effective processes remain collision-only entries.
- The converter enables the generated electron-collision reactions in the
  selected phase's reaction selection when needed.

When invoked from the command line with `--insert`, the utility also tries
to load the generated mechanism to validate the selected phase.

For example, if `Ar` and `Ar+` are available in the phase, an imported
ionization can be written as:

```yaml
electron-collisions:
- equation: Ar + e => Ar+ + e + e
  type: electron-collision-plasma
  kind: ionization
  product: Ar^+
  threshold: 15.8
  energy-levels: [15.8, 16.0, 17.0]
  cross-sections: [0.0, 2.02e-22, 1.34e-21]
```

If `Ar+` is not present, the same process is retained as a collision-only
entry with `target: Ar` and `kind: ionization`, without an `equation`.

:::{note}
When using `--insert`, the output mechanism is serialized from the selected
Cantera phase. Other phases from the original mechanism are not guaranteed to
be included in the output file.
:::

## Mapping LXCat products to mechanism species

LXCat product labels may differ from the species names used in a Cantera
mechanism. Use the optional `--mapping` argument to supply a YAML mapping:

```yaml
products:
  "Ar*(11.5eV)": Ar_excited
  "N2(A3Sigma)": "N2(A)"
```

Then run:

```bash
lxcat2yaml \
    --input=my-cross-sections.xml \
    --mech=plasma-mechanism.yaml \
    --phase=plasma \
    --insert \
    --mapping=products.yaml \
    --output=plasma-mechanism-with-cross-sections.yaml
```

Only the `products` mapping is supported. It changes the species names used
in generated reaction equations, but never changes the original LXCat label
stored in `product`. The mapped species must already be declared in the
selected Cantera phase; the converter does not create new species.

## Cross-section and threshold conventions

LXCat energies and cross sections are expected in eV and m², respectively.
Units are checked when specified in the XML; if omitted, these units are
assumed. Energy levels must be strictly increasing, and energies and cross
sections must be finite and non-negative.

An explicitly specified LXCat threshold is written as `threshold`. When no
threshold is provided, this field is omitted so that Cantera can apply its
normal threshold-inference rules.

For **excitation and ionization**, an explicit threshold causes points at or
below the threshold to be replaced by a point with zero cross section at the
threshold; higher-energy points are preserved. For **elastic, effective, and
attachment** collisions, the tabulated cross sections are retained without
this threshold normalization.

A cross-section table containing only one point can be extended to two
points: for elastic, effective, and attachment processes, a point at 0 eV
is added with the same cross section; for excitation and ionization, the
added point has zero cross section at the explicit threshold, or at 0 eV if
no threshold is supplied. The measured energy must be strictly above the
inserted point. Otherwise, the one-point table is rejected.

The converter does not resample cross sections onto the EEDF energy grid or
add high-energy extrapolation points. Interpolation, extrapolation, and rate
integration are handled by Cantera when the mechanism is loaded. In
particular, the default extrapolation policy is not overridden by the
converter.

Common electron labels (`e`, `e-`, and `Electron`) and LXCat charge notation,
such as `Ar^+` and `O^-`, are normalized for reaction generation. The `product` field retains the LXCat
spelling (for example, `Ar^+`), while an equation uses the Cantera species name
(for example, `Ar+`). Product names still need to match the species in the
selected phase, possibly through a mapping.

## Command-line options

```{eval-rst}
.. argparse::
   :module: cantera.lxcat2yaml
   :func: create_argparser
   :prog: lxcat2yaml
```
