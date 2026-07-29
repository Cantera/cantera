# Generated Julia API

```{caution}
The generated Julia API is an experimental part of Cantera and may be changed or
removed without notice.
```

The Julia API is written in Julia and supports Julia 1.9 (and newer) on all platforms
that support both Julia and the Cantera C++ library. It calls `libcantera` directly
through the [CLib API](../clib/index) using Julia's built-in `ccall`, so it requires
neither a Julia compiler toolchain of its own nor the Cantera Python module.

The Julia API implementation draws on two parts:

1. [](sec-sourcegen-julia-generation), which translates the CLib specifications into
   corresponding `ccall` wrappers.

1. The Julia API included in `interfaces/julia`, which contains hand-coded API wrappers
   that extend generated code to form the actual user interface.

(sec-sourcegen-julia-install)=
## Building the Julia Interface

The generated bindings are scaffolded as part of the regular build, so
[building the main Cantera library](sec-compiling) with default options is sufficient:

```bash
scons build
```

Because the bindings are scaffolded from the Doxygen tag file that sourcegen parses
(see [](sourcegen)), [Doxygen](https://www.doxygen.org) must be installed even when no
documentation is being built.

The Julia package is not registered in Julia's General registry, so it is used from the
source tree by activating its project. Julia locates the shared library and the
mechanism files through two environment variables:

```bash
export CANTERA_LIBRARY_PATH=/path/to/cantera/lib   # or a conda environment
export CANTERA_DATA=/path/to/cantera/data          # for gri30.yaml and friends
```

Neither is required when Cantera was built in-tree: the interface falls back to the
sibling `build/lib` and `data` directories of the source checkout.

The dependencies are instantiated, and the test suite is invoked, by running

```bash
julia --project=interfaces/julia -e 'using Pkg; Pkg.instantiate()'
julia --project=interfaces/julia interfaces/julia/test/runtests.jl
```

The HTML API reference is built with [Documenter](https://documenter.juliadocs.org) by
running

```bash
scons julia_docs
```

which writes its output to `build/doc/html/julia` so that it is deployed alongside the
rest of the Cantera documentation. Building the reference loads the compiled library in
order to introspect docstrings, so it requires a successful `scons build` first.

(sec-sourcegen-julia-generation)=
### Julia Code Generation

Source generation for the Julia interface is fully integrated into the build process by
`interfaces/julia/SConscript`. Julia files used by the Julia API can be generated for
informational purposes by running the following command from the root folder of the
Cantera source code:

```bash
sourcegen --api=julia --output=build/julia
```

Generated Julia files are placed in the `src/generated` subfolder of the output folder,
that is, `build/julia/src/generated`. Note that this step requires installation of
sourcegen via `python -m pip install -e interfaces/sourcegen`.

The generator emits one file per CLib specification — `libctthermo.jl` for
`ctthermo.yaml`, and so on — plus a `_manifest.jl` listing them, which is what
`LibCantera.jl` iterates over to `include` the bindings. Because the file names follow
from the specifications, `SConscript` can declare them as build targets without running
the generator first.

Two aspects of the generator are worth noting for developers:

- The output directory is emptied on every run. A binding file for a specification that
  is no longer generated would otherwise linger and shadow the current CLib.
- Generated files are not committed to the repository, since they can be recreated on
  demand. A checkout that has never been built therefore has no `src/generated`
  directory, and `using Cantera` fails until `scons build` has run.

Like the other sourcegen sub-packages, the Julia generator is configured by
`interfaces/sourcegen/src/sourcegen/julia/config.yaml`. In addition to the standard
`ignore_files` and `ignore_funcs` keys described in [](sec-sourcegen-details), it
defines a `c_type_crosswalk` mapping each C type used by the CLib to the Julia type used
in the `ccall` signature, for example `double*` to `Ptr{Float64}`. Generation aborts
with an error naming the offending type if the CLib introduces a type that has no entry,
rather than emitting a binding that would be wrong at the ABI level.

## Julia Interface Overview

The public API is a conventional Julia package rooted at `interfaces/julia/src`, with
`Cantera.jl` as the module entry point. Its structure follows Cantera's C++ class
hierarchy rather than the flat CLib namespace:

`LibCantera.jl`
: Internal module: locates `libcantera` and includes the generated bindings.

`errors.jl`
: Translation of CLib error sentinels into `CanteraError` exceptions.

`handles.jl`
: Wrapper supertype, phase views, and string/array marshalling helpers.

`solution.jl`
: `Solution`, the primary user-facing object.

`thermo.jl`, `kinetics.jl`, `transport.jl`
: Property accessors for the three managers of a phase.

`reaction.jl`, `func1.jl`
: `Reaction` and `Func1` objects.

`reactor.jl`, `reactornet.jl`, `connectors.jl`
: Zero-dimensional reactor networks.

`onedim.jl`
: One-dimensional flames.

`multiphase.jl`, `rdiag.jl`, `solutionarray.jl`
: Multiphase equilibrium, reaction-path diagrams, and collections of states.

`utils.jl`
: Physical constants and library-level functions.

Public names are exported from the file that defines them, so adding a method requires
touching only one file.

## Implementation Notes

The following are the aspects of the interface that are least obvious from reading the
source, and the ones most likely to be disturbed by a well-intentioned change.

### Library discovery happens at run time, not precompile time

`LibCantera.libcantera` is a `Ref{String}` that is assigned in the module's `__init__`
function, and the generated wrappers dereference it as `libcantera[]` on every call.
This indirection is deliberate. Julia caches a precompiled image of the module, so a
library path resolved at precompile time would be frozen into that image, and a later
change to `CANTERA_LIBRARY_PATH` would be silently ignored. Resolving the path in
`__init__` keeps it correct for the session that actually runs.

Only the source tree's own `data` directory is registered by `Cantera.__init__`.
`CANTERA_DATA` and the conda `share/cantera/data` directory are already added by
Cantera's `Application::setDefaultDirectories` before any input file is read, and adding
them again from Julia would be redundant.

### Errors are sentinel values, not exceptions

The CLib cannot propagate C++ exceptions, so it signals failure through return values:
functions returning `int32_t` return a negative value, and functions returning `double`
return `DERR` (-999.999). Every call therefore passes through `check` (for integers) or
`checkd` (for doubles), which retrieve the message with `ct_getCanteraError` and throw a
`CanteraError`. A new hand-written wrapper that omits these will report failures as
implausible numbers rather than as errors.

The one complication is that a negative `int32_t` is not always an error: CLib functions
that look up a name return a negative `npos` value when the name is absent. Wrappers
such as `species_index` must therefore test for that case *before* calling `check`.

### Strings and arrays use a size-query protocol

CLib string getters are called twice: once with a null buffer to obtain the required
length, and again with a buffer of that size. `get_string` encapsulates this, taking a
closure so that any getter's leading arguments can be bound at the call site.

Array getters take a caller-allocated buffer, filled by `get_array!` or by `get_array`,
which allocates first. Most array properties expose both a normal and an in-place
`foo!(out, gas)` variant, the latter to avoid an allocation per call in loops.
Non-`Float64` input vectors are converted at the boundary by `as_f64`, since a `ccall`
needs a contiguous `Vector{Float64}`.

### Handle lifetime is the interface's main hazard

Every CLib object is an `Int32` handle into a C++-side cabinet, and Julia's garbage
collector knows nothing about the object it refers to. Each wrapper is therefore a
`mutable struct` with a finalizer that calls the matching `*_del` function through an
idempotent `close!`, guarded by a `closed` flag. These `close!` methods deliberately
swallow exceptions: they run on the finalizer path, where throwing is not useful.

Two consequences shape the type definitions:

- The `ThermoPhase`, `Kinetics` and `Transport` views do not own their handles — they
  borrow them from a `Solution`. Each is parametrized on the type of its owner and
  stores a reference to it, so that a view keeps the `Solution` alive even if the
  variable holding the `Solution` goes out of scope first. Without that reference, the
  `Solution` could be finalized while a view is still in use, leaving the view pointing
  at a deleted cabinet entry.
- Reactors clone the phase they are constructed from, as the C++ core expects. A
  reactor's state is consequently *not* observable through the `Solution` passed to
  its constructor — that object keeps its original state. The reactor's own phase is
  reached with `reactor_phase`, and scalar reactor state has direct accessors such as
  `temperature(reactor)`.

Not every phase has a kinetics or transport manager. `Solution` records a `-1` for a
manager that could not be created and raises a descriptive `CanteraError` when one is
requested, rather than passing an invalid handle across the boundary.

### Repetitive accessors are metaprogrammed

Families of accessors that differ only in the CLib function they call — the partial
molar properties, for instance — are generated by `@eval` over a list of name/function
pairs. Docstrings do not survive this transformation automatically, so the loops carry
an explicit `@doc` for each generated method; a generated accessor without one is
invisible in the API reference.

### One-dimensional flames replicate the Python solver's staging

`solve!` for a `FreeFlame` does not simply call `sim1D_solve`. It reproduces the staged
strategy of the Python interface's `Sim1D.solve(auto=True)`: solving on a sequence of
progressively finer uniform grids, re-seeding the initial guess at each stage,
attempting an energy-enabled solve and falling back to a fixed-temperature solve, then
running a refinement pass. It also reproduces the `DomainTooNarrow` check, doubling the
domain width when the temperature gradient at either edge is significant compared to the
average gradient. This staging is what makes ignition cases converge from a cold start,
and its absence is felt as solver failures rather than as wrong answers.

## Extending the Julia API

Like all of Cantera, we welcome [contributions](CONTRIBUTING) to the Julia interface.
Contributors should review the [general](sec-style-general) coding standards; Julia
code in Cantera also observes the 88-character line length limit used across the other
interfaces, and every hand-written source file carries the license preamble.

Adding a property to the interface usually means working from the bottom up:

1. If the underlying CLib function does not yet exist, add a recipe to the appropriate
   specification in `interfaces/sourcegen/src/sourcegen/headers`, as described in
   [](clib-extensions). The Julia binding then appears automatically on the next build.

1. Add the hand-written wrapper to the corresponding `interfaces/julia/src` file, using
   `check`/`checkd` and the marshalling helpers rather than calling `ccall` directly,
   and export it from that same file.

1. Add a test to the matching `interfaces/julia/test/test_*.jl` file. Tests are written
   against the Python interface's behavior; `interfaces/julia/test/reference_values.md`
   records where the expected values come from.

Because the low-level bindings are generated, a missing wrapper is far more often a
missing CLib recipe than a missing Julia method. Running
`sourcegen --api=julia --output=build/julia` and inspecting the generated file is the
quickest way to confirm whether a function crosses the boundary at all.
