# Julia Interface Documentation

```{caution}
The Julia API is an experimental part of Cantera and may be changed or removed without
notice.
```

The Julia interface is in preview and still missing many features needed for parity
with the C++ and Python interfaces. The interface is written in Julia and relies on
[code generation](../develop/julia-extensions) based on [Cantera CLib](../clib/index).
It calls `libcantera` through Julia's built-in `ccall`.

## Installation and Examples

To install, follow the steps outlined in the
[development reference](sec-sourcegen-julia-install). Examples can be run with

```bash
julia --project=interfaces/julia interfaces/julia/examples/X.jl
```

where X is the name of an example, for example

```bash
julia --project=interfaces/julia interfaces/julia/examples/basic.jl
```

A gallery of the same examples is available under [](/examples/julia/index).

## API Reference

The docstrings in `interfaces/julia/src` document every type and function exported by
`Cantera.jl`. They are also available from the Julia REPL help mode by typing `?`
followed by a name.

```{include} api-reference.md
```

```{note}
For current information on the state of the Julia API, refer to
GitHub issues and enhancement requests, specifically:

- [**Open Pull Requests**](https://github.com/Cantera/cantera/pulls?q=is%3Apr+state%3Aopen+label%3Ajulia)
- [**Open Issues**](https://github.com/Cantera/cantera/issues?q=is%3Aissue+state%3Aopen+label%3Ajulia)
- [**Open Enhancements**](https://github.com/Cantera/enhancements/issues?q=is%3Aissue+state%3Aopen+label%3Ajulia)
```
