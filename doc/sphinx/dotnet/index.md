# .NET Interface Documentation

```{caution}
The .NET API is an experimental part of Cantera and may be changed or removed without
notice.
```

The .NET Interface is in preview and still missing many features needed for parity with
the C++ and Python interfaces. The interface is written in C# and relies on
[code generation](../develop/dotnet-extensions) based on [Cantera CLib](../clib/index).

## Installation and Examples

To install, follow steps outlined in the
[development reference](sec-sourcegen-dotnet-install).
Examples can be run with

```bash
dotnet run --project examples/X
```

where X is the name of the directory containing the example project, for example

```bash
dotnet run --project examples/SoundSpeed`
```

```{note}
For current information on the state of the .NET API, refer to GitHub issue
[Cantera/enhancements#156](https://github.com/Cantera/enhancements/issues/156).
```

```{toctree}
:hidden:
:maxdepth: 1
```
