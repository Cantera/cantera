# Develop

(sec-compiling)=
## Compiling Cantera from Source

If you're interested in contributing new features to Cantera, or you want to try the
latest and development version, you will need to compile Cantera from source. The first
step is to make sure you have all the [](compiling/compilation-reqs) installed. Then,
you can [download the Cantera source code](compiling/source-code). Finally, you can
determine the appropriate configuration options and [compile
Cantera](compiling/configure-build) on your computer.

The following additional references may also be useful:

- [](compiling/dependencies.md)
- [](compiling/config-options)
- [](compiling/special-cases)

```{toctree}
:caption: Compiling Cantera from Source
:hidden:
:maxdepth: 1

compiling/compilation-reqs
compiling/source-code
compiling/configure-build
compiling/dependencies
compiling/config-options
compiling/special-cases
```

## How Cantera Works

```{caution}
This section is a work in progress.
```

- [](reactor-integration)

```{toctree}
:caption: How Cantera Works
:hidden:
:maxdepth: 1

reactor-integration
```

## Adding New Features to Cantera

- [](CONTRIBUTING)
- [](style-guidelines)
- [](vscode-tips)
- [](writing-tests)
- [](running-tests)
- [](writing-examples)
- [](doc-formatting)
- [](continuous-integration)

```{toctree}
:caption: Adding New Features to Cantera
:hidden:

CONTRIBUTING
style-guidelines
vscode-tips
writing-tests
running-tests
writing-examples
doc-formatting
```

## Generated API Code

- [](sourcegen)
- [](sourcegen-config)
- [](clib-extensions)
- [](dotnet-extensions)
- [](yaml-extensions)

```{toctree}
:caption: Generated API Code
:hidden:

sourcegen
sourcegen-config
clib-extensions
dotnet-extensions
yaml-extensions
```

(sec-distributing)=
## Releasing & Distributing Cantera

- [](distribution-packages/release-checklist)
- [](distribution-packages/pypi-sdist-wheel)
- [](distribution-packages/conda)
- [](distribution-packages/ubuntu-ppa)
- [](distribution-packages/windows-and-macos.md)

```{toctree}
:caption: Releasing & Distributing Cantera
:hidden:
:maxdepth: 1

distribution-packages/release-checklist
distribution-packages/pypi-sdist-wheel
distribution-packages/conda
distribution-packages/ubuntu-ppa
distribution-packages/windows-and-macos
```
