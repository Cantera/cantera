# Generated .NET API

```{caution}
The generated .NET API is an experimental part of Cantera and may be changed or
removed without notice.
```

The .NET API is written in C# and supports .NET 8 (and newer) on all platforms
that support both .NET and the Cantera C++ library. The .NET interface requires an
installation of the .NET 8.0 SDK to build.

The .NET API implementation draws on two parts:

1. [](sec-sourcegen-dotnet-generation), which translates CLib headers into
   corresponding C# implementations.

1. The .NET API included in `interfaces/dotnet`, which contains hand-coded API wrappers
   that extend generated code to form the actual user interface.

(sec-sourcegen-dotnet-install)=
## Building the .NET Interface

After [building the main Cantera library](sec-compiling) with default options, switch
to the `interfaces/dotnet` directory and run

```bash
dotnet build
```

The .NET test suite is invoked by running

```bash
dotnet test
```

In order to force re-import of generated code from sourcegen, run

```bash
dotnet clean
```

(sec-sourcegen-dotnet-generation)=
### C# Code Generation

Source generation for the .NET interface is fully integrated into the build process.
C# files used by the .NET API can be generated for informational purposes by running the
following command from the root folder of the Cantera source code:

```bash
sourcegen --api=csharp --output=build/csharp
```

Generated C# files are placed in the output folder `build/csharp`. Note that this
step requires installation of sourcegen via
`python -m pip install -e interfaces/sourcegen`.

## .NET Source Generator Overview

The primary C# project is `Cantera.csproj`. This project uses P/Invoke extensively with
the native Cantera library via the Cantera C interface (CLib), and wraps the low-level
interfaces with classes and concepts familiar to a .NET developer. As part of the build
process, it invokes [sourcegen](sourcegen) to scaffold the interop code and some of the
code for the wrapper objects, such as simple properties which can mapped directly to
CLib getter and setter functions. `Cantera.csproj` targets .NET 8. This project will be
released as a NuGet package.

`Cantera.Tests.csproj` contains the unit tests for the Cantera .NET library and targets
.Net 8.

The examples directory contains separate projects for each Cantera example. These will
be packaged as a template for use with the .NET CLI, allowing developers to create a
solution with all the examples using `dotnet new`.

## Extending the .NET API

Like all of Cantera, we welcome [contributions](CONTRIBUTING) to the .NET interface.
Contributors should review [general](sec-style-general) and
[C#-specific](sec-style-csharp) coding standards. Please note that contributing to the
.NET interface requires a relatively advanced understanding of C#.
This includes being comfortable with `unsafe` blocks and interop concepts such as raw
pointers, understanding how to use `Span<T>` and similar memory APIs, and having a solid
grip on object-oriented library design.

Discussions about contributing to the .NET interface assumes the developer has installed
the .NET SDK and is using the text editor of their choice, optionally with a plugin for
C# development installed. Visual Studio is _not_ required, and although contributors
may choose to use it, submissions should not depend on Visual Studio or require
the use of Windows.
