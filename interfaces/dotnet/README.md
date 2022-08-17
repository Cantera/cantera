# Cantera – .NET Interface

This directory and its children contain the interface to Cantera for .NET languages
such as C# and F#. It is written in C# and supports .NET Standard 2.0
(for the primary project) and .NET 6 (and newer) on all platforms that support both
.NET and the Cantera C++ library.

### Project Layout

The primary C# project is Cantera.csproj. This project uses P/Invoke extensively with
the native Cantera library via the Cantera C interface (CLib), and wraps the low-level
interfaces with classes and concepts familiar to a .NET developer. As part of the build
process, it invokes [sourcegen](/interfaces/sourcegen)
to scaffold the interop code and some of the code for the wrapper objects, such as
simple properties which can mapped directly to CLib getter and setter functions.
Cantera.csproj targets .NET Standard 2.0 and .NET 6. This project will be released as
a NuGet package.

Cantera.Tests.csproj contains the unit tests for the Cantera .NET library and targets
.Net 6.

The examples directory contains separate projects for each Cantera example. These will
be packaged as a template for use with the .NET CLI, allowing developers to create a
solution with all the examples using `dotnet new`.

### Contributing

Like all of Cantera, we welcome contributions to the .NET interface. Contributors should
read the [code of conduct](/CODE_OF_CONDUCT.md) and check out our general and
C#-specific coding standards in [CONTRIBUTING.md](/CONTRIBUTING.md). Please note that
contributing to the .NET interface requires a relatively advanced understanding of C#.
This includes being comfortable with `unsafe` blocks and interop concepts such as
raw pointers, understanding how to use `Span<T>` and similar memory APIs, and having
a solid grip on object-oriented library design.

### Status

The .NET Interface is in preview and still missing many features
needed for parity with the C++ and Python interfaces. Interested contributors can see
the status of the project at [Enhancement #156](https://github.com/Cantera/enhancements/issues/156).
