# Adding Features to the .NET API

```{warning}
The .NET API is an experimental part of Cantera and
may be changed or removed without notice.
```

The (original) version of the .NET API builds on the traditional CLib API.

## Background

```{note}
Information below was created for the initial version of **sourcegen**, which was
specific to the .NET interface. While the information is still accurate, future versions
of the `csharp` option will be based on the doxygen XML tree.
```

The `sourcegen.generate_source` function crudely parses the CLib header files and generates intermediate objects which represent the functions:
* header file path
* funcs: list of `Func` objects containing
    * return type (`string`)
    * name (`string`)
    * params: list of function arguments (`ArgList`)
    * optional annotations

`sourcegen.generate_source` then delegates the source generation to a language-specific sub-package. The sub-package creates the source files in a location appropriate to the build for that language’s build process.

`sourcegen.generate_source` takes two mandatory arguments:
1. The name of the language-specific module.
1. The full path of the directory where generated files shall be placed.

The language-specific sub-packages export a class that derives from `SourceGenerator`:
1. An initializer which takes the output directory and the parsed config object\
    The initializer stores these for later use.
1. A `generate_source` function which takes the a list of header files with their parsed functions\
    The generate_source function uses this information to generate syntactically correct source code in the destination language

The implementation details for the source generator class is of course highly dependent on the build process of the destination language.

### The config file

Each sub-package can contain a yaml-based config file named `config.yaml`. The core script recognizes two special keys:
1. `ignore_files`: a list of header file names\
    These files will be ignored entirely from source generation, for example because they cannot be parsed directly or contain functionality that is not planned for implementation in the destination language.
1. `ignore_funcs`: a mapping of header file names to lists of function names\
    The listed functions contained within those files will not be scaffolded, for example because they cannot be translated automatically and need to be written by hand in the destination language.

The config file may contain additional values for use by the language-specific sub-package.
