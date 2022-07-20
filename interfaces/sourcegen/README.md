## sourcegen – Python based source generator for creating Cantera interfaces for other languages

A common script “parses” the CLib header files and generates intermediate objects which represent the functions:
* return type (string)
* name (string)
* params: list of
    * return type
    * name

The script delegates the source generation to a language-specific module. The module creates the source files in a location appropriate to the build for that languages build process.

The common script takes two mandatory arguments:
1. The name of the language-specific module.
1. The full path of the directory where generated files shall be placed.

The language-specific modules export a class named “SourceGenerator”:
1. A constructor which takes the output directory and the parsed config object\
    The constructor stores these for later use.
1. A “generate_source” function which will takes the included file and its list of parsed functions\
    The generate_source function generates source files for each file it receives as it is called
    and/or builds up data structures containing source symbols to be generated later.
1. A “finalize” function which takes no arguments
    The finalize function performs any final code generation tasks and necessary clean up.

The implementation details for the SourceGenerator classes is of course highly dependent on the build process of the destination language.

In addition, each module can contain a yaml-based config file named “config.yaml”. The core script recoginizes a key called “ignore”, which is a map. Each key in the “ignore” map is the name of a clib header file. The value is an array which is either empty or contains strings representing function names defined in the header file. The core script interprets the ignore map as follows:
* for each key followed by an empty array, the entire file is excluded from code generation
* for each key followed-by a non-empty array, only the functions referenced in the array are excluded from code generation.

The config file may contain additional values for use by the language-specific module.
