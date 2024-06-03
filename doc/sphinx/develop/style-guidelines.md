# Code Style Guidelines

The following style guidelines are recommended for all new code added to Cantera.
Following these guidelines will simplify the review process for pull requests and make
it easier for others to understand your code in the context of Cantera as a whole.

## General Style

* Try to follow the style of surrounding code, and use variable names that follow
  existing patterns. Pay attention to indentation and spacing.
* Configure your editor to use 4 spaces per indentation level, and **never to use
  tabs**.
* Avoid introducing trailing whitespace
* Limit line lengths to 88 characters when possible
* Write comments to explain non-obvious operations
* Use whitespace to improve code readability. Examples:
  * after commas
  * before and after binary operators (`&&`/`||`/...) and comparisons (`<`/`>`/`==`/...)
  * before and after equality signs `=` unless used for the assignment of a default
    parameter.
  * For mathematical operators (`+`/`-`/`*`/`/` except `^`), whitespace should be added
    around the operators with the lowest priority (examples: `x + y + z`, `x*2 - 1`, or
    `(a+b) * (a-b)`).
  * For additional guidance, refer to [Python
    PEP-8](https://peps.python.org/pep-0008/#whitespace-in-expressions-and-statements),
    where recommendations can be extrapolated to other programming languages
* Do not go out of your way to change formatting in otherwise unmodified code
* Write 'for example', 'such as', or 'that is' instead of using the Latin abbreviations
  'i.e.' and 'e.g.'.

## C++

* Avoid defining non-trivial functions in header files
* Header files should include an 'include guard'
* Protected and private member variable names are generally prefixed with `m_`. For most
  classes, member variables should not be public.
* Class names use `InitialCapsNames`
* Methods use `camelCaseNames`
* Do not indent the contents of namespaces
* Code should follow the C++17 standard, with minimum required compiler versions GCC
  7.0, Clang 4.0, MSVC 14.14 (Visual Studio 2017 version 15.7) and Intel 19.0.
* Cantera moves frequently used C++ standard namespace types and functions into the
  declarative region, meaning that the `std` scope resolution can be omitted. This
  applies to the following: `string`, `vector`, `map`, `set`, `pair`, `shared_ptr`,
  `make_shared`, `unique_ptr`, `make_unique` and `function`. Example: use `string`
  instead of `std::string`; a `using namespace std;` declaration is not required.
* Avoid manual memory management (that is, `new` and `delete`), preferring to use
  standard library containers, as well as `unique_ptr` and `shared_ptr` when dynamic
  allocation is required.
* Portions of Boost which are "header only" may be used. If possible, include Boost
  header files only within .cpp files rather than other header files to avoid
  unnecessary increases in compilation time. Boost should not be added to the public
  interface unless its existence and use is optional. This keeps the number of
  dependencies low for users of Cantera. In these cases, `CANTERA_API_NO_BOOST` should
  be used to conditionally remove Boost dependencies.
* While Cantera does not specifically follow these rules, the following style guides are
  useful references for possible style choices and the rationales behind them.
  * [The Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html)
  * [GeoSoft C++ Programming Style Guidelines](http://geosoft.no/development/cppstyle.html)
* For any new code, do *not* use the obsolete `doublereal` and `integer` typedefs for
  the basic types `double` and `int`, but also do not go out of your way to change uses
  of these in otherwise unmodified code.
* Initialize member variables with their declarations, when possible, rather than using
  constructor-based initialization.

### Doxygen Comments

* All classes, member variables, and methods should use
  [Doxygen-style comments](https://www.doxygen.nl/manual/docblocks.html).
* Comments should provide brief and/or detailed descriptions. For example, comment
  blocks starting with `/**` or `//!` use the autobrief feature (comments are split into
  brief and detailed descriptions at the first dot `'.'`). For short comments, the C++
  style `//!` is preferred; do not use `///` or `/*!` comment styles in new code.
* [Doxygen commands](https://www.doxygen.nl/manual/commands.html) should use the `@`
  prefix instead of `\` in order to better differentiate from LaTeX input.
* Whenever appropriate, classes and functions should be added to
  [Doxygen groupings](https://www.doxygen.nl/manual/grouping.html) using the `@ingroup`
  command. Alternatively, entire code sections can be added using the `@addtogroup`
  command, where grouped classes and functions are bracketed by `@{` and `@}`.
* If applicable, new features should reference literature using the `@cite` command,
  with BibTeX-style entries added to `cantera.bib`.
* Indicate the version added for new functions and classes with an annotation like
  `@since New in %Cantera X.Y` where `X.Y` is the next Cantera version. This notation
  should also be used to indicate significant changes in behavior.

## Python

### Style Guide

* Style generally follows [PEP8](https://www.python.org/dev/peps/pep-0008/)
* The minimum Python version that Cantera supports is Python 3.8, so code should only
  use features added in Python 3.8 or earlier
* Please use double quotes in all new Python code

### Sphinx comments

* Cantera Python documentation is based on the Python documentation generator
  [Sphinx](https://www.sphinx-doc.org/en/master/index.html)
* All classes, member variables, and methods should include
  [Python docstrings](https://peps.python.org/pep-0257/#what-is-a-docstring)
* New classes and global functions need to be added to one of the pages in
  `doc/sphinx/python` so they will appear in the API reference. For a Cython class
  defined in `.pyx` file, the argument list needs to be repeated as part of the
  `.. autoclass::` declaration; for a function or pure Python class, the signature is
  automatically read by Sphinx.
* Docstrings should use annotations compatible with
  [autodoc](https://www.sphinx-doc.org/en/master/tutorial/automatic-doc-generation.html).
  For guidance, refer to existing Cantera documentation or online tutorials (see
  [example](https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html))
* Indicate the version added for new functions and classes with an annotation like
  `.. versionadded:: X.Y` where `X.Y` is the next Cantera version. Significant changes
  in behavior should be indicated with `.. versionchanged:: X.Y`.
* To document an attribute of a Cython class, include a docstring *below* the member
  declaration (in the `.pxd` file). For example:
  ```cython
  cdef class ReactorBase:
      cdef public dict node_attr
      """ Attributes of a node """
  ```

## C#

* C# coding conventions should follow https://docs.microsoft.com/en-us/dotnet/csharp/fundamentals/coding-style/coding-conventions
* All identifiers should follow the naming conventions in the above document including
  * Prefixing with `_` for private instance fields (`_foo`, unlike C++)
  * Prefixing with `s_` for private static fields (`s_bar`), `t_` for private
    `[ThreadStatic]` fields (`t_baz`).
  * Initial caps names for class methods (`DoSomething()`, unlike C++)
* Give the opening brace of a statement block its own line (unlike C++), except empty
  blocks, which may be written as an `{ }` (for example, a constructor which calls
  a base-class constructor only).
* Use only one statement per line.
* Always use statement blocks (`{ ... }`) for the bodies of statements that can take
  either a statement block or a single statement (`if`, `for`, etc.)
* Use file-scoped namespaces in each new file.
* Do not take any extra Nuget dependencies in the `Cantera.csproj` project.
* Use C# XML Doc-Comments on types and members, including at least the `<summary>` tag.
  Always include a doc comment for types, but for members with self-explanatory names,
  you may omit the doc comment and suppress the build error that would be thrown with
  `#pragma warning disable/restore CS1591`.
  * C# doc-comments use `///`, unlike Cantera's preferred use of `//!` for C++
* Do not expose any code requiring the `unsafe` keyword via a public API
  (pointers, the `fixed` statement, etc.). Pointers are used for the high-performance
  interop layer with the native Cantera library, but such access should have a
  “safe” wrapper, such as a `Span<T>` or a managed array.
* Do not allow exceptions to pass uncaught out of a callback invoked from native code,
  as the interop layer cannot marshall exceptions between managed and native code,
  and the process will crash. Use `CallbackException.Register()` within a catch-all
  block to log the exception for later throwing back in managed code.
* The primary API for accessing Cantera is the `Application` class, which handles
  required static initialization of the library. When exposing a new wrapper for CLib
  functionality, do not expose a public constructor. Rather, mark the constructor
  `internal` and wrap it in an appropriate factory method in the `Application` class
  (`public static CreateFoo(string filename) { ... }`).
