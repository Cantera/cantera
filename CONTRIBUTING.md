# Contributing to Cantera

* For significant changes, please start a discussion on the Cantera
  Users' Group or create an issue on the [Cantera/enhancements](https://github.com/Cantera/enhancements/issues/new/choose) repository
  on GitHub to plan your modifications so that they can be implemented
  efficiently and in a way that doesn't conflict with any other planned
  future development
* Fork the `Cantera/cantera` repository on Github
* Clone your new repository or add it as a remote to an existing repository
* Check out the existing `main` branch, then start a new feature branch for
  your work
* When making changes, write code that is consistent with the surrounding code
  (see the [style guidelines](#style-guidelines) below)
* Add tests for any new features that you are implementing to either the
  GoogleTest-based test suite or the Python test suite.
* Add examples that highlight new capabilities, or update existing
  examples to make use of new features.
* As you make changes, commit them to your feature branch
  * Configure Git with your name and e-mail address before making any commits
  * Use descriptive commit messages (summary line of no more than 72 characters,
    followed by a blank line and a more detailed summary, if any)
  * Make related changes in a single commit, and unrelated changes in separate
    commits
  * Make sure that your commits do not include any undesired files, such as files
    produced as part of the build process or other temporary files.
  * Use Git's history-rewriting features (such as `git rebase -i`; see
    https://help.github.com/articles/about-git-rebase/) to organize your commits
    and squash "fixup" commits and reversions.
  * Do not merge your branch with `main`. If needed, you should rebase your branch
    onto the most recent `HEAD` commit of `main`.
  * Periodically run the test suite (`scons test`) to make sure that your
    changes are not causing any test failures.
* Push the changes on your new feature branch to your forked copy of the
  `Cantera/cantera` repository on GitHub.

* Submit a Pull Request on Github, from your forked copy. Check the results
  of the continuous-integration tests run using GitHub Actions and resolve
  any issues that arise.
* Additional discussion of good Git & Github workflow is provided at
  http://matplotlib.org/devel/gitwash/development_workflow.html and
  https://docs.scipy.org/doc/numpy-1.15.0/dev/gitwash/development_workflow.html
* Cantera is licensed under a [BSD
  license](https://github.com/Cantera/cantera/blob/main/License.txt) which
  allows others to freely modify the code, and if your Pull Request is accepted,
  then that code will be release under this license as well. The copyright for
  Cantera is held collectively by the contributors. If you have made a
  significant contribution, please add your name to the `AUTHORS` file.

# Style Guidelines

* Try to follow the style of surrounding code, and use variable names that
  follow existing patterns. Pay attention to indentation and spacing.
* Configure your editor to use 4 spaces per indentation level, and **never to
  use tabs**.
* Avoid introducing trailing whitespace
* Limit line lengths to 88 characters when possible
* Write comments to explain non-obvious operations
* Use whitespaces to improve code readability. Examples: after commas; before and
  after binary operators (`&&`/`||`/...), and comparisons (`<`/`>`/`==`/...); before and
  after equality signs `=` unless used for the assignment of a default parameter. For
  mathematical operators (`+`/`-`/`*`/`/` except `^`), whitespace should be added around
  the operators with the lowest priority (examples: `x + y + z`, `x*2 - 1`, or
  `(a+b) * (a-b)`). For additional guidance, refer to
  [Python PEP-8](https://peps.python.org/pep-0008/#whitespace-in-expressions-and-statements),
  where recommendations can be extrapolated to other programming languages
* Do not go out of your way to change formatting in otherwise unmodified code
* Write 'for example', 'such as', or 'that is' instead of using the Latin
  abbreviations 'i.e.' and 'e.g.'.

## C++

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
  with BibTeX-style entries added to `cantera.bib`. Equations can be added using
  LaTeX input bracketed by `@f[` and `@f]`. In-line math expressions are enclosed by
  a pair of `@f$` directives, for example `@f$ \sin(x) @f$`.
* Indicate the version added for new functions and classes with an annotation like
  `@since New in %Cantera X.Y` where `X.Y` is the next Cantera version. This notation
  should also be used to indicate significant changes in behavior.

### Style Guide

* Avoid defining non-trivial functions in header files
* Header files should include an 'include guard'
* Protected and private member variable names are generally prefixed with
  `m_`. For most classes, member variables should not be public.
* Class names use `InitialCapsNames`
* Methods use `camelCaseNames`
* Do not indent the contents of namespaces
* Code should follow the C++17 standard, with minimum required compiler versions
  GCC 7.0, Clang 4.0, MSVC 14.14 (Visual Studio 2017 version 15.7) and Intel 19.0.
* Cantera moves frequently used C++ standard namespace types and functions into the
  declarative region, meaning that the `std` scope resolution can be omitted. This
  applies to the following: `string`, `vector`, `map`, `set`, `pair`, `shared_ptr`,
  `make_shared`, `unique_ptr`, `make_unique` and `function`. Example: use `string`
  instead of `std::string`; a `using namespace std;` declaration is not required.
* Avoid manual memory management (that is, `new` and `delete`), preferring to use
  standard library containers, as well as `unique_ptr` and `shared_ptr` when dynamic
  allocation is required.
* Portions of Boost which are "header only" may be used. If possible, include
  Boost header files only within .cpp files rather than other header files to
  avoid unnecessary increases in compilation time. Boost should not be added
  to the public interface unless its existence and use is optional. This keeps
  the number of dependencies low for users of Cantera. In these cases,
  `CANTERA_API_NO_BOOST` should be used to conditionally remove Boost dependencies.
* While Cantera does not specifically follow these rules, the following style
  guides are useful references for possible style choices and the rationales behind them.
  * The Google C++ Style Guide: https://google.github.io/styleguide/cppguide.html
  * http://geosoft.no/development/cppstyle.html
* For any new code, do *not* use the `doublereal` and `integer` typedefs for the
  basic types `double` and `int`, but also do not go out of your way to change
  uses of these in otherwise unmodified code.
* Initialize member variables with their declarations, when possible, rather than using
  constructor-based initialization.

## Python

### Sphinx comments

* Cantera Python documentation is based on the Python documentation generator
  [Sphinx](https://www.sphinx-doc.org/en/master/index.html)
* All classes, member variables, and methods should include
  [Python docstrings](https://peps.python.org/pep-0257/#what-is-a-docstring)
* Docstrings should use annotations compatible with
  [automatic documentation generation from code](https://www.sphinx-doc.org/en/master/tutorial/automatic-doc-generation.html).
  For guidance, refer to existing Cantera documentation or online tutorials (see
  [example](https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html))
* Indicate the version added for new functions and classes with an annotation like
  `.. versionadded:: X.Y` where `X.Y` is the next Cantera version. Significant changes
  in behavior should be indicated with `.. versionchanged:: X.Y`.

### Style Guide

* Style generally follows PEP8 (https://www.python.org/dev/peps/pep-0008/)
* Code in `.py` and `.pyx` files needs to be written to work with Python 3
* The minimum Python version that Cantera supports is Python 3.8, so code should only
  use features added in Python 3.8 or earlier
* Please use double quotes in all new Python code

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
  (pointers, the `fixed` statement, etc). Pointers are used for the high-performance
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
