(sec-python-type-annotations)=
# Python Type Annotations

Cantera employs type annotations in the Python interface primarily to support its use
within third-party Python scripts and libraries, and thus aims for full type coverage
of the external interface. These annotations indicate the intended types for all inputs
and outputs without restricting the runtime behavior of the code.

A secondary benefit of type annotations is to facilitate static analysis of Cantera's
Python code in order to catch potential typing-related errors using tools like
[mypy](https://mypy.readthedocs.io/en/latest/index.html),
[pyright](https://microsoft.github.io/pyright/), [ty](https://docs.astral.sh/ty/), and
[pyrefly](https://pyrefly.org/). These are only able to analyze pure Python code, so
they are of limited utility due to Cantera's extensive use of Cython syntax.

## Adding Type Annotations

Annotations should be added directly to all pure Python code (`.py` files). For Cython
code (`.pyx` files), annotation support is limited and optional; instead, type stubs
(`.pyi` files) must be added and maintained to document the external interface. Any
changes to the externally-facing API must include explicit type annotations. Type
annotations of implementation details is optional but encouraged for pure Python code,
and may be necessary in some situations to ensure the static type checks pass.

Ensure the syntax and features support the lowest-support Python version (currently
3.12), referencing the [Python](https://typing.python.org/en/latest/) and
[mypy](https://mypy.readthedocs.io/en/latest/index.html) documentation for guidance and
current best practices. Some additional recommendations which are more specific to
Cantera include:

* Type aliases should generally be prepended with an underscore.
* Aliases and special functions such as parametric `TypeGuard`s which will be used in
  many parts of the code should be placed in `_types.py`.
* Aliases representing any object which can be coerced into an object of type `{}`
  should be named `{}Like`, in keeping with `ArrayLike` from `NumPy`.
  * Examples: `ArrayLike`, `CompositionLike`, `_Func1Like`
* Prefer false positives to false negatives. In other words, if a type is too
  restrictive for a valid, if uncommon, usage, it should be switched for a more broad
  alternative.
  * Example: The `_TransportModel` alias could be a set of string literals which
    enumerate all built-in transport models, but this would incorrectly flag any custom
    extensions. It is therefore generalized to a `str` type.
  * Exception: For input collections with a strict format it is often better to
    restrict the type beyond what might be accepted at runtime. For example:
      * A `TypedDict` provides better documentation of a dictionary of options than
        something like `dict[str, str | float | dict[str, str | float]]`.
      * For the thermodynamic state setters, `tuple[float, float, CompositionLike]` is
        more informative than the alternative `Sequence[float | CompositionLike]`.
      * Clever use of `Protocol`s may also help to achieve optimal typing.
* The following should be avoided when possible:
  * Use of `Any` should be rare, and is typically reserved for dictionaries with
    complex or flexible contents (such as `kwargs`).
  * Use of `type: ignore` directives should be considered temporary workarounds and
    targeted for removal in future pull requests.
  * In the context of a pull request it is good practice to insert these as a
    placeholder for items where the developer requires assistance, allowing the CI
    job to progress and catch any other issues. Developers are encouraged to add a
    code comment elaborating on why it was deemed necessary.

## Verifying Type Annotations

The `type-checking` CI job runs a set of checks on the type annotations, and developers
should ensure these checks also pass locally prior to pushing. In this section we will
discuss each command, its purpose, and how to run it locally.

### Local Setup

Due to the use of compiled Cython code, most checks are performed on the built library
rather than running directly against the source code. It is thus recommended to add the
build directory to the environment variables following [](sec-using-build-dir):

::::{tab-set}
:::{tab-item} Linux
:sync: linux
```sh
export LD_LIBRARY_PATH=$(pwd)/build/lib
export PYTHONPATH=$(pwd)/build/python
```
:::

:::{tab-item} macOS
:sync: macos
```sh
export DYLD_LIBRARY_PATH=$(pwd)/build/lib
export PYTHONPATH=$(pwd)/build/python
```
:::

:::{tab-item} Windows (PowerShell)
:sync: windows
```pwsh
$Env:PYTHONPATH = (Get-Location).Path + '\build\python'
```
:::
::::

Currently Cantera primarily uses the `mypy` type checker with an additional check from
`pyright`. They can be installed with

```sh
conda install mypy>=1.19.0 lxml pyright
```

or

```sh
pip install mypy[reports]>=1.19.0 pyright
```

Some of the optional external dependencies provide type hints which are also employed
by Cantera when available, and thus must be installed for the type checks to pass:

```sh
conda install graphviz pandas pandas-stubs pint typing_extensions
```

### Static Type Correctness

This is the standard type-checking pass for a Python library which scans the raw
source code for type-related errors such as missing annotations, attempts to access
nonexistant methods/attributes, or empty collections (such as lists or dictionaries)
whose type could not be inferred.

This check can be run without first building Cantera by navigating to the
`interfaces/cython` directory and executing:

```sh
mypy -p cantera
```

Which will analyze the entire package and print the results to `stdout`. You can also
analyze specific files with `mypy [filename]`. Mypy should report no issues and will
print relevant information if any are found.

When not inside the cython directory, the `-p` option will attempt to check the built
library instead. Note that the configuration settings for `mypy` are included within
the `pyproject.toml` file in the `interfaces/cython` folder, which from the
root directory can be passed in on the command line like so:

```sh
mypy --config-file interfaces/cython/pyproject.toml -p cantera
```

### External Interface Type Coverage

These checks examine the public interface looking for missing or underspecified
annotations.

The first check uses Pyright's `--verifytypes` feature to print a quick summary of the
type coverage. Ideally there will be 100% coverage with known types; however, many
types originate from external libraries which may not themselves be fully defined.
These are omitted from the results by adding the `--ignoreexternal` option:

```sh
pyright --ignoreexternal --verifytypes cantera
```

Note that this check will additionally print various static type errors prior to the
coverage report. This is expected due to differences between `mypy` and `pyright`, and
it is not necessary to address these `pyright`-specific errors at this time.

The second check is a more thorough coverage report, similar in format to the unit test
coverage reports, generated using Mypy. This report highlights warnings and errors on a
per-line basis. While the CI generates a Cobertura-format XML file, the HTML format is
more straightforward for developers to view:

```sh
mypy --config-file interfaces/cython/pyproject.toml --html-report type_check/ -p cantera
```

This will generate a static HTML site with the index at `type_check/index.html`. As
with the static type check, this can also be run from `interfaces/cython` without first
building the code.

### Runtime Type Stub Correctness

The final category of type checks utilizes Mypy's `stubtest` utility to import Cantera
and compare the runtime code against the type annotations provided by the type stubs.
Because the Cython code cannot be analyzed statically, this provides a means to verify
the stubs are both correctly implemented and cover the public interface.

Because it performs runtime analysis, this can only be run on the built library. The
basic command is:

```sh
stubtest --mypy-config-file interfaces/cython/pyproject.toml \
--ignore-disjoint-bases --allowlist interfaces/cython/.mypyignore \
--concise cantera
```

Where the `--concise` option may be omitted to make `stubtest` print more details for
any identified issues. The `--allowlist` option points to a file containing a manually
curated list of methods to skip analyzing due to known and presently-unfixable errors;
for example, Cython methods have signatures which cannot be inspected at runtime, with
input parameters showing up as `(*args, **kwargs)`.

If `stubtest` reports any errors and it's not clear how to fix them, the allow list
can be automatically updated with the following commands:

```sh
# Remove everything from the allowlist except the pattern matches:
sed -i '16,$d' interfaces/cython/.mypyignore

# Generate a new allowlist and append it to the file:
stubtest --mypy-config-file interfaces/cython/pyproject.toml \
--ignore-disjoint-bases --allowlist interfaces/cython/.mypyignore \
--generate-allowlist cantera >> interfaces/cython/.mypyignore
```

This will ensure a consistent ordering. At present it is acceptable to add all stub
errors to the list to be addressed in later pull requests. Note that once a stub error
is fixed or otherwise rendered unnecessary, `stubtest` will raise an error until it is
removed from the allowlist.
