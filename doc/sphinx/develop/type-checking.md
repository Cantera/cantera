(sec-python-type-annotations)=
# Python Type Annotations

Cantera employs type annotations in the Python interface primarily to support its use
within third-party Python scripts and libraries, and thus aims for full type coverage
of the external interface. These annotations indicate the intended types for all inputs
and outputs without restricting the runtime behavior of the code.

A secondary benefit of type annotations is to facilitate static analysis of Python code
in order to catch potential dynamic typing-related errors using tools like
[mypy](https://mypy.readthedocs.io/en/latest/index.html),
[pyright](https://microsoft.github.io/pyright/), [ty](https://docs.astral.sh/ty/), and
[pyrefly](https://pyrefly.org/). These are only able to analyze pure Python code, so
they are of limited utility due to Cantera's extensive use of Cython syntax.

## Adding Type Annotations

Annotations should be added directly to all pure Python code (`.py` files). For Cython
code (`.pyx` files), annotation support is limited and optional; instead, type stubs
(`.pyi` files) must be added and maintained to document the external interface.

Any changes to the externally-facing API must include explicit type annotations. Types
should be as narrow as feasible to fit the desired usage. For example, when an input
should be a tuple of two float values, prefer using `tuple[float, float]` over
`Iterable[float]` even if the method would technically work with any iterable. Usage of
`Any` types and `type: ignore` directives require an appropriate justification. Type
annotations of implementation details is optional but encouraged for pure Python code,
and may be necessary in some situations to ensure the static type checks pass.

Ensure all type annotations to be added are compatible with the lowest-supported Python
version (currently 3.12).

## Verifying Type Annotations

The `type-checking` CI job runs a set of checks on the type annotations, and developers
should ensure these checks also pass locally prior to pushing. In this section we will
discuss each command, its purpose, and how to run it locally.

### Local Setup

Currently Cantera primarily uses the `mypy` type checker with an additional check from
`pyright`, both of which can be installed via `pip` or `conda`. Note that
`mypy>=1.19.0` is required, and additional reporting formats can be generated if the
optional dependencies are installed with `mypy[reports]`.

Some of the optional external dependencies provide type hints which are also employed
by Cantera when available. It is therefore recommended to install them when developing
within the Python interface:

```yaml
dependencies:
- graphviz
- pandas
- pandas-stubs
- pint
- typing_extensions
```

Due to the use of compiled Cython code, most checks are performed on the built library
rather than running directly against the source code. It is thus recommended to follow
[](sec-using-build-dir).

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
