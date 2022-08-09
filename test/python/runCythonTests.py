"""
Unit tests for Cantera's Cython-based Python module.

This script is used by SCons to run the Python test suite, where ``pytest`` is used
as the test runner. All arguments are handled by SCons configurations.

This script gathers all Python tests, runs them, and prints a report::

    python runCythonTests.py

Extra command line arguments can be used to run subsets of the test suite, for example
all tests from ``test_thermo.py`` and ``test_kinetics.py``::

    python runCythonTests.py thermo kinetics

As an alternative, tests can be run using ``pytest`` directly, as illustrated for the
following examples (run from Cantera's root folder):

all tests::

    pytest test/python

all tests from ``test_transport.py``::

    pytest test/python/test_transport.py

all tests from the ``test_reactor.TestIdealGasReactor`` class::

    pytest test/python/test_reactor.py::TestIdealGasReactor

a single test::

    pytest test/python/test_onedim.py::TestDiffusionFlame::test_mixture_averaged
"""

import sys
import os
from pathlib import Path

CANTERA_ROOT = Path(__file__).parents[2]

try:
    import pytest
except ImportError:
        print("\n* ERROR: The Cantera Python test suite requires "
            "the Python package 'pytest'.")
        print("* ERROR: Use pip or conda to install 'pytest', "
            "which will enable this feature.")
        sys.exit(21)  # test/SConscript has special handling for this error code

import cantera

if __name__ == "__main__":
    print("\n* INFO: using Cantera module found at this location:")
    print(f"*       '{cantera.__file__}'")
    print(f"* INFO: Cantera version: {cantera.__version__}")
    print(f"* INFO: Git commit: {cantera.__git_commit__}\n")
    sys.stdout.flush()

    subset_start = 1
    fast_fail = False
    show_long = False
    verbose = False
    coverage = False
    if "fast_fail" in sys.argv:
        fast_fail = True
        subset_start += 1
    if "show_long" in sys.argv:
        show_long = True
        subset_start += 1
    if "verbose" in sys.argv:
        verbose = True
        subset_start += 1
    if "coverage" in sys.argv:
        coverage = True
        subset_start += 1

    base = CANTERA_ROOT / "test" / "python"
    subsets = []
    for name in sys.argv[subset_start:]:
        subsets.append(str(base.joinpath(f"test_{name}.py")))

    if not subsets:
        subsets.append(str(base))

    pytest_args = ["-raP", "--junitxml=pytest.xml"]
    if show_long:
        pytest_args += ["--durations=50"]
    if fast_fail:
        pytest_args.insert(0, "-x")
    if coverage:
        pytest_args.extend([
            "--cov=cantera",
            "--cov-config=test/python/coverage.ini",
            "--cov-report=xml:build/pycov.xml"
        ])
    if verbose:
        pytest_args.insert(0, "-v")
    else:
        pytest_args.append("--log-level=ERROR")

    ret_code = pytest.main(pytest_args + subsets)
    sys.exit(ret_code)
