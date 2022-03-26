"""
Unit tests for Cantera's Cython-based Python module.

This script gathers all the tests defined 'cantera.test' module, runs them,
and prints a report. Extra command line arguments can be used to run subsets
of the test suite, for example:

all tests from 'test_thermo.py' and 'test_kinetics.py':

    python runCythonTests.py thermo kinetics

all tests from the 'test_reactor.TesTIdealGasReactor' class:

    python runCythonTests.py reactor.TestIdealGasReactor

a single test:

    python runCythonTests.py onedim.TestDiffusionFlame.test_mixture_averaged
"""

import sys
import os
from pathlib import Path

CANTERA_ROOT = Path(__file__).parents[2]
os.chdir(str(CANTERA_ROOT / "test" / "work"))

try:
    import pytest
except ImportError:
        print("\n* ERROR: The Cantera Python test suite requires "
            "the Python package 'pytest'.")
        print("* ERROR: Use pip or conda to install 'pytest', "
            "which will enable this feature.")
        sys.exit(21)  # test/SConscript has special handling for this error code

import cantera
import cantera.test

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
    if "fast_fail" in sys.argv:
        fast_fail = True
        subset_start += 1
    if "show_long" in sys.argv:
        show_long = True
        subset_start += 1
    if "verbose" in sys.argv:
        verbose = True
        subset_start += 1

    base = Path(cantera.__file__).parent.joinpath("test")
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
    if verbose:
        pytest_args.insert(0, "-v")
    else:
        pytest_args.append("--log-level=ERROR")

    ret_code = pytest.main(pytest_args + subsets)
    sys.exit(ret_code)
