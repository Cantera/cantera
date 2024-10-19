import cantera
import os
import sys
from pathlib import Path
import pytest
import tempfile

try:
    from ruamel import yaml
except ImportError:
    import ruamel_yaml as yaml

pytest.register_assert_rewrite("pint.testing")

TEST_DATA_PATH = Path(__file__).parents[1] / "data"
CANTERA_DATA_PATH = Path(cantera.__file__).parent / "data"


def pytest_addoption(parser):
    parser.addoption(
        "--save-reference", action="store", default=None,
        help="Save the reference output files for specific tests. "
             "Options: diffusion, counterflow_premixed, counterflow_premixed_nonideal, "
             "combustor, wall"
    )

def pytest_configure(config):
    config.addinivalue_line("markers", "slow_test: mark test as slow")

def pytest_collection_modifyitems(config, items):
    if os.environ.get("CT_SKIP_SLOW", "0") == "1":
        skip_slow = pytest.mark.skip(reason="slow test")
        for item in items:
            if "slow_test" in item.keywords:
                item.add_marker(skip_slow)

@pytest.fixture
def allow_deprecated():
    cantera.suppress_deprecation_warnings()
    yield
    cantera.make_deprecation_warnings_fatal()


@pytest.fixture
def has_temperature_derivative_warnings():
    with pytest.warns(UserWarning, match="ddTScaledFromStruct"):
        # test warning raised for BlowersMasel and TwoTempPlasma derivatives
        yield

@pytest.fixture(scope="session")
def test_data_path():
    return TEST_DATA_PATH

@pytest.fixture(scope="session")
def cantera_data_path():
    return CANTERA_DATA_PATH

@pytest.fixture(scope="module", autouse=True)
def cantera_setup():
    """
    Fixture to set up Cantera environment for the entire test session.
    """
    # Add data directories
    cantera.add_directory(TEST_DATA_PATH)
    cantera.add_directory(CANTERA_DATA_PATH)
    cantera.print_stack_trace_on_segfault()
    cantera.CanteraError.set_stack_trace_depth(20)
    cantera.make_deprecation_warnings_fatal()

    # Yield control to tests
    yield


@pytest.fixture(scope="class", autouse=True)
def test_work_path(request, cantera_setup):
    """
    Fixture to create a working directory for a test class.
    This will only run for class-based tests.

    The check on the request.cls attribute is to ensure that this fixture
    does not run for function-based tests. There is likely a better way to
    do this, perhaps by disabling autouse and using the fixture explicitly
    in the test classes that need it.
    """

    if request.cls is None:
        # If this is not a class-based test, do nothing.
        yield
        return

    root_dir = Path(__file__).parents[2].resolve()
    if (root_dir / "SConstruct").is_file():
        work_path = root_dir / "test" / "work" / "python"
        using_tempfile = False
        try:
            work_path.mkdir(exist_ok=True)
        except FileNotFoundError:
            work_path = Path(tempfile.mkdtemp())
            using_tempfile = True
    else:
        work_path = Path(tempfile.mkdtemp())
        using_tempfile = True

    cantera.add_directory(work_path)

    # Assign to the test class
    request.cls.test_work_path = work_path
    request.cls.using_tempfile = using_tempfile

    yield

    # Teardown: Remove the working directory if it was a temporary one
    if using_tempfile:
        try:
            for f in work_path.glob("*.*"):
                f.unlink()
            work_path.rmdir()
        except FileNotFoundError:
            pass
