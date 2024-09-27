import cantera
from os import environ
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


def pytest_configure(config):
    config.addinivalue_line("markers", "slow_test: mark test as slow")

def pytest_collection_modifyitems(config, items):
    if environ.get("CT_SKIP_SLOW", "0") == "1":
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
def load_yaml():
    """
    Fixture to load YAML files safely.
    """
    def _load(yml_file):
        try:
            yaml_parser = yaml.YAML(typ="safe")
            with open(yml_file, "rt", encoding="utf-8") as stream:
                return yaml_parser.load(stream)
        except yaml.constructor.ConstructorError:
            # Ensure that the loader remains backward-compatible with legacy
            # ruamel.yaml versions (prior to 0.17.0).
            with open(yml_file, "rt", encoding="utf-8") as stream:
                return yaml.safe_load(stream)
    return _load


@pytest.fixture(scope="session", autouse=True)
def cantera_setup():
    """
    Fixture to set up Cantera environment for the entire test session.
    """
    # Add data directories
    cantera.add_directory(TEST_DATA_PATH)
    cantera.add_directory(CANTERA_DATA_PATH)
    cantera.print_stack_trace_on_segfault()
    cantera.CanteraError.set_stack_trace_depth(20)

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

    cantera.make_deprecation_warnings_fatal()
    cantera.add_directory(work_path)

    # Assign to the test class
    request.cls.test_work_path = work_path
    request.cls.using_tempfile = using_tempfile
    request.cls.test_data_path = TEST_DATA_PATH
    request.cls.cantera_data_path = CANTERA_DATA_PATH

    yield

    # Teardown: Remove the working directory if it was a temporary one
    if using_tempfile:
        try:
            for f in work_path.glob("*.*"):
                f.unlink()
            work_path.rmdir()
        except FileNotFoundError:
            pass
