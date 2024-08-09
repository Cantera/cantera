import pytest

pytest.register_assert_rewrite("pint.testing")

def pytest_addoption(parser):
    # this option is meant for debugging usage only
    parser.addoption(
        "--diagnose", action="store_true",
        default=False,
        help="Run only diagnostic tests, debugging purposes only.")


def pytest_configure(config):
    config.addinivalue_line("markers",
                            "diagnose: mark to only run this test.")


def pytest_collection_modifyitems(config, items):
    # run only diagnostic tests
    if config.getoption("--diagnose"):
        reason = "running only diagnostic tests."
        skip_diagnose = pytest.mark.skip(reason=reason)
        for item in items:
            if "diagnose" not in item.keywords:
                item.add_marker(skip_diagnose)
        return
