import pytest

# Loads fixtures defined in utilities.py for use in tests
pytest_plugins = ["utilities"]

pytest.register_assert_rewrite("pint.testing")
