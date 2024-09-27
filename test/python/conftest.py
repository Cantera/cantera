import sys
import os
import pytest


# Loads fixtures defined in utilities.py for use in tests
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
pytest_plugins = ["utilities"]

pytest.register_assert_rewrite("pint.testing")
