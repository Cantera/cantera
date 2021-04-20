import numpy as np

import cantera as ct
from . import utilities
from .utilities import unittest

from cantera._cantera import _py_to_any_to_py

class TestPyToAnyValue(utilities.CanteraTest):

    def check_conversion(self, value, check_type=None):
        out, held_type = _py_to_any_to_py(value)
        self.assertEqual(out, value)
        if check_type is not None:
            self.assertEqual(held_type, check_type)

    def check_inexact_conversion(self, value, check_type=None):
        out, held_type = _py_to_any_to_py(value)
        if isinstance(value, np.ndarray):
            self.assertEqual(out, value.tolist())
        else:
            self.assertEqual(out, list(value))
        if check_type is not None:
            self.assertEqual(held_type, check_type)

    def check_raises(self, value, ee, regex):
        with self.assertRaisesRegex(ee, regex):
            _py_to_any_to_py(value)

    def test_none(self):
        out, held_type = _py_to_any_to_py(None)
        self.assertEqual(out, None)
        self.assertEqual(held_type, "void")

    def test_set(self):
        # Sets are converted to lists
        self.check_inexact_conversion({'a', 'b'}, 'vector<string>')

    def test_empty_list(self):
        self.check_conversion([])

    def test_empty_ndarray(self):
        self.check_inexact_conversion(np.ndarray((0,)))

    def test_empty_dict(self):
        self.check_conversion({})

    def test_scalar_string(self):
        self.check_conversion('spam', 'string')

    def test_scalar_int(self):
        self.check_conversion(3, 'long int')

    def test_scalar_float(self):
        self.check_conversion(3.1415, 'double')

    def test_scalar_bool(self):
        self.check_conversion(True, 'bool')

    def test_list_string(self):
        self.check_conversion(['spam', 'eggs'], 'vector<string>')

    def test_list_int(self):
        self.check_conversion([1, 2, 3], 'vector<long int>')

    def test_list_float(self):
        self.check_conversion([1., 2., 3.], 'vector<double>')

    def test_list_bool(self):
        self.check_conversion([True, False], 'vector<bool>')

    def test_list_various(self):
        self.check_conversion([True, 'spam', 3, 4., {'foo': 'bar'}],
                              'vector<AnyValue>')

    def test_tuple(self):
        self.check_inexact_conversion((True, 'spam', 3, 4.), 'vector<AnyValue>')

    def test_ndarray1(self):
        self.check_inexact_conversion(np.random.randn(10), 'vector<double>')

    def test_ndarray2(self):
        self.check_inexact_conversion(np.random.randn(3, 2), 'vector<vector<double>>')

    def test_ndarray3(self):
        # Each inner AnyValue holds a vector<vector<double>>
        self.check_inexact_conversion(np.random.randn(3, 2, 4), 'vector<AnyValue>')

    def test_nested_string(self):
        self.check_conversion([['spam', 'eggs'], ['foo', 'bar']],
                              'vector<vector<string>>')

    def test_nested_int(self):
        self.check_conversion([[1, 2, 3], [4, 5, 6]], 'vector<vector<long int>>')

    def test_nested_float(self):
        self.check_conversion([[1., 2., 3.], [4., 5., 6.]], 'vector<vector<double>>')

    def test_nested_bool(self):
        self.check_conversion([[True, False], [False, True]], 'vector<vector<bool>>')

    def test_multi_dict(self):
        vv = [{'a': [['spam', 'eggs'], ['foo', 'bar']], 'b': {'c': 4}}, {'d': 3}]
        self.check_conversion(vv, 'vector<AnyMap>')

    def test_dict(self):
        self.check_conversion({'a': 1, 'b': 2., 'c': 'eggs', 'd': True}, 'AnyMap')

    def test_nested_dict(self):
        self.check_conversion({'a': 1, 'b': 2., 'c': {'d': 'eggs'}}, 'AnyMap')

    def test_unconvertible(self):
        class Foo: pass
        self.check_raises(Foo(), ct.CanteraError, "Unable to convert")

    def test_unconvertible2(self):
        self.check_raises([3+4j, 1-2j], ct.CanteraError, "Unable to convert")
