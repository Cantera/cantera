import numpy as np

import cantera as ct
from . import utilities
from .utilities import unittest

from cantera._cantera import _py_to_any_to_py

class TestPyToAnyValue(utilities.CanteraTest):

    def check_conversion(self, value):
        out = _py_to_any_to_py(value)
        self.assertEqual(out, value)

    def check_inexact_conversion(self, value):
        out = _py_to_any_to_py(value)
        if isinstance(value, np.ndarray):
            self.assertEqual(out, value.tolist())
        else:
            self.assertEqual(out, list(value))

    def check_raises(self, value, ee, regex):
        with self.assertRaisesRegex(ee, regex):
            _py_to_any_to_py(value)

    def test_none(self):
        # None is converted to []
        out = _py_to_any_to_py(None)
        self.assertEqual(out, [])

    def test_set(self):
        self.check_raises({'a', 'b'}, NotImplementedError, "Python sets")

    def test_empty_list(self):
        self.check_conversion([])

    def test_empty_ndarray(self):
        self.check_inexact_conversion(np.ndarray((0,)))

    def test_empty_dict(self):
        self.check_conversion({})

    def test_scalar_string(self):
        self.check_conversion('spam')

    def test_scalar_int(self):
        self.check_conversion(3)

    def test_scalar_float(self):
        self.check_conversion(3.1415)

    def test_scalar_bool(self):
        self.check_conversion(True)

    def test_list_string(self):
        self.check_conversion(['spam', 'eggs'])

    def test_list_int(self):
        self.check_conversion([1, 2, 3])

    def test_list_float(self):
        self.check_conversion([1., 2., 3.])

    def test_list_bool(self):
        self.check_conversion([True, False])

    def test_list_various(self):
        self.check_conversion([True, 'spam', 3, 4.])

    def test_tuple(self):
        self.check_inexact_conversion((True, 'spam', 3, 4.))

    def test_ndarray1(self):
        self.check_inexact_conversion(np.random.randn(10))

    def test_ndarray2(self):
        self.check_inexact_conversion(np.random.randn(3, 2))

    def test_ndarray3(self):
        self.check_raises(np.random.randn(3, 2, 4),
            NotImplementedError, 'Cannot process float')

    def test_nested_string(self):
        self.check_conversion([['spam', 'eggs'], ['foo', 'bar']])

    def test_nested_int(self):
        self.check_conversion([[1, 2, 3], [4, 5, 6]])

    def test_nested_float(self):
        self.check_conversion([[1., 2., 3.], [4., 5., 6.]])

    def test_nested_bool(self):
        self.check_conversion([[True, False], [False, True]])

    def test_raises_string(self):
        self.check_raises([[['spam', 'eggs'], ['foo', 'bar']]],
            NotImplementedError, 'Cannot process string')

    def test_raises_int(self):
        self.check_raises([[[1, 2, 3], [4, 5, 6]]],
            NotImplementedError, 'Cannot process integer')

    def test_raises_float(self):
        self.check_raises([[[1., 2., 3.], [4., 5., 6.]]],
            NotImplementedError, 'Cannot process float')

    def test_raises_bool(self):
        self.check_raises([[[True, False], [False, True]]],
            NotImplementedError, 'Cannot process boolean')

    def test_inhomogeneous(self):
        self.check_raises([[1, 2], ['a', 'b']], NotImplementedError, 'inhomogeneous')

    def test_ragged(self):
        self.check_raises([[1, 2, 3], [4]], NotImplementedError, 'ragged')

    def test_dict(self):
        self.check_conversion({'a': 1, 'b': 2., 'c': 'eggs', 'd': True})

    def test_nested_dict(self):
        self.check_conversion({'a': 1, 'b': 2., 'c': {'d': 'eggs'}})
