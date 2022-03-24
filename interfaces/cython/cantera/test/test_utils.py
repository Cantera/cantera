import numpy as np

import cantera as ct
from . import utilities

from cantera._cantera import _py_to_any_to_py


class TestUnitSystem(utilities.CanteraTest):

    def test_default(self):
        units = ct.UnitSystem().units
        checks = {
            "activation-energy": "J / kmol",
            "current": "A",
            "energy": "J",
            "length": "m",
            "mass": "kg",
            "pressure": "Pa",
            "quantity": "kmol",
            "temperature": "K",
            "time": "s",
        }
        for dim, unit in units.items():
            self.assertIn(dim, checks)
            self.assertEqual(checks[dim], unit)

    def test_cgs(self):
        system = ct.UnitSystem({
            "length": "cm", "mass": "g", "time": "s",
            "quantity": "mol", "pressure": "dyn / cm^2", "energy": "erg",
            "activation-energy": "cal / mol"})
        units = system.units
        checks = {
            "activation-energy": "cal / mol",
            "current": "A",
            "energy": "erg",
            "length": "cm",
            "mass": "g",
            "pressure": "dyn / cm^2",
            "quantity": "mol",
            "temperature": "K",
            "time": "s",
        }
        for dim, unit in units.items():
            self.assertIn(dim, checks)
            self.assertEqual(checks[dim], unit)

    def test_activation_energy(self):
        system = ct.UnitSystem({"activation-energy": "eV"})
        units = system.units
        self.assertEqual(units["activation-energy"], "eV")

        system = ct.UnitSystem({"activation-energy": "K"})
        units = system.units
        self.assertEqual(units["activation-energy"], "K")

    def test_raises(self):
        with self.assertRaisesRegex(ct.CanteraError, "non-unity conversion factor"):
            ct.UnitSystem({"temperature": "2 K"})
        with self.assertRaisesRegex(ct.CanteraError, "non-unity conversion factor"):
            ct.UnitSystem({"current": "2 A"})


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
        self.check_conversion(None, "void")

    def test_set(self):
        # Sets are converted to lists
        self.check_inexact_conversion({"a", "b"}, "vector<string>")

    def test_empty_list(self):
        self.check_conversion([])

    def test_empty_ndarray(self):
        self.check_inexact_conversion(np.ndarray((0,)))

    def test_empty_dict(self):
        self.check_conversion({})

    def test_scalar_string(self):
        self.check_conversion("spam", "string")

    def test_scalar_int(self):
        self.check_conversion(3, "long int")

    def test_scalar_float(self):
        self.check_conversion(3.1415, "double")

    def test_scalar_bool(self):
        self.check_conversion(True, "bool")

    def test_list_string(self):
        self.check_conversion(["spam", "eggs"], "vector<string>")

    def test_list_int(self):
        self.check_conversion([1, 2, 3], "vector<long int>")

    def test_list_float(self):
        self.check_conversion([1., 2., 3.], "vector<double>")

    def test_list_bool(self):
        self.check_conversion([True, False], "vector<bool>")

    def test_list_various(self):
        self.check_conversion([True, "spam", 3, 4., {"foo": "bar"}],
                              "vector<AnyValue>")

    def test_tuple(self):
        self.check_inexact_conversion((True, "spam", 3, 4.), "vector<AnyValue>")

    def test_ndarray1(self):
        self.check_inexact_conversion(np.random.randn(10), "vector<double>")

    def test_ndarray2(self):
        self.check_inexact_conversion(np.random.randn(3, 2), "vector<vector<double>>")

    def test_ndarray3(self):
        # Each inner AnyValue holds a vector<vector<double>>
        self.check_inexact_conversion(np.random.randn(3, 2, 4), "vector<AnyValue>")

    def test_nested_string(self):
        self.check_conversion([["spam", "eggs"], ["foo", "bar"]],
                              "vector<vector<string>>")

    def test_nested_int(self):
        self.check_conversion([[1, 2, 3], [4, 5, 6]], "vector<vector<long int>>")

    def test_nested_float(self):
        self.check_conversion([[1., 2., 3.], [4., 5., 6.]], "vector<vector<double>>")

    def test_nested_bool(self):
        self.check_conversion([[True, False], [False, True]], "vector<vector<bool>>")

    def test_multi_dict(self):
        vv = [{"a": [["spam", "eggs"], ["foo", "bar"]], "b": {"c": 4}}, {"d": 3}]
        self.check_conversion(vv, "vector<AnyMap>")

    def test_dict(self):
        self.check_conversion({"a": 1, "b": 2., "c": "eggs", "d": True}, "AnyMap")

    def test_nested_dict(self):
        self.check_conversion({"a": 1, "b": 2., "c": {"d": "eggs"}}, "AnyMap")

    def test_unconvertible(self):
        self.check_raises(object(), ct.CanteraError, "Unable to convert")

    def test_unconvertible2(self):
        self.check_raises([3+4j, 1-2j], ct.CanteraError, "Unable to convert")
