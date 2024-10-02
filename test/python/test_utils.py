import numpy as np
from pathlib import Path
import pytest
from pytest import approx

import cantera as ct
from cantera._utils import _py_to_any_to_py, _py_to_anymap_to_py


class TestUnitSystem:

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
            assert dim in checks
            assert checks[dim] == unit

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
            assert dim in checks
            assert checks[dim] == unit

    def test_activation_energy(self):
        system = ct.UnitSystem({"activation-energy": "eV"})
        units = system.units
        assert units["activation-energy"] == "eV"

        system = ct.UnitSystem({"activation-energy": "K"})
        units = system.units
        assert units["activation-energy"] == "K"

    def test_raises(self):
        with pytest.raises(ct.CanteraError, match="non-unity conversion factor"):
            ct.UnitSystem({"temperature": "2 K"})
        with pytest.raises(ct.CanteraError, match="non-unity conversion factor"):
            ct.UnitSystem({"current": "2 A"})

    def test_convert_to_default(self):
        system = ct.UnitSystem()
        assert system.convert_to("3 cm", "m") == 0.03
        assert system.convert_to(4, "mm") == 4000.0
        assert system.convert_to("3 cm", ct.Units("m")) == 0.03
        assert system.convert_to(4, ct.Units("m")) == 4

    def test_convert_activation_energy(self):
        system = ct.UnitSystem()
        assert system.convert_activation_energy_to("3 J/mol", "J/kmol") == 3000
        assert system.convert_activation_energy_to(4, "J/mol") == 0.004

    def test_convert_rate_coeff(self):
        system = ct.UnitSystem({"length": "cm"})
        assert system.convert_rate_coeff_to(11, ct.Units("m^3/kmol/s")) == approx(11e-6)
        assert system.convert_rate_coeff_to("22 m^3/mol/s", "m^3/kmol/s") == approx(22e3)

    def test_convert_to_custom(self):
        system = ct.UnitSystem({"length": "cm", "mass": "g"})
        assert system.convert_to(10000, "m^2") == 1.0
        assert system.convert_to(500, "kg") == 0.5

    def test_convert_to_array(self):
        system = ct.UnitSystem({"length": "km"})
        x = np.array(((3, 4), (0.5, 2.0), (1.0, 0.0)))
        assert system.convert_to(x, "m") == approx(1000 * x)

    def test_convert_activation_energy_to_array(self):
        system = ct.UnitSystem({"activation-energy": "J/mol"})
        x = np.array(((3, 4), (0.5, 2.0), (1.0, 0.0)))
        assert system.convert_activation_energy_to(x, "J/kmol") == approx(1000 * x)

    def test_convert_rate_coeff_to_array(self):
        system = ct.UnitSystem({"length": "cm"})
        x = np.array(((3, 4), (0.5, 2.0), (1.0, 0.0)))
        assert system.convert_rate_coeff_to(x, "m^2/kmol/s") == approx(0.0001 * x)

    def test_convert_to_sequence(self):
        system = ct.UnitSystem({"length": "km"})
        x = [("3000 mm", 4), (0.5, 2.0), 1.0]
        x_m = system.convert_to(x, "m")
        assert x_m[0][0] == 3.0
        assert x_m[1][1] == 2000.0
        assert x_m[2] == 1000.0

    def test_convert_activation_energy_to_sequence(self):
        system = ct.UnitSystem({"activation-energy": "J/mol"})
        x = [("3000 K", 4), (0.5, 2.0), 1.0]
        x_m = system.convert_activation_energy_to(x, "J/kmol")
        assert x_m[0][0] == approx(3000 * ct.gas_constant)
        assert x_m[1][1] == 2000.0
        assert x_m[2] == 1000.0

    def test_convert_rate_coeff_to_sequence(self):
        system = ct.UnitSystem({"length": "cm"})
        x = [("3000 mm^3/kmol/s", 4), (0.5, 2.0), 1.0]
        x_m = system.convert_rate_coeff_to(x, ct.Units("m^3/kmol/s"))
        assert x_m[0][0] == approx(3000 / 1e9)
        assert x_m[1][1] == approx(2 / 1e6)
        assert x_m[2] == approx(1e-6)

    def test_convert_errors(self):
        system = ct.UnitSystem()
        with pytest.raises(ct.CanteraError):
            system.convert_to("eggs", "J/kmol")

        with pytest.raises(TypeError):
            system.convert_to("5 cm", True)

        with pytest.raises(TypeError):
            system.convert_to(None, "J/kmol")

        with pytest.raises(TypeError):
            system.convert_to(5, 6)

        with pytest.raises(ct.CanteraError):
            system.convert_activation_energy_to(4, "m^3/s")

        with pytest.raises(TypeError):
            system.convert_activation_energy_to(5, True)

        with pytest.raises(ct.CanteraError):
            system.convert_activation_energy_to("spam spam eggs", "K")

        with pytest.raises(TypeError):
            system.convert_activation_energy_to({"spam": 5}, "K")

        with pytest.raises(TypeError):
            system.convert_rate_coeff_to(11, ["m^3"])

        with pytest.raises(TypeError):
            system.convert_rate_coeff_to({"spam": 13}, "m^6/kmol^2/s")

class TestPyToAnyValue:

    def check_conversion(self, value, check_type=None):
        out, held_type = _py_to_any_to_py(value)
        assert out == value
        if check_type is not None:
            assert held_type == check_type

    def check_inexact_conversion(self, value, check_type=None):
        out, held_type = _py_to_any_to_py(value)
        if isinstance(value, np.ndarray):
            assert out == value.tolist()
        else:
            assert out == list(value)
        if check_type is not None:
            assert held_type == check_type

    def check_raises(self, value, ee, regex):
        with pytest.raises(ee, match=regex):
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


@pytest.fixture(scope='class')
def setup_any_map(request):
    request.cls.data = _py_to_anymap_to_py(request.cls.orig_data)

@pytest.mark.usefixtures("setup_any_map")
class TestAnyMap:
    orig_data = {
            "units": {"length": "mm", "energy": "kJ"},
            "group1": {
                "a": 5000,
                "b": "12 MJ",
                "c": "8000 K",
                "d": [16, "10 cm^2"]
            },
            "group2": {
                "units": {"mass": "g"},
                "x": 1300
            }
        }

    def test_units_simple(self):
        assert self.data['group1'].convert('a', 'm') == 5.0
        assert self.data['group1'].convert('b', 'J') == 12e6
        assert self.data['group1'].convert('d', 'm^2') == [16e-6, 10e-4]

    def test_units_activation_energy(self):
        assert self.data['group1'].convert_activation_energy('a', 'J/kmol') == 5e6
        assert (self.data['group1'].convert_activation_energy('c', 'J/kmol')
                == approx(8000 * ct.gas_constant))

    def test_units_nested(self):
        assert self.data['group2'].convert('x', 'J/kg') == 1300 * 1e6

    def test_set_quantity(self):
        params = ct.AnyMap()
        params.set_quantity('spam', [2, 3, 4], 'Gg')
        params.set_quantity('eggs', 10, ct.Units('kg/m^3'))
        params.set_activation_energy('beans', 5, 'K')

        converted = _py_to_anymap_to_py(params)
        assert converted['spam'] == [2e6, 3e6, 4e6]
        assert converted['eggs'] == approx(10)
        assert converted['beans'] == approx(5 * ct.gas_constant)

        # Unit conversions are deferred
        outer = ct.AnyMap()
        outer['units'] = {'mass': 'g', 'activation-energy': 'K'}
        outer['inner'] = params
        converted = _py_to_anymap_to_py(outer)
        assert converted['inner']['spam'] == approx([2e9, 3e9, 4e9])
        assert converted['inner']['eggs'] == approx(10e3)
        assert converted['inner']['beans'] == approx(5)

        outer.set_quantity('cheese', {'gouda': 5.5}, 'J/kg')
        with pytest.raises(ct.CanteraError):
            _py_to_anymap_to_py(outer)

        outer.set_activation_energy('cheese', 12, 'V/m')
        with pytest.raises(ct.CanteraError):
            _py_to_anymap_to_py(outer)

class TestListDataFiles:

    def test_list_data_files(self):
        data_files = ct.list_data_files()
        assert "gri30.yaml" in data_files
        assert str(Path("example_data/oxygen-plasma-itikawa.yaml")) in data_files
