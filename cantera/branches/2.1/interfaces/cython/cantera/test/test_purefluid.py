import cantera as ct
from . import utilities

class TestPureFluid(utilities.CanteraTest):
    def setUp(self):
        self.water = ct.Water()

    def test_critical_properties(self):
        self.assertNear(self.water.critical_pressure, 22.089e6)
        self.assertNear(self.water.critical_temperature, 647.286)
        self.assertNear(self.water.critical_density, 317.0)

    def test_set_state(self):
        self.water.PX = 101325, 0.5
        self.assertNear(self.water.P, 101325)
        self.assertNear(self.water.X, 0.5)

        self.water.TX = 500, 0.8
        self.assertNear(self.water.T, 500)
        self.assertNear(self.water.X, 0.8)
