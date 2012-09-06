import cantera as ct
from . import utilities

class TestPureFluid(utilities.CanteraTest):
    def setUp(self):
        self.water = ct.Water()

    def test_critProperties(self):
        self.assertNear(self.water.critPressure, 22.089e6)
        self.assertNear(self.water.critTemperature, 647.286)
        self.assertNear(self.water.critDensity, 317.0)

    def test_setState(self):
        self.water._setState_Psat(101325, 0.5)
        self.assertNear(self.water.pressure, 101325)
        self.assertNear(self.water.vaporFraction, 0.5)

        self.water._setState_Tsat(500, 0.8)
        self.assertNear(self.water.temperature, 500)
        self.assertNear(self.water.vaporFraction, 0.8)
