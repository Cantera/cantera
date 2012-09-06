import cantera as ct
from . import utilities

class TestOnedim(utilities.CanteraTest):
    def test_instantiate(self):
        gas = ct.Solution('h2o2.xml')

        flame = ct.FreeFlame(gas)

    def test_badInstantiate(self):
        solid = ct.Solution('diamond.xml', 'diamond')
        with self.assertRaises(TypeError):
            flame = ct.FreeFlame(solid)

    def test_instantiateSurface(self):
        gas = ct.Solution('diamond.xml', 'gas')
        solid = ct.Solution('diamond.xml', 'diamond')
        interface = ct.Solution('diamond.xml', 'diamond_100', (gas, solid))

        surface = ct.ReactingSurface1D()
        surface.setKinetics(interface)

