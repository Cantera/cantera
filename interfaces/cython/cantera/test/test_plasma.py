import cantera as ct
from . import utilities
import numpy as np

class TestPlasmaPhase(utilities.CanteraTest):
    def setUp(self):
        self.phase = ct.Plasma('oxygen-plasma.yaml', transport_model=None)

    def test_set_get_electron_energy_grid(self):
        grid = np.linspace(0.01, 10, num=9)
        self.phase.electron_energy_grid = grid
        self.assertArrayNear(grid, self.phase.electron_energy_grid)
