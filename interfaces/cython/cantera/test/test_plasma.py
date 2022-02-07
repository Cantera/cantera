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

    def test_isotropic_velocity_electron_energy_distribution(self):
        grid = np.linspace(0.01, 10, num=9)
        self.phase.electron_energy_grid = grid
        self.phase.Te = 2e5
        mean_electron_energy = 3.0 / 2.0 * (self.phase.Te * ct.gas_constant /
                               (ct.avogadro * ct.electron_charge))
        self.assertNear(mean_electron_energy , self.phase.mean_electron_energy)

    def test_user_specified_electron_energy_distribution(self):
        grid = np.linspace(0, 1, num=2)
        distrb = np.linspace(0, 1, num=2)
        self.phase.set_electron_energy_distribution(grid, distrb)
        self.assertArrayNear(grid, self.phase.electron_energy_grid)
        self.assertArrayNear(distrb, self.phase.electron_energy_distribution)
        self.assertNear(self.phase.mean_electron_energy, 0.2)
        electron_temp = 2.0 / 3.0 * (self.phase.mean_electron_energy *
                        ct.avogadro * ct.electron_charge / ct.gas_constant)
        self.assertNear(self.phase.Te, electron_temp)
