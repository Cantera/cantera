import numpy as np

import cantera as ct
from . import utilities
import copy


class TestPlasma(utilities.CanteraTest):
    def setUp(self):
        self.gas = ct.PlasmaPhase(infile='oxygen_plasma.yaml')
        self.gas.TPX = 1000, ct.one_atm, 'O2:1.0'
        self.gas.electric_field = 1e5
        self.gas.set_electron_energy_grid(np.linspace(0, 9.95, 200))
        self.gas.set_initial_mean_electron_energy(0.5)

    def test_electron_properties(self):
        self.assertNear(self.gas.electron_temperature, 13113, 1e-3)
        self.assertNear(self.gas.electron_mobility, 0.3985, 1e-3)
        self.assertNear(self.gas.electron_diffusivity, 0.5279, 1e-3)
        self.assertNear(self.gas.plasma_process_rate_coefficient(5), 1.55e-16, rtol=1e-5, atol=1e-17)
        self.assertNear(self.gas.plasma_process_reverse_rate_coefficient(5), 4.16e-16, rtol=1e-5, atol=1e-17)
        self.assertNear(self.gas.electron_total_collision_frequency, 3.433e11, 1e-3)
        self.assertNear(self.gas.electron_power_gain, 3.9811e9, 1e-3)
        self.assertNear(self.gas.electron_elastic_power_loss, 2.4114e7, 1e-3)
        self.assertNear(self.gas.electron_total_power_loss, 3.88e9, 1e-3)
        self.assertNear(self.gas.mean_electron_energy, 1.6949, 1e-3)
        self.assertNear(self.gas.electric_field, 1e5, 1e-3)

    def test_change_electric_field_freq(self):
        Te0 = self.gas.electron_temperature
        self.gas.set_reuse_EEDF(True)
        self.gas.electric_field_freq = 1e9
        Te = self.gas.electron_temperature
        self.gas.electric_field_freq = 2e9
        self.assertLess(Te, Te0)
        self.assertLess(self.gas.electron_temperature, Te)
        self.assertNear(Te, Te0, 1e-3)
        self.assertNear(self.gas.electron_temperature, Te, 1e-3)
        self.gas.electric_field_freq = 0.0
        self.gas.set_reuse_EEDF(False)

    def test_change_gas_temperature(self):
        Te0 = self.gas.electron_temperature
        self.gas.TP = 1100, ct.one_atm * 1.1
        # The gas temperature is important only when E/N is small
        self.assertNotEqual(self.gas.electron_temperature, Te0)
        self.assertNear(self.gas.electron_temperature, Te0, 1e-4)
        self.gas.electric_field = 1.0
        Te = self.gas.electron_temperature
        self.gas.TP = 1000, ct.one_atm
        self.assertLess(0.05, abs(self.gas.electron_temperature - Te) / Te)
        # set conditions back
        self.gas.electric_field = 1e5
        self.gas.TPX = 1000, ct.one_atm, 'O2:1.0'

    def test_change_electric_field_strength(self):
        grid = np.linspace(0.0, 2.0, num=1000)
        self.gas.set_electron_energy_grid(grid)
        self.gas.electric_field = 1.0
        # The electron temperature approach gas temperature when E is small
        self.assertNear(self.gas.electron_temperature, self.gas.T, 1e-3)
        # set conditions back
        self.gas.electric_field = 1e5
