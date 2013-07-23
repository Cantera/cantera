from .utilities import unittest
import numpy as np

import cantera as ct
from . import utilities

class TestTransport(utilities.CanteraTest):
    def setUp(self):
        self.phase = ct.Solution('h2o2.xml')
        self.phase.X = [0.1, 1e-4, 1e-5, 0.2, 2e-4, 0.3, 1e-6, 5e-5, 0.4]
        self.phase.TP = 800, 2*ct.one_atm

    def test_scalar_properties(self):
        self.assertTrue(self.phase.viscosity > 0.0)
        self.assertTrue(self.phase.thermal_conductivity > 0.0)

    def test_mixtureAveraged(self):
        self.assertEqual(self.phase.transport_model, 'Mix')
        Dkm1 = self.phase.mix_diff_coeffs
        Dkm1b = self.phase.mix_diff_coeffs_mole
        Dkm1c = self.phase.mix_diff_coeffs_mass
        Dbin1 = self.phase.binary_diff_coeffs

        self.phase.transport_model = 'Multi'
        Dkm2 = self.phase.mix_diff_coeffs
        Dkm2b = self.phase.mix_diff_coeffs_mole
        Dkm2c = self.phase.mix_diff_coeffs_mass
        Dbin2 = self.phase.binary_diff_coeffs
        self.assertArrayNear(Dkm1, Dkm2)
        self.assertArrayNear(Dkm1b, Dkm2b)
        self.assertArrayNear(Dkm1c, Dkm2c)
        self.assertArrayNear(Dbin1, Dbin2)
        self.assertArrayNear(Dbin1, Dbin1.T)

    def test_multiComponent(self):
        self.assertRaises(Exception,
                          lambda: self.phase.Multi_diff_coeffs)

        self.assertArrayNear(self.phase.thermal_diff_coeffs,
                             np.zeros(self.phase.n_species))

        self.phase.transport_model = 'Multi'
        self.assertTrue(all(self.phase.multi_diff_coeffs.flat >= 0.0))
        self.assertTrue(all(self.phase.thermal_diff_coeffs.flat != 0.0))


class TestDustyGas(utilities.CanteraTest):
    def setUp(self):
        self.phase = ct.DustyGas('h2o2.xml')
        self.phase.TPX = 500.0, ct.one_atm, "O2:2.0, H2:1.0, H2O:1.0"
        self.phase.porosity = 0.2
        self.phase.tortuosity = 0.3
        self.phase.mean_pore_radius = 1e-4
        self.phase.mean_particle_diameter = 5e-4
        self.Dref = self.phase.multi_diff_coeffs

    def test_porosity(self):
        self.phase.porosity = 0.4
        D = self.phase.multi_diff_coeffs
        self.assertArrayNear(self.Dref * 2, D)

    def test_tortuosity(self):
        self.phase.tortuosity = 0.6
        D = self.phase.multi_diff_coeffs
        self.assertArrayNear(self.Dref * 0.5, D)

    # The other parameters don't have such simple relationships to the diffusion
    # coefficients, so we can't test them as easily

    def test_molar_fluxes(self):
        T1, rho1, Y1 = self.phase.TDY
        self.phase.TPX = 500.0, ct.one_atm, "O2:2.0, H2:1.001, H2O:0.999"

        T2, rho2, Y2 = self.phase.TDY

        fluxes0 = self.phase.molar_fluxes(T1, T1, rho1, rho1, Y1, Y1, 1e-4)
        self.assertArrayNear(fluxes0, np.zeros(self.phase.n_species))

        fluxes1 = self.phase.molar_fluxes(T1, T2, rho1, rho2, Y1, Y2, 1e-4)
        kH2 = self.phase.species_index('H2')
        kH2O = self.phase.species_index('H2O')
        self.assertTrue(fluxes1[kH2] < 0)
        self.assertTrue(fluxes1[kH2O] > 0)

        # Not sure why the following condition is not satisfied:
        # self.assertNear(sum(fluxes1) / sum(abs(fluxes1)), 0.0)
