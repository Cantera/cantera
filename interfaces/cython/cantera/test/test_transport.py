import numpy as np

import cantera as ct
from . import utilities
import copy


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
        with self.assertRaises(ct.CanteraError):
            self.phase.multi_diff_coeffs

        self.assertArrayNear(self.phase.thermal_diff_coeffs,
                             np.zeros(self.phase.n_species))

        self.phase.transport_model = 'Multi'
        self.assertTrue(all(self.phase.multi_diff_coeffs.flat >= 0.0))
        self.assertTrue(all(self.phase.thermal_diff_coeffs.flat != 0.0))

    def test_add_species_mix(self):
        S = {s.name: s for s in ct.Species.listFromFile('gri30.xml')}

        base = ['H', 'H2', 'OH', 'O2', 'AR']
        extra = ['H2O', 'CH4']

        state = 500, 2e5, 'H2:0.4, O2:0.29, CH4:0.01, H2O:0.3'

        gas1 = ct.Solution(thermo='IdealGas', species=[S[s] for s in base+extra])
        gas1.transport_model = 'Mix'
        gas1.TPX = state

        gas2 = ct.Solution(thermo='IdealGas', species=[S[s] for s in base])
        gas2.transport_model = 'Mix'
        for s in extra:
            gas2.add_species(S[s])
        gas2.TPX = state

        self.assertNear(gas1.viscosity, gas2.viscosity)
        self.assertNear(gas1.thermal_conductivity, gas2.thermal_conductivity)
        self.assertArrayNear(gas1.binary_diff_coeffs, gas2.binary_diff_coeffs)
        self.assertArrayNear(gas1.mix_diff_coeffs, gas2.mix_diff_coeffs)

    def test_add_species_multi(self):
        S = {s.name: s for s in ct.Species.listFromFile('gri30.xml')}

        base = ['H', 'H2', 'OH', 'O2', 'AR', 'N2']
        extra = ['H2O', 'CH4']

        state = 500, 2e5, 'H2:0.3, O2:0.28, CH4:0.02, H2O:0.3, N2:0.1'

        gas1 = ct.Solution(thermo='IdealGas', species=[S[s] for s in base+extra])
        gas1.transport_model = 'Multi'
        gas1.TPX = state

        gas2 = ct.Solution(thermo='IdealGas', species=[S[s] for s in base])
        gas2.transport_model = 'Multi'
        for s in extra:
            gas2.add_species(S[s])
        gas2.TPX = state

        self.assertNear(gas1.thermal_conductivity, gas2.thermal_conductivity)
        self.assertArrayNear(gas1.multi_diff_coeffs, gas2.multi_diff_coeffs)


class TestTransportGeometryFlags(utilities.CanteraTest):
    phase_data = """
units(length="cm", time="s", quantity="mol", act_energy="cal/mol")

ideal_gas(name="test",
    elements="O  H  E",
    species="H2  H  H2O  OHp  E",
    initial_state=state(temperature=300.0, pressure=OneAtm),
    transport='Mix'
)

species(name="H2",
    atoms=" H:2 ",
    thermo=const_cp(t0=1000, h0=51.7, s0=19.5, cp0=8.41),
    transport=gas_transport(
                geom="{H2}",
                diam=2.92, well_depth=38.00, polar=0.79, rot_relax=280.00)
)

species(name="H",
    atoms=" H:1 ",
    thermo=const_cp(t0=1000, h0=51.7, s0=19.5, cp0=8.41),
    transport=gas_transport(
                geom="{H}",
                diam=2.05, well_depth=145.00)
)
species(name="H2O",
    atoms=" H:2  O:1 ",
    thermo=const_cp(t0=1000, h0=51.7, s0=19.5, cp0=8.41),
    transport=gas_transport(
                geom="{H2O}",
                diam=2.60, well_depth=572.40, dipole=1.84, rot_relax=4.00)
)

species(name="OHp",
    atoms=" H:1  O:1  E:-1",
    thermo=const_cp(t0=1000, h0=51.7, s0=19.5, cp0=8.41),
    transport=gas_transport(
                geom="{OHp}",
                diam=2.60, well_depth=572.40, dipole=1.84, rot_relax=4.00)
)

species(name="E",
    atoms=" E:1 ",
    thermo=const_cp(t0=1000, h0=0, s0=0, cp0=0),
    transport=gas_transport(
                geom="{E}", diam=0.01, well_depth=1.00)
)
"""
    def test_bad_geometry(self):
        good = {'H':'atom', 'H2':'linear', 'H2O':'nonlinear', 'OHp':'linear',
                'E':'atom'}
        ct.Solution(source=self.phase_data.format(**good))

        bad = [{'H':'linear'}, {'H':'nonlinear'}, {'H2':'atom'},
               {'H2':'nonlinear'}, {'H2O':'atom'}, {'OHp':'atom'},
               {'OHp':'nonlinear'}, {'E':'linear'}]
        for geoms in bad:
            test = copy.copy(good)
            test.update(geoms)
            with self.assertRaises(ct.CanteraError):
                ct.Solution(source=self.phase_data.format(**test))


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


class TestTransportData(utilities.CanteraTest):
    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        cls.gas = ct.Solution('h2o2.xml')
        cls.gas.X = 'H2O:0.6, H2:0.4'

    def test_read(self):
        tr = self.gas.species('H2').transport
        self.assertEqual(tr.geometry, 'linear')
        self.assertNear(tr.diameter, 2.92e-10)
        self.assertNear(tr.well_depth, 38.0 * ct.boltzmann)
        self.assertNear(tr.polarizability, 0.79e-30)
        self.assertNear(tr.rotational_relaxation, 280)

    def test_set_customary_units(self):
        tr1 = ct.GasTransportData()
        tr1.set_customary_units('nonlinear', 2.60, 572.40, 1.84, 0.0, 4.00)
        tr2 = self.gas.species('H2O').transport
        self.assertEqual(tr1.geometry, tr2.geometry)
        self.assertNear(tr1.diameter, tr2.diameter)
        self.assertNear(tr1.well_depth, tr2.well_depth)
        self.assertNear(tr1.dipole, tr2.dipole)
        self.assertNear(tr1.polarizability, tr2.polarizability)
        self.assertNear(tr1.rotational_relaxation, tr2.rotational_relaxation)
