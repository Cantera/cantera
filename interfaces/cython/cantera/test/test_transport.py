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

    def test_unityLewis(self):
        self.phase.transport_model = 'UnityLewis'
        alpha = self.phase.thermal_conductivity/(self.phase.density*self.phase.cp)
        Dkm_prime = self.phase.mix_diff_coeffs

        Dkm = self.phase.mix_diff_coeffs_mass

        eps = np.spacing(1) # Machine precision
        self.assertTrue(all(np.diff(Dkm) < 2*eps))
        self.assertNear(Dkm[0], alpha)
        self.assertTrue(all(np.diff(Dkm_prime) < 2*eps))
        self.assertNear(Dkm_prime[0], alpha)

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

    def test_mixDiffCoeffsChange(self):
        # This test is mainly to make code coverage in GasTransport.cpp
        # consistent by always covering the path where the binary diffusion
        # coefficients need to be updated
        Dkm1 = self.phase.mix_diff_coeffs_mole
        self.phase.TP = self.phase.T + 1, None
        Dkm2 = self.phase.mix_diff_coeffs_mole
        self.assertTrue(all(Dkm2 > Dkm1))

        Dkm1 = self.phase.mix_diff_coeffs_mass
        self.phase.TP = self.phase.T + 1, None
        Dkm2 = self.phase.mix_diff_coeffs_mass
        self.assertTrue(all(Dkm2 > Dkm1))

        Dkm1 = self.phase.mix_diff_coeffs
        self.phase.TP = self.phase.T + 1, None
        Dkm2 = self.phase.mix_diff_coeffs
        self.assertTrue(all(Dkm2 > Dkm1))

    def test_CK_mode(self):
        mu_ct = self.phase.viscosity
        self.phase.transport_model = 'CK_Mix'
        self.assertEqual(self.phase.transport_model, 'CK_Mix')
        mu_ck = self.phase.viscosity
        # values should be close, but not identical
        self.assertGreater(abs(mu_ct - mu_ck) / mu_ct, 1e-8)
        self.assertLess(abs(mu_ct - mu_ck) / mu_ct, 1e-2)

    def test_ionGas(self):
        # IonGasTransport gives the same result for a mixture
        # without ionized species
        self.phase.transport_model = 'Ion'
        Dkm1 = self.phase.mix_diff_coeffs
        Dbin1 = self.phase.binary_diff_coeffs

        self.phase.transport_model = 'Mix'
        Dkm2 = self.phase.mix_diff_coeffs
        Dbin2 = self.phase.binary_diff_coeffs
        self.assertArrayNear(Dkm1, Dkm2)
        self.assertArrayNear(Dbin1, Dbin2)

    def test_multiComponent(self):
        with self.assertRaisesRegex(ct.CanteraError, 'NotImplementedError'):
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

    def test_species_visosities(self):
        for species_name in self.phase.species_names:
            # check that species viscosity matches overall for single-species
            # state
            self.phase.X = {species_name: 1}
            self.phase.TP = 800, 2*ct.one_atm
            visc = self.phase.viscosity
            self.assertNear(self.phase[species_name].species_viscosities[0],
                            visc)
            # and ensure it doesn't change with pressure
            self.phase.TP = 800, 5*ct.one_atm
            self.assertNear(self.phase[species_name].species_viscosities[0],
                            visc)


class TestIonTransport(utilities.CanteraTest):
    def setUp(self):
        self.p = ct.one_atm
        self.T = 2237
        self.gas = ct.Solution('ch4_ion.cti')
        self.gas.X = 'O2:0.7010, H2O:0.1885, CO2:9.558e-2'
        self.gas.TP = self.T, self.p
        self.kN2 = self.gas.species_index("N2")
        self.kH3Op = self.gas.species_index("H3O+")

    def test_binary_diffusion(self):
        bdiff = self.gas.binary_diff_coeffs[self.kN2][self.kH3Op]
        self.assertNear(bdiff, 4.258e-4, 1e-4)  # Regression test

    def test_mixture_diffusion(self):
        mdiff = self.gas.mix_diff_coeffs[self.kH3Op]
        self.assertNear(mdiff, 5.057e-4, 1e-4)  # Regression test

    def test_O2_anion_mixture_diffusion(self):
        mdiff = self.gas['O2-'].mix_diff_coeffs[0]
        self.assertNear(mdiff, 2.784e-4, 1e-3)  # Regression test

    def test_mobility(self):
        mobi = self.gas.mobilities[self.kH3Op]
        self.assertNear(mobi, 2.623e-3, 1e-4)  # Regression test

    def test_update_temperature(self):
        bdiff = self.gas.binary_diff_coeffs[self.kN2][self.kH3Op]
        mdiff = self.gas.mix_diff_coeffs[self.kH3Op]
        mobi = self.gas.mobilities[self.kH3Op]
        self.gas.TP = 0.9 * self.T, self.p
        self.assertTrue(bdiff != self.gas.binary_diff_coeffs[self.kN2][self.kH3Op])
        self.assertTrue(mdiff != self.gas.mix_diff_coeffs[self.kH3Op])
        self.assertTrue(mobi != self.gas.mobilities[self.kH3Op])


class TestIonTransportYAML(TestIonTransport):
    def setUp(self):
        self.p = ct.one_atm
        self.T = 2237
        self.gas = ct.Solution('ch4_ion.yaml')
        self.gas.X = 'O2:0.7010, H2O:0.1885, CO2:9.558e-2'
        self.gas.TP = self.T, self.p
        self.kN2 = self.gas.species_index("N2")
        self.kH3Op = self.gas.species_index("H3O+")


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
            with self.assertRaisesRegex(ct.CanteraError, 'invalid geometry'):
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


class TestWaterTransport(utilities.CanteraTest):
    """
    Comparison values are taken from the NIST Chemistry WebBook. Agreement is
    limited by the use of a different equation of state here (Reynolds) than
    in the Webbook (IAPWS95), as well as a different set of coefficients for
    the transport property model. Differences are largest in the region near
    the critical point.
    """
    @classmethod
    def setUpClass(cls):
        cls.water = ct.Water()

    def check_viscosity(self, T, P, mu, rtol):
        self.water.TP = T, P
        self.assertNear(self.water.viscosity, mu, rtol)

    def check_thermal_conductivity(self, T, P, k, rtol):
        self.water.TP = T, P
        self.assertNear(self.water.thermal_conductivity, k, rtol)

    def test_viscosity_liquid(self):
        self.check_viscosity(400, 1e6, 2.1880e-4, 1e-3)
        self.check_viscosity(400, 8e6, 2.2061e-4, 1e-3)
        self.check_viscosity(620, 1.6e7, 6.7489e-5, 2e-3)
        self.check_viscosity(620, 2.8e7, 7.5684e-5, 2e-3)

    def test_thermal_conductivity_liquid(self):
        self.check_thermal_conductivity(400, 1e6, 0.68410, 1e-3)
        self.check_thermal_conductivity(400, 8e6, 0.68836, 1e-3)
        self.check_thermal_conductivity(620, 1.6e7, 0.45458, 2e-3)
        self.check_thermal_conductivity(620, 2.8e7, 0.49705, 2e-3)

    def test_viscosity_vapor(self):
        self.check_viscosity(600, 1e6, 2.1329e-5, 1e-3)
        self.check_viscosity(620, 5e6, 2.1983e-5, 1e-3)
        self.check_viscosity(620, 1.5e7, 2.2858e-5, 2e-3)

    def test_thermal_conductivity_vapor(self):
        self.check_thermal_conductivity(600, 1e6, 0.047636, 1e-3)
        self.check_thermal_conductivity(620, 5e6, 0.055781, 1e-3)
        self.check_thermal_conductivity(620, 1.5e7, 0.10524, 2e-3)

    def test_viscosity_supercritical(self):
        self.check_viscosity(660, 2.2e7, 2.7129e-5, 2e-3)
        self.check_viscosity(660, 2.54e7, 3.8212e-5, 1e-2)
        self.check_viscosity(660, 2.8e7, 5.3159e-5, 1e-2)

    def test_thermal_conductivity_supercritical(self):
        self.check_thermal_conductivity(660, 2.2e7, 0.14872, 1e-2)
        self.check_thermal_conductivity(660, 2.54e7, 0.35484, 2e-2)
        self.check_thermal_conductivity(660, 2.8e7, 0.38479, 1e-2)

class TestIAPWS95WaterTransport(utilities.CanteraTest):
    """
    Water transport properties test using the IAPWS95 equation of state. This
    results in better comparisons with data from the NIST Webbook.
    """
    @classmethod
    def setUpClass(cls):
        cls.water = ct.Solution('thermo-models.yaml', 'liquid-water')

    def check_viscosity(self, T, P, mu, rtol):
        self.water.TP = T, P
        self.assertNear(self.water.viscosity, mu, rtol)

    def check_thermal_conductivity(self, T, P, k, rtol):
        self.water.TP = T, P
        self.assertNear(self.water.thermal_conductivity, k, rtol)

    def test_viscosity_liquid(self):
        self.check_viscosity(400, 1e6, 2.1880e-4, 2e-4)
        self.check_viscosity(400, 8e6, 2.2061e-4, 2e-4)
        self.check_viscosity(620, 1.6e7, 6.7489e-5, 1e-4)
        self.check_viscosity(620, 2.8e7, 7.5684e-5, 1e-4)

    def test_thermal_conductivity_liquid(self):
        self.check_thermal_conductivity(400, 1e6, 0.68410, 1e-4)
        self.check_thermal_conductivity(400, 8e6, 0.68836, 1e-4)
        self.check_thermal_conductivity(620, 1.6e7, 0.45458, 1e-4)
        self.check_thermal_conductivity(620, 2.8e7, 0.49705, 1e-4)

    def test_viscosity_supercritical(self):
        self.check_viscosity(660, 2.2e7, 2.7129e-5, 1e-4)
        self.check_viscosity(660, 2.54e7, 3.8212e-5, 1e-4)
        self.check_viscosity(660, 2.8e7, 5.3159e-5, 1e-4)

    def test_thermal_conductivity_supercritical(self):
        self.check_thermal_conductivity(660, 2.2e7, 0.14872, 1e-4)
        self.check_thermal_conductivity(660, 2.54e7, 0.35484, 1e-4)
        self.check_thermal_conductivity(660, 2.8e7, 0.38479, 1e-4)


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


class TestIonGasTransportData(utilities.CanteraTest):
    def setUp(self):
        self.gas = ct.Solution('ch4_ion.cti')

    def test_read_ion(self):
        tr = self.gas.species('N2').transport
        self.assertNear(tr.dispersion_coefficient, 2.995e-50)
        self.assertNear(tr.quadrupole_polarizability, 3.602e-50)

    def test_set_customary_units(self):
        tr1 = ct.GasTransportData()
        tr1.set_customary_units('linear', 3.62, 97.53, 1.76,
                                dispersion_coefficient = 2.995,
                                quadrupole_polarizability = 3.602)
        tr2 = self.gas.species('N2').transport
        self.assertNear(tr1.dispersion_coefficient, tr2.dispersion_coefficient)
        self.assertNear(tr1.quadrupole_polarizability, tr2.quadrupole_polarizability)
