import copy
import numpy as np
import pytest
from pytest import approx

import cantera as ct


@pytest.fixture(scope='function')
def setup_transport(request):
    # Create instance-level data for each test
    request.cls.phase = ct.Solution('h2o2.yaml')
    request.cls.phase.X = [0.1, 1e-4, 1e-5, 0.2, 2e-4, 0.3, 1e-6, 5e-5, 1e-6, 0.4]
    request.cls.phase.TP = 800, 2*ct.one_atm

@pytest.mark.usefixtures("setup_transport")
class TestTransport:

    def test_scalar_properties(self):
        assert self.phase.viscosity > 0.0
        assert self.phase.thermal_conductivity > 0.0

    def test_unityLewis(self):
        self.phase.transport_model = 'unity-Lewis-number'
        alpha = self.phase.thermal_conductivity/(self.phase.density*self.phase.cp)
        Dkm_prime = self.phase.mix_diff_coeffs

        Dkm = self.phase.mix_diff_coeffs_mass

        eps = np.spacing(1) # Machine precision
        assert all(np.diff(Dkm) < 2*eps)
        assert Dkm[0] == approx(alpha)
        assert all(np.diff(Dkm_prime) < 2*eps)
        assert Dkm_prime[0] == approx(alpha)

    def test_mixtureAveraged(self):
        assert self.phase.transport_model == 'mixture-averaged'
        Dkm1 = self.phase.mix_diff_coeffs
        Dkm1b = self.phase.mix_diff_coeffs_mole
        Dkm1c = self.phase.mix_diff_coeffs_mass
        Dbin1 = self.phase.binary_diff_coeffs

        self.phase.transport_model = 'multicomponent'
        Dkm2 = self.phase.mix_diff_coeffs
        Dkm2b = self.phase.mix_diff_coeffs_mole
        Dkm2c = self.phase.mix_diff_coeffs_mass
        Dbin2 = self.phase.binary_diff_coeffs
        assert Dkm1 == approx(Dkm2)
        assert Dkm1b == approx(Dkm2b)
        assert Dkm1c == approx(Dkm2c)
        assert Dbin1 == approx(Dbin2)
        assert Dbin1 == approx(Dbin1.T)

    def test_mixDiffCoeffsChange(self):
        # This test is mainly to make code coverage in GasTransport.cpp
        # consistent by always covering the path where the binary diffusion
        # coefficients need to be updated
        Dkm1 = self.phase.mix_diff_coeffs_mole
        self.phase.TP = self.phase.T + 1, None
        Dkm2 = self.phase.mix_diff_coeffs_mole
        assert all(Dkm2 > Dkm1)

        Dkm1 = self.phase.mix_diff_coeffs_mass
        self.phase.TP = self.phase.T + 1, None
        Dkm2 = self.phase.mix_diff_coeffs_mass
        assert all(Dkm2 > Dkm1)

        Dkm1 = self.phase.mix_diff_coeffs
        self.phase.TP = self.phase.T + 1, None
        Dkm2 = self.phase.mix_diff_coeffs
        assert all(Dkm2 > Dkm1)

    def test_CK_mode(self):
        mu_ct = self.phase.viscosity
        self.phase.transport_model = 'mixture-averaged-CK'
        assert self.phase.transport_model == 'mixture-averaged-CK'
        mu_ck = self.phase.viscosity
        # values should be close, but not identical
        assert abs(mu_ct - mu_ck) / mu_ct > 1e-8
        assert abs(mu_ct - mu_ck) / mu_ct < 1e-2

    def test_ionGas(self):
        # IonGasTransport gives the same result for a mixture
        # without ionized species
        self.phase.transport_model = 'ionized-gas'
        Dkm1 = self.phase.mix_diff_coeffs
        Dbin1 = self.phase.binary_diff_coeffs

        self.phase.transport_model = 'mixture-averaged'
        Dkm2 = self.phase.mix_diff_coeffs
        Dbin2 = self.phase.binary_diff_coeffs
        assert Dkm1 == approx(Dkm2)
        assert Dbin1 == approx(Dbin2)

    def test_multiComponent(self):
        with pytest.raises(NotImplementedError):
            self.phase.multi_diff_coeffs

        assert self.phase.thermal_diff_coeffs == approx(np.zeros(self.phase.n_species))

        self.phase.transport_model = 'multicomponent'
        assert all(self.phase.multi_diff_coeffs.flat >= 0.0)
        assert all(self.phase.thermal_diff_coeffs.flat != 0.0)

    def test_add_species_mix(self, cantera_data_path):
        yaml = (cantera_data_path / "gri30.yaml").read_text()
        S = {s.name: s for s in ct.Species.list_from_yaml(yaml, "species")}

        base = ['H', 'H2', 'OH', 'O2', 'AR']
        extra = ['H2O', 'CH4']

        state = 500, 2e5, 'H2:0.4, O2:0.29, CH4:0.01, H2O:0.3'

        gas1 = ct.Solution(thermo='ideal-gas', species=[S[s] for s in base+extra])
        gas1.transport_model = 'mixture-averaged'
        gas1.TPX = state

        gas2 = ct.Solution(thermo='ideal-gas', species=[S[s] for s in base])
        gas2.transport_model = 'mixture-averaged'
        for s in extra:
            gas2.add_species(S[s])
        gas2.TPX = state

        assert gas1.viscosity == approx(gas2.viscosity)
        assert gas1.thermal_conductivity == approx(gas2.thermal_conductivity)
        assert gas1.binary_diff_coeffs == approx(gas2.binary_diff_coeffs)
        assert gas1.mix_diff_coeffs == approx(gas2.mix_diff_coeffs)

    def test_add_species_multi(self, cantera_data_path):
        yaml = (cantera_data_path / "gri30.yaml").read_text()
        S = {s.name: s for s in ct.Species.list_from_yaml(yaml, "species")}

        base = ['H', 'H2', 'OH', 'O2', 'AR', 'N2']
        extra = ['H2O', 'CH4']

        state = 500, 2e5, 'H2:0.3, O2:0.28, CH4:0.02, H2O:0.3, N2:0.1'

        gas1 = ct.Solution(thermo='ideal-gas', species=[S[s] for s in base+extra])
        gas1.transport_model = 'multicomponent'
        gas1.TPX = state

        gas2 = ct.Solution(thermo='ideal-gas', species=[S[s] for s in base])
        gas2.transport_model = 'multicomponent'
        for s in extra:
            gas2.add_species(S[s])
        gas2.TPX = state

        assert gas1.thermal_conductivity == approx(gas2.thermal_conductivity)
        assert gas1.multi_diff_coeffs == approx(gas2.multi_diff_coeffs)

    def test_species_visosities(self):
        for species_name in self.phase.species_names:
            # check that species viscosity matches overall for single-species
            # state
            self.phase.X = {species_name: 1}
            self.phase.TP = 800, 2*ct.one_atm
            visc = self.phase.viscosity
            assert self.phase[species_name].species_viscosities[0] == approx(visc)
            # and ensure it doesn't change with pressure
            self.phase.TP = 800, 5*ct.one_atm
            assert self.phase[species_name].species_viscosities[0] == approx(visc)

    def test_transport_polynomial_fits_viscosity(self):
        visc1_h2o = self.phase['H2O'].species_viscosities[0]
        mu_poly_h2o = self.phase.get_viscosity_polynomial(self.phase.species_index("H2O"))
        visc1_h2 = self.phase['H2'].species_viscosities[0]
        mu_poly_h2 = self.phase.get_viscosity_polynomial(self.phase.species_index('H2'))
        self.phase.set_viscosity_polynomial(self.phase.species_index('H2'), mu_poly_h2o)
        visc2_h2 = self.phase['H2'].species_viscosities[0]
        self.phase.set_viscosity_polynomial(self.phase.species_index('H2'), mu_poly_h2)
        visc3_h2 = self.phase['H2'].species_viscosities[0]
        assert visc1_h2o != visc1_h2
        assert visc1_h2o == visc2_h2
        assert visc1_h2 == visc3_h2

    def test_transport_polynomial_fits_conductivity(self):
        self.phase.X = {'O2': 1}
        cond1_o2 = self.phase.thermal_conductivity
        lambda_poly_o2 = self.phase.get_thermal_conductivity_polynomial(self.phase.species_index("O2"))
        self.phase.X = {"H2": 1}
        cond1_h2 = self.phase.thermal_conductivity
        lambda_poly_h2 = self.phase.get_thermal_conductivity_polynomial(self.phase.species_index('H2'))
        self.phase.set_thermal_conductivity_polynomial(self.phase.species_index('H2'), lambda_poly_o2)
        cond2_h2 = self.phase.thermal_conductivity
        self.phase.set_thermal_conductivity_polynomial(self.phase.species_index('H2'), lambda_poly_h2)
        cond3_h2 = self.phase.thermal_conductivity
        assert cond1_o2 != cond1_h2
        assert cond1_o2 == cond2_h2
        assert cond1_h2 == cond3_h2

    def test_transport_polynomial_fits_diffusion(self):
        D12 = self.phase.binary_diff_coeffs[1, 2]
        D23 = self.phase.binary_diff_coeffs[2, 3]
        bd_poly_12 = self.phase.get_binary_diff_coeffs_polynomial(1, 2)
        bd_poly_23 = self.phase.get_binary_diff_coeffs_polynomial(2, 3)
        self.phase.set_binary_diff_coeffs_polynomial(1, 2, bd_poly_23)
        self.phase.set_binary_diff_coeffs_polynomial(2, 3, bd_poly_12)
        D12mod = self.phase.binary_diff_coeffs[1, 2]
        D23mod = self.phase.binary_diff_coeffs[2, 3]
        self.phase.set_binary_diff_coeffs_polynomial(1, 2, bd_poly_12)
        self.phase.set_binary_diff_coeffs_polynomial(2, 3, bd_poly_23)
        D12new = self.phase.binary_diff_coeffs[1, 2]
        D23new = self.phase.binary_diff_coeffs[2, 3]
        assert D12 != D23
        assert D12 == D23mod
        assert D23 == D12mod
        assert D12 == D12new
        assert D23 == D23new


@pytest.fixture(scope='function')
def setup_ion_transport(request):
    # Create instance-level data for each test
    request.cls.p = ct.one_atm
    request.cls.T = 2237
    request.cls.gas = ct.Solution('ch4_ion.yaml')
    request.cls.gas.X = 'O2:0.7010, H2O:0.1885, CO2:9.558e-2'
    request.cls.gas.TP = request.cls.T, request.cls.p
    request.cls.kN2 = request.cls.gas.species_index("N2")
    request.cls.kH3Op = request.cls.gas.species_index("H3O+")

@pytest.mark.usefixtures("setup_ion_transport")
class TestIonTransport:

    def test_binary_diffusion(self):
        bdiff = self.gas.binary_diff_coeffs[self.kN2][self.kH3Op]
        assert bdiff == approx(4.258e-4, rel=1e-4)  # Regression test

    def test_mixture_diffusion(self):
        mdiff = self.gas.mix_diff_coeffs[self.kH3Op]
        assert mdiff == approx(5.057e-4, rel=1e-4)  # Regression test

    def test_O2_anion_mixture_diffusion(self):
        mdiff = self.gas['O2-'].mix_diff_coeffs[0]
        assert mdiff == approx(2.784e-4, rel=1e-3)  # Regression test

    def test_mobility(self):
        mobi = self.gas.mobilities[self.kH3Op]
        assert mobi == approx(2.623e-3, rel=1e-4)  # Regression test

    def test_update_temperature(self):
        bdiff = self.gas.binary_diff_coeffs[self.kN2][self.kH3Op]
        mdiff = self.gas.mix_diff_coeffs[self.kH3Op]
        mobi = self.gas.mobilities[self.kH3Op]
        self.gas.TP = 0.9 * self.T, self.p
        assert bdiff != self.gas.binary_diff_coeffs[self.kN2][self.kH3Op]
        assert mdiff != self.gas.mix_diff_coeffs[self.kH3Op]
        assert mobi != self.gas.mobilities[self.kH3Op]


class TestTransportGeometryFlags:
    species_data = """
    - name: H2
      composition: {{H: 2}}
      thermo: &dummy-thermo
        {{model: constant-cp, T0: 1000, h0: 51.7, s0: 19.5, cp0: 8.41}}
      transport:
        model: gas
        geometry: {H2}
        diameter: 2.92
        well-depth: 38.00
        polarizability: 0.79
        rotational-relaxation: 280.0
    - name: H
      composition: {{H: 1}}
      thermo: *dummy-thermo
      transport:
        model: gas
        geometry: {H}
        diameter: 2.05
        well-depth: 145.00
    - name: H2O
      composition: {{H: 2, O: 1}}
      thermo: *dummy-thermo
      transport:
        model: gas
        geometry: {H2O}
        diameter: 2.60
        well-depth: 572.40
        dipole: 1.84
        rotational-relaxation: 4.0
    - name: OHp
      composition: {{H: 1, O: 1, E: -1}}
      thermo: *dummy-thermo
      transport:
        model: gas
        geometry: {OHp}
        diameter: 2.60
        well-depth: 572.40
        dipole: 1.84
        rotational-relaxation: 4.0
    - name: E
      composition: {{E: 1}}
      thermo: *dummy-thermo
      transport:
        model: gas
        geometry: {E}
        diameter: 0.01
        well-depth: 1.0
    """

    def test_bad_geometry(self):
        good = {'H':'atom', 'H2':'linear', 'H2O':'nonlinear', 'OHp':'linear',
                'E':'atom'}
        ct.Species.list_from_yaml(self.species_data.format(**good))

        bad = [{'H':'linear'}, {'H':'nonlinear'}, {'H2':'atom'},
               {'H2':'nonlinear'}, {'H2O':'atom'}, {'OHp':'atom'},
               {'OHp':'nonlinear'}, {'E':'linear'}]
        for geoms in bad:
            test = copy.copy(good)
            test.update(geoms)
            with pytest.raises(ct.CanteraError, match='invalid geometry'):
                ct.Species.list_from_yaml(self.species_data.format(**test))


@pytest.fixture(scope='function')
def setup_dusty_gas(request):
    # Create instance-level data for each test
    request.cls.phase = ct.DustyGas("h2o2.yaml")
    request.cls.phase.TPX = 500.0, ct.one_atm, "O2:2.0, H2:1.0, H2O:1.0"
    request.cls.phase.porosity = 0.2
    request.cls.phase.tortuosity = 0.3
    request.cls.phase.mean_pore_radius = 1e-4
    request.cls.phase.mean_particle_diameter = 5e-4
    request.cls.Dref = request.cls.phase.multi_diff_coeffs

@pytest.mark.usefixtures("setup_dusty_gas")
class TestDustyGas:

    def test_porosity(self):
        self.phase.porosity = 0.4
        D = self.phase.multi_diff_coeffs
        assert self.Dref * 2 == approx(D)

    def test_tortuosity(self):
        self.phase.tortuosity = 0.6
        D = self.phase.multi_diff_coeffs
        assert self.Dref * 0.5 == approx(D)

    # The other parameters don't have such simple relationships to the diffusion
    # coefficients, so we can't test them as easily

    def test_molar_fluxes(self):
        T1, rho1, Y1 = self.phase.TDY
        self.phase.TPX = 500.0, ct.one_atm, "O2:2.0, H2:1.001, H2O:0.999"

        T2, rho2, Y2 = self.phase.TDY

        fluxes0 = self.phase.molar_fluxes(T1, T1, rho1, rho1, Y1, Y1, 1e-4)
        assert fluxes0 == approx(np.zeros(self.phase.n_species))

        fluxes1 = self.phase.molar_fluxes(T1, T2, rho1, rho2, Y1, Y2, 1e-4)
        kH2 = self.phase.species_index('H2')
        kH2O = self.phase.species_index('H2O')
        assert fluxes1[kH2] < 0
        assert fluxes1[kH2O] > 0

        # Not sure why the following condition is not satisfied:
        #assert sum(fluxes1) == approx(0.0)
        #assert sum(fluxes1) / sum(abs(fluxes1)) == approx(0.0)

    def test_thermal_conductivity(self):
        gas1 = ct.Solution("h2o2.yaml", transport_model="multicomponent")
        gas1.TPX = self.phase.TPX

        assert self.phase.thermal_conductivity == gas1.thermal_conductivity


@pytest.fixture(scope='function')
def water(request):
    request.cls.water = ct.Water()

@pytest.mark.usefixtures("water")
class TestWaterTransport:
    """
    Comparison values are taken from the NIST Chemistry WebBook. Agreement is
    limited by the use of a different equation of state here (Reynolds) than
    in the Webbook (IAPWS95), as well as a different set of coefficients for
    the transport property model. Differences are largest in the region near
    the critical point.
    """

    def check_viscosity(self, T, P, mu, rtol):
        self.water.TP = T, P
        assert self.water.viscosity == approx(mu, rel=rtol)

    def check_thermal_conductivity(self, T, P, k, rtol):
        self.water.TP = T, P
        assert self.water.thermal_conductivity == approx(k, rel=rtol)

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


@pytest.fixture(scope='function')
def setup_IAPWS95Water(request):
    # Create instance-level data for each test
    request.cls.water = ct.Solution('thermo-models.yaml', 'liquid-water')

@pytest.mark.usefixtures("setup_IAPWS95Water")
class TestIAPWS95WaterTransport:
    """
    Water transport properties test using the IAPWS95 equation of state. This
    results in better comparisons with data from the NIST Webbook.
    """

    def check_viscosity(self, T, P, mu, rtol):
        self.water.TP = T, P
        assert self.water.viscosity == approx(mu, rel=rtol)

    def check_thermal_conductivity(self, T, P, k, rtol):
        self.water.TP = T, P
        assert self.water.thermal_conductivity == approx(k, rel=rtol)

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


@pytest.fixture(scope='function')
def setup_transport_data(request):
    # Create instance-level data for each test
    request.cls.gas = ct.Solution("h2o2.yaml")
    request.cls.gas.X = 'H2O:0.6, H2:0.4'

@pytest.mark.usefixtures("setup_transport_data")
class TestTransportData:

    def test_read(self):
        tr = self.gas.species('H2').transport
        assert tr.geometry == 'linear'
        assert tr.diameter == approx(2.92e-10)
        assert tr.well_depth == approx(38.0 * ct.boltzmann)
        assert tr.polarizability == approx(0.79e-30)
        assert tr.rotational_relaxation == approx(280)

    def test_set_customary_units(self):
        tr1 = ct.GasTransportData()
        tr1.set_customary_units('nonlinear', 2.60, 572.40, 1.84, 0.0, 4.00)
        tr2 = self.gas.species('H2O').transport
        assert tr1.geometry == tr2.geometry
        assert tr1.diameter == approx(tr2.diameter)
        assert tr1.well_depth == approx(tr2.well_depth)
        assert tr1.dipole == approx(tr2.dipole)
        assert tr1.polarizability == approx(tr2.polarizability)
        assert tr1.rotational_relaxation == approx(tr2.rotational_relaxation)


@pytest.fixture(scope='function')
def setup_ion_transport_data(request):
    # Create instance-level data for each test
    request.cls.gas = ct.Solution("ch4_ion.yaml")

@pytest.mark.usefixtures("setup_ion_transport_data")
class TestIonGasTransportData:

    def test_read_ion(self):
        tr = self.gas.species('N2').transport
        assert tr.dispersion_coefficient == approx(2.995e-50)
        assert tr.quadrupole_polarizability == approx(3.602e-50)

    def test_set_customary_units(self):
        tr1 = ct.GasTransportData()
        tr1.set_customary_units('linear', 3.62, 97.53, 1.76,
                                dispersion_coefficient = 2.995,
                                quadrupole_polarizability = 3.602)
        tr2 = self.gas.species('N2').transport
        assert tr1.dispersion_coefficient == approx(tr2.dispersion_coefficient)
        assert tr1.quadrupole_polarizability == approx(tr2.quadrupole_polarizability)
