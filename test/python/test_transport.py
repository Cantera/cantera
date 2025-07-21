import copy
import numpy as np
import pytest
from pytest import approx

import cantera as ct
from .utilities import compareProfiles

class TestTransport:

    @pytest.fixture(scope='function')
    def phase(self):
        phase = ct.Solution('h2o2.yaml')
        phase.X = [0.1, 1e-4, 1e-5, 0.2, 2e-4, 0.3, 1e-6, 5e-5, 1e-6, 0.4]
        phase.TP = 800, 2*ct.one_atm
        return phase

    def test_scalar_properties(self, phase):
        assert phase.viscosity > 0.0
        assert phase.thermal_conductivity > 0.0

    def test_unityLewis(self, phase):
        phase.transport_model = 'unity-Lewis-number'
        alpha = phase.thermal_conductivity/(phase.density*phase.cp)
        Dkm_prime = phase.mix_diff_coeffs

        Dkm = phase.mix_diff_coeffs_mass

        eps = np.spacing(1) # Machine precision
        assert all(np.diff(Dkm) < 2*eps)
        assert Dkm[0] == approx(alpha)
        assert all(np.diff(Dkm_prime) < 2*eps)
        assert Dkm_prime[0] == approx(alpha)

    def test_mixtureAveraged(self, phase):
        assert phase.transport_model == 'mixture-averaged'
        Dkm1 = phase.mix_diff_coeffs
        Dkm1b = phase.mix_diff_coeffs_mole
        Dkm1c = phase.mix_diff_coeffs_mass
        Dbin1 = phase.binary_diff_coeffs
        assert all(phase.thermal_diff_coeffs.flat != 0.0)
        assert np.sum(phase.thermal_diff_coeffs) == approx(0.0)

        phase.transport_model = 'multicomponent'
        Dkm2 = phase.mix_diff_coeffs
        Dkm2b = phase.mix_diff_coeffs_mole
        Dkm2c = phase.mix_diff_coeffs_mass
        Dbin2 = phase.binary_diff_coeffs
        assert Dkm1 == approx(Dkm2)
        assert Dkm1b == approx(Dkm2b)
        assert Dkm1c == approx(Dkm2c)
        assert Dbin1 == approx(Dbin2)
        assert Dbin1 == approx(Dbin1.T)

    def test_multicomponent(self, phase):
        phase.transport_model = 'multicomponent'
        Dmulti = phase.multi_diff_coeffs

        assert Dmulti.shape == (phase.n_species, phase.n_species)
        assert all(Dmulti.diagonal() == 0)

        # The sum of mass fluxes should equal zero for any consistent vector of mole
        # fraction gradients (that is, where the sum of mole fraction gradients is zero)
        W = phase.molecular_weights
        W_mean = phase.mean_molecular_weight
        rho = phase.density

        net_flux = 0.0
        flux_scale = 0.0
        grad_X = np.array([0.1, 0.2, -0.3, 0.4, -0.6, 0.2, 0.1, -0.2, 0.5, -0.4])
        rng = np.random.default_rng(123)
        for _ in range(10):
            rng.shuffle(grad_X)
            for k in range(phase.n_species):
                flux_k = rho * W[k] / W_mean**2 * sum(W * Dmulti[k,:] * grad_X)
                net_flux += flux_k
                flux_scale += abs(flux_k)

            assert net_flux == approx(0.0, abs=1e-14 * flux_scale)

    def test_mixDiffCoeffsChange(self, phase):
        # This test is mainly to make code coverage in GasTransport.cpp
        # consistent by always covering the path where the binary diffusion
        # coefficients need to be updated
        Dkm1 = phase.mix_diff_coeffs_mole
        phase.TP = phase.T + 1, None
        Dkm2 = phase.mix_diff_coeffs_mole
        assert all(Dkm2 > Dkm1)

        Dkm1 = phase.mix_diff_coeffs_mass
        phase.TP = phase.T + 1, None
        Dkm2 = phase.mix_diff_coeffs_mass
        assert all(Dkm2 > Dkm1)

        Dkm1 = phase.mix_diff_coeffs
        phase.TP = phase.T + 1, None
        Dkm2 = phase.mix_diff_coeffs
        assert all(Dkm2 > Dkm1)

    def test_CK_mode(self, phase):
        mu_ct = phase.viscosity
        cond_ct = phase.thermal_conductivity
        diff_ct = phase.binary_diff_coeffs
        err_ct = phase.transport_fitting_errors
        phase.transport_model = 'mixture-averaged-CK'
        assert phase.transport_model == 'mixture-averaged-CK'
        mu_ck = phase.viscosity
        cond_ck = phase.thermal_conductivity
        diff_ck = phase.binary_diff_coeffs

        err_ck = phase.transport_fitting_errors
        # values should be close, but not identical
        assert mu_ck == approx(mu_ct, rel=1e-2)
        assert mu_ck != approx(mu_ct, rel=1e-8)
        assert cond_ck == approx(cond_ct, rel=1e-2)
        assert cond_ck != approx(cond_ct, rel=1e-8)
        assert diff_ck == approx(diff_ct, rel=1e-2)
        for (i, j), Dij in np.ndenumerate(diff_ck):
            assert Dij != approx(diff_ct[i,j], rel=1e-8), (i, j)

        # Cantera's fits should be an improvement in all cases
        for key in err_ct:
            assert err_ct[key] < err_ck[key]

    def test_ionized_gas_with_no_ions(self, phase):
        # IonGasTransport gives the same result for a mixture
        # without ionized species
        phase.transport_model = 'ionized-gas'
        Dkm1 = phase.mix_diff_coeffs
        Dbin1 = phase.binary_diff_coeffs

        phase.transport_model = 'mixture-averaged'
        Dkm2 = phase.mix_diff_coeffs
        Dbin2 = phase.binary_diff_coeffs
        assert Dkm1 == approx(Dkm2)
        assert Dbin1 == approx(Dbin2)

    def test_ionized_low_T(self):
        """ (C10H8, O2-) interaction exercises low T* range of Stockmayer potential """
        phase = ct.Solution('ET_test.yaml')
        kO2m = phase.species_index("O2^-")
        kNaphthalene = phase.species_index("C10H8")
        # Regression test values
        phase.TP = 300, ct.one_atm
        Dbin = phase.binary_diff_coeffs
        assert Dbin[kO2m, kNaphthalene] == approx(2.18902175e-06)
        phase.TP = 350, ct.one_atm
        Dbin = phase.binary_diff_coeffs
        assert Dbin[kO2m, kNaphthalene] == approx(2.92899733e-06)

    def test_multiComponent(self, phase):
        with pytest.raises(NotImplementedError):
            phase.multi_diff_coeffs

        phase.transport_model = 'multicomponent'
        assert all(phase.multi_diff_coeffs.flat >= 0.0)
        assert all(phase.thermal_diff_coeffs.flat != 0.0)

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

    def test_species_visosities(self, phase):
        for species_name in phase.species_names:
            # check that species viscosity matches overall for single-species
            # state
            phase.X = {species_name: 1}
            phase.TP = 800, 2*ct.one_atm
            visc = phase.viscosity
            assert phase[species_name].species_viscosities[0] == approx(visc)
            # and ensure it doesn't change with pressure
            phase.TP = 800, 5*ct.one_atm
            assert phase[species_name].species_viscosities[0] == approx(visc)

    def test_transport_polynomial_fits_viscosity(self, phase):
        with pytest.raises(ct.CanteraError, match='IndexError'):
            phase.get_viscosity_polynomial(58)
        visc1_h2o = phase['H2O'].species_viscosities[0]
        mu_poly_h2o = phase.get_viscosity_polynomial(phase.species_index("H2O"))
        visc1_h2 = phase['H2'].species_viscosities[0]
        mu_poly_h2 = phase.get_viscosity_polynomial(phase.species_index('H2'))
        phase.set_viscosity_polynomial(phase.species_index('H2'), mu_poly_h2o)
        visc2_h2 = phase['H2'].species_viscosities[0]
        phase.set_viscosity_polynomial(phase.species_index('H2'), mu_poly_h2)
        visc3_h2 = phase['H2'].species_viscosities[0]
        assert visc1_h2o != visc1_h2
        assert visc1_h2o == visc2_h2
        assert visc1_h2 == visc3_h2

    def test_transport_polynomial_fits_conductivity(self, phase):
        phase.X = {'O2': 1}
        cond1_o2 = phase.thermal_conductivity
        lambda_poly_o2 = phase.get_thermal_conductivity_polynomial(phase.species_index("O2"))
        phase.X = {"H2": 1}
        cond1_h2 = phase.thermal_conductivity
        lambda_poly_h2 = phase.get_thermal_conductivity_polynomial(phase.species_index('H2'))
        phase.set_thermal_conductivity_polynomial(phase.species_index('H2'), lambda_poly_o2)
        cond2_h2 = phase.thermal_conductivity
        phase.set_thermal_conductivity_polynomial(phase.species_index('H2'), lambda_poly_h2)
        cond3_h2 = phase.thermal_conductivity
        assert cond1_o2 != cond1_h2
        assert cond1_o2 == cond2_h2
        assert cond1_h2 == cond3_h2

    def test_transport_polynomial_fits_diffusion(self, phase):
        D12 = phase.binary_diff_coeffs[1, 2]
        D23 = phase.binary_diff_coeffs[2, 3]
        bd_poly_12 = phase.get_binary_diff_coeffs_polynomial(1, 2)
        bd_poly_23 = phase.get_binary_diff_coeffs_polynomial(2, 3)
        phase.set_binary_diff_coeffs_polynomial(1, 2, bd_poly_23)
        phase.set_binary_diff_coeffs_polynomial(2, 3, bd_poly_12)
        D12mod = phase.binary_diff_coeffs[1, 2]
        D23mod = phase.binary_diff_coeffs[2, 3]
        phase.set_binary_diff_coeffs_polynomial(1, 2, bd_poly_12)
        phase.set_binary_diff_coeffs_polynomial(2, 3, bd_poly_23)
        with pytest.raises(ct.CanteraError, match='IndexError'):
            phase.set_binary_diff_coeffs_polynomial(2, 33, bd_poly_23)
        D12new = phase.binary_diff_coeffs[1, 2]
        D23new = phase.binary_diff_coeffs[2, 3]
        assert D12 != D23
        assert D12 == D23mod
        assert D23 == D12mod
        assert D12 == D12new
        assert D23 == D23new

    def test_transport_polynomial_fits_collision_integrals(self, phase):
        kO2 = phase.species_index("O2")
        kH2O = phase.species_index("H2O")  # unique poly because of dipole moment
        phase.transport_model = 'multicomponent'
        coll_polys_H2O = phase.get_collision_integral_polynomials(kH2O, kH2O)
        coll_polys_O2 = phase.get_collision_integral_polynomials(kO2, kO2)

        def get_cond(species):
            phase.TPX = 400, 2 * ct.one_atm, {species: 1.0}
            return phase.thermal_conductivity

        cond1_O2 = get_cond("O2")
        cond1_OH = get_cond("OH")
        phase.set_collision_integral_polynomial(kO2, kO2, *coll_polys_H2O, actualT=True)
        assert get_cond("O2") != cond1_O2  # different
        assert get_cond("OH") == cond1_OH  # unchanged; normally shares poly with O2

        phase.set_collision_integral_polynomial(kO2, kO2, *coll_polys_O2, actualT=False)
        assert get_cond("O2") == cond1_O2  # back to original


class TestIonTransport:

    @pytest.fixture(scope='function')
    def gas(self):
        gas = ct.Solution('ch4_ion.yaml')
        gas.TPX = 2237, ct.one_atm, 'O2:0.7010, H2O:0.1885, CO2:9.558e-2'
        return gas

    def test_binary_diffusion(self, gas):
        N2_idx = gas.species_index("N2")
        H3Op_idx = gas.species_index("H3O+")

        bdiff = gas.binary_diff_coeffs[N2_idx][H3Op_idx]
        assert bdiff == approx(4.258e-4, rel=1e-4)  # Regression test

    def test_mixture_diffusion(self, gas):
        H3Op_idx = gas.species_index("H3O+")

        mdiff = gas.mix_diff_coeffs[H3Op_idx]
        assert mdiff == approx(5.057e-4, rel=1e-4)  # Regression test

    def test_O2_anion_mixture_diffusion(self, gas):
        mdiff = gas['O2-'].mix_diff_coeffs[0]
        assert mdiff == approx(2.784e-4, rel=1e-3)  # Regression test

    def test_mobility(self, gas):
        H3Op_idx = gas.species_index("H3O+")

        mobi = gas.mobilities[H3Op_idx]
        assert mobi == approx(2.623e-3, rel=1e-4)  # Regression test

    def test_update_temperature(self, gas):
        N2_idx = gas.species_index("N2")
        H3Op_idx = gas.species_index("H3O+")

        bdiff = gas.binary_diff_coeffs[N2_idx][H3Op_idx]
        mdiff = gas.mix_diff_coeffs[H3Op_idx]
        mobi = gas.mobilities[H3Op_idx]
        gas.TP = 0.9 * gas.T, gas.P
        assert bdiff != gas.binary_diff_coeffs[N2_idx][H3Op_idx]
        assert mdiff != gas.mix_diff_coeffs[H3Op_idx]
        assert mobi != gas.mobilities[H3Op_idx]


@pytest.mark.parametrize(
    "key,value,message",
    [
        ("H_geom", "linear", "invalid geometry"),
        ("H_geom", "nonlinear", "invalid geometry"),
        ("H2_geom", "atom", "invalid geometry"),
        ("H2_geom", "nonsense", "invalid geometry"),
        ("H2_geom", "nonlinear", "invalid geometry"),
        ("H2O_geom", "atom", "invalid geometry"),
        ("OHp_geom", "atom", "invalid geometry"),
        ("OHp_geom", "nonlinear", "invalid geometry"),
        ("E_geom", "linear", "invalid geometry"),
        ("H2_well", -33.4, "negative well depth.*H2"),
        ("H2O_diam", 0.0, "negative or zero diameter.*H2O"),
        ("H2O_dipole", -1.84, "negative dipole moment.*H2O"),
        ("H2_polar", -0.79, "negative polarizability.*H2"),
        ("OHp_rot", -4, "negative rotation relaxation number.*OHp"),
        ("H2_disp", -3.1, "negative dispersion coefficient.*H2"),
        ("H2_quad", -3.1, "negative quadrupole polarizability.*H2"),
    ]
)
def test_bad_transport_input(key, value, message):
    """ Check that invalid transport inputs raise appropriate exceptions """
    # Default parameters are valid
    subs = {"H_geom":"atom", "H2_geom":"linear", "H2O_geom":"nonlinear",
            "OHp_geom":"linear", "E_geom":"atom", "H2_well": 38.0, "H2O_diam": 2.60,
            "H2O_dipole": 1.84, "H2_polar": 0.79, "OHp_rot": 4.0, "H2_disp": 2.995,
            "H2_quad": 3.602}
    subs[key] = value

    species_data = """
    - name: H2
      composition: {{H: 2}}
      thermo: &dummy-thermo
        {{model: constant-cp, T0: 1000, h0: 51.7, s0: 19.5, cp0: 8.41}}
      transport:
        model: gas
        geometry: {H2_geom}
        diameter: 2.92
        well-depth: {H2_well}
        polarizability: {H2_polar}
        rotational-relaxation: 280.0
        dispersion-coefficient: {H2_disp}
        quadrupole-polarizability: {H2_quad}
    - name: H
      composition: {{H: 1}}
      thermo: *dummy-thermo
      transport:
        model: gas
        geometry: {H_geom}
        diameter: 2.05
        well-depth: 145.00
    - name: H2O
      composition: {{H: 2, O: 1}}
      thermo: *dummy-thermo
      transport:
        model: gas
        geometry: {H2O_geom}
        diameter: {H2O_diam}
        well-depth: 572.40
        dipole: {H2O_dipole}
        rotational-relaxation: 4.0
    - name: OHp
      composition: {{H: 1, O: 1, E: -1}}
      thermo: *dummy-thermo
      transport:
        model: gas
        geometry: {OHp_geom}
        diameter: 2.60
        well-depth: 572.40
        dipole: 1.84
        rotational-relaxation: {OHp_rot}
    - name: E
      composition: {{E: 1}}
      thermo: *dummy-thermo
      transport:
        model: gas
        geometry: {E_geom}
        diameter: 0.01
        well-depth: 1.0
    """.format(**subs)

    with pytest.raises(ct.CanteraError, match=message):
        ct.Species.list_from_yaml(species_data)

@pytest.mark.parametrize(
    "model",
    [
        "mixture-averaged",
        "mixture-averaged-CK",
        "ionized-gas",
        pytest.param("multicomponent",
                     marks=pytest.mark.xfail(reason="See Issue #1823"))
     ]
)
def test_single_species_transport(model):
    """
    A phase with only one species defined should have the same transport
    properties as a pure species state in a multi-species phase definition.
    """
    yaml_ref = """
    phases:
    - name: gas
      thermo: ideal-gas
      species:
      - gri30.yaml/species: [H2, O2, H2O2, OH]
    """
    ref = ct.Solution(yaml=yaml_ref, transport_model=model)
    single = ct.Solution(thermo='ideal-gas', species=[ref.species('H2O2')],
                         transport_model=model)
    single.TPX = ref.TPX = 500, 5 * ct.one_atm, 'H2O2:1.0'
    assert single.min_temp == ref.min_temp
    assert single.max_temp == ref.max_temp
    # assert single.viscosity == approx(ref.viscosity)
    assert single.thermal_conductivity == approx(ref.thermal_conductivity)
    k = ref.species_index('H2O2')
    assert single.mix_diff_coeffs[0] == approx(ref.binary_diff_coeffs[k,k])

class TestDustyGas:

    @pytest.fixture
    def phase(self):
        phase = ct.DustyGas("h2o2.yaml")
        phase.TPX = 500.0, ct.one_atm, "O2:2.0, H2:1.0, H2O:1.0"
        phase.porosity = 0.2
        phase.tortuosity = 0.3
        phase.mean_pore_radius = 1e-4
        phase.mean_particle_diameter = 5e-4
        return phase

    @pytest.fixture
    def Dref(self, phase):
        return phase.multi_diff_coeffs

    def test_porosity(self, phase, Dref):
        phase.porosity = 0.4
        D = phase.multi_diff_coeffs
        assert Dref * 2 == approx(D)

    def test_tortuosity(self, phase, Dref):
        phase.tortuosity = 0.6
        D = phase.multi_diff_coeffs
        assert Dref * 0.5 == approx(D)

    # The other parameters don't have such simple relationships to the diffusion
    # coefficients, so we can't test them as easily

    def test_molar_fluxes(self, phase):
        T1, rho1, Y1 = phase.TDY
        phase.TPX = 500.0, ct.one_atm, "O2:2.0, H2:1.001, H2O:0.999"

        T2, rho2, Y2 = phase.TDY

        fluxes0 = phase.molar_fluxes(T1, T1, rho1, rho1, Y1, Y1, 1e-4)
        assert fluxes0 == approx(np.zeros(phase.n_species))

        fluxes1 = phase.molar_fluxes(T1, T2, rho1, rho2, Y1, Y2, 1e-4)
        H2_idx = phase.species_index('H2')
        H2O_idx = phase.species_index('H2O')
        assert fluxes1[H2_idx] < 0
        assert fluxes1[H2O_idx] > 0

        # Not sure why the following condition is not satisfied:
        #assert sum(fluxes1) == approx(0.0)
        #assert sum(fluxes1) / sum(abs(fluxes1)) == approx(0.0)

    def test_thermal_conductivity(self, phase):
        gas1 = ct.Solution("h2o2.yaml", transport_model="multicomponent")
        gas1.TPX = phase.TPX

        assert phase.thermal_conductivity == gas1.thermal_conductivity



class TestWaterTransport:
    """
    Comparison values are taken from the NIST Chemistry WebBook. Agreement is
    limited by the use of a different equation of state here (Reynolds) than
    in the Webbook (IAPWS95), as well as a different set of coefficients for
    the transport property model. Differences are largest in the region near
    the critical point.
    """

    @pytest.fixture(scope='class')
    def water(self):
        return ct.Water()

    @pytest.mark.parametrize("T, P, mu, rtol", [
        (400, 1e6, 2.1880e-4, 1e-3),
        (400, 8e6, 2.2061e-4, 1e-3),
        (620, 1.6e7, 6.7489e-5, 2e-3),
        (620, 2.8e7, 7.5684e-5, 2e-3),
    ])
    def test_viscosity_liquid(self, water, T, P, mu, rtol):
        water.TP = T, P
        assert water.viscosity == approx(mu, rel=rtol)

    @pytest.mark.parametrize("T, P, mu, rtol", [
        (600, 1e6, 2.1329e-5, 1e-3),
        (620, 5e6, 2.1983e-5, 1e-3),
        (620, 1.5e7, 2.2858e-5, 2e-3),
    ])
    def test_viscosity_vapor(self, water, T, P, mu, rtol):
        water.TP = T, P
        assert water.viscosity == approx(mu, rel=rtol)

    @pytest.mark.parametrize("T, P, mu, rtol", [
        (660, 2.2e7, 2.7129e-5, 2e-3),
        (660, 2.54e7, 3.8212e-5, 1e-2),
        (660, 2.8e7, 5.3159e-5, 1e-2),
    ])
    def test_viscosity_supercritical(self, water, T, P, mu, rtol):
        water.TP = T, P
        assert water.viscosity == approx(mu, rel=rtol)

    @pytest.mark.parametrize("T, D, mu, rtol", [
        (647.43, 280.34, 3.7254e-05, 6e-3),
        (647.43, 318.89, 4.2286e-05, 6e-3),
        (648.23, 301.34, 3.9136e-05, 6e-3),
        (648.23, 330.59, 4.2102e-05, 6e-3)
    ])
    def test_viscosity_near_critical(self, water, T, D, mu, rtol):
        water.TD = T, D
        assert water.viscosity == approx(mu, rel=rtol)

    @pytest.mark.parametrize("T, P, k, rtol", [
        (400, 1e6, 0.68410, 1e-3),
        (400, 8e6, 0.68836, 1e-3),
        (620, 1.6e7, 0.45458, 2e-3),
        (620, 2.8e7, 0.49705, 2e-3),
    ])
    def test_thermal_conductivity_liquid(self, water, T, P, k, rtol):
        water.TP = T, P
        assert water.thermal_conductivity == approx(k, rel=rtol)

    @pytest.mark.parametrize("T, P, k, rtol", [
        (600, 1e6, 0.047636, 1e-3),
        (620, 5e6, 0.055781, 1e-3),
        (620, 1.5e7, 0.10524, 2e-3),
    ])
    def test_thermal_conductivity_vapor(self, water, T, P, k, rtol):
        water.TP = T, P
        assert water.thermal_conductivity == approx(k, rel=rtol)

    @pytest.mark.parametrize("T, P, k, rtol", [
        (660, 2.2e7, 0.14872, 1e-2),
        (660, 2.54e7, 0.35484, 2e-2),
        (660, 2.8e7, 0.38479, 1e-2),
    ])
    def test_thermal_conductivity_supercritical(self, water, T, P, k, rtol):
        water.TP = T, P
        assert water.thermal_conductivity == approx(k, rel=rtol)


class TestIAPWS95WaterTransport:
    """
    Water transport properties test using the IAPWS95 equation of state. This
    results in better comparisons with data from the NIST Webbook.
    """

    @pytest.fixture(scope='class')
    def water(self):
        return ct.Solution('thermo-models.yaml', 'liquid-water')

    @pytest.mark.parametrize("T, P, mu, rtol", [
        (400, 1e6, 2.1880e-4, 2e-4),
        (400, 8e6, 2.2061e-4, 2e-4),
        (620, 1.6e7, 6.7489e-5, 1e-4),
        (620, 2.8e7, 7.5684e-5, 1e-4),
    ])
    def test_viscosity_liquid(self, water, T, P, mu, rtol):
        water.TP = T, P
        assert water.viscosity == approx(mu, rel=rtol)

    @pytest.mark.parametrize("T, P, mu, rtol", [
        (660, 2.2e7, 2.7129e-5, 1e-4),
        (660, 2.54e7, 3.8212e-5, 1e-4),
        (660, 2.8e7, 5.3159e-5, 1e-4),
    ])
    def test_viscosity_supercritical(self, water, T, P, mu, rtol):
        water.TP = T, P
        assert water.viscosity == approx(mu, rel=rtol)

    @pytest.mark.parametrize("T, P, k, rtol", [
        (400, 1e6, 0.68410, 1e-4),
        (400, 8e6, 0.68836, 1e-4),
        (620, 1.6e7, 0.45458, 1e-4),
        (620, 2.8e7, 0.49705, 1e-4),
    ])
    def test_thermal_conductivity_liquid(self, water, T, P, k, rtol):
        water.TP = T, P
        assert water.thermal_conductivity == approx(k, rel=rtol)

    @pytest.mark.parametrize("T, P, k, rtol", [
        (660, 2.2e7, 0.14872, 1e-4),
        (660, 2.54e7, 0.35484, 1e-4),
        (660, 2.8e7, 0.38479, 1e-4),
    ])
    def test_thermal_conductivity_supercritical(self, water, T, P, k, rtol):
        water.TP = T, P
        assert water.thermal_conductivity == approx(k, rel=rtol)


class TestTransportData:

    @pytest.fixture(scope='class')
    def gas(self):
        gas = ct.Solution("h2o2.yaml")
        return gas

    def test_read(self, gas):
        tr = gas.species('H2').transport
        assert tr.geometry == 'linear'
        assert tr.diameter == approx(2.92e-10)
        assert tr.well_depth == approx(38.0 * ct.boltzmann)
        assert tr.polarizability == approx(0.79e-30)
        assert tr.rotational_relaxation == approx(280)

    def test_set_customary_units(self, gas):
        tr1 = ct.GasTransportData()
        tr1.set_customary_units('nonlinear', 2.60, 572.40, 1.84, 0.0, 4.00)
        tr2 = gas.species('H2O').transport
        assert tr1.geometry == tr2.geometry
        assert tr1.diameter == approx(tr2.diameter)
        assert tr1.well_depth == approx(tr2.well_depth)
        assert tr1.dipole == approx(tr2.dipole)
        assert tr1.polarizability == approx(tr2.polarizability)
        assert tr1.rotational_relaxation == approx(tr2.rotational_relaxation)


class TestIonGasTransportData:
    @pytest.fixture(scope='class')
    def gas(self):
        return ct.Solution("ch4_ion.yaml")

    def test_read_ion(self, gas):
        tr = gas.species('N2').transport
        assert tr.dispersion_coefficient == approx(2.995e-50)
        assert tr.quadrupole_polarizability == approx(3.602e-50)

    def test_set_customary_units(self, gas):
        tr1 = ct.GasTransportData()
        tr1.set_customary_units('linear', 3.62, 97.53, 1.76,
                                dispersion_coefficient = 2.995,
                                quadrupole_polarizability = 3.602)
        tr2 = gas.species('N2').transport
        assert tr1.dispersion_coefficient == approx(tr2.dispersion_coefficient)
        assert tr1.quadrupole_polarizability == approx(tr2.quadrupole_polarizability)

    def test_serialization(self, gas):
        data = gas.species('N2').transport.input_data
        assert data['dispersion-coefficient'] == approx(2.995)
        assert data['quadrupole-polarizability'] == approx(3.602)

        assert 'dispersion-coefficient' not in gas.species('CO2').transport.input_data


class TestHighPressureGasTransport:
    """
    Note: to re-create the reference file:
    (1) Set PYTHONPATH to build/python.
    (2) Go into test/python directory and run:
        pytest --save-reference=high_pressure_transport test_transport.py::TestHighPressureGasTransport::test_high_pressure_transport
        pytest --save-reference=high_pressure_chung_transport test_transport.py::TestHighPressureGasTransport::test_high_pressure_chung_transport
    (3) Compare the reference files created in the current working directory with
        the ones in test/data and replace them if needed.
    """

    def test_failure_for_species_with_no_properties(self):
        """
        All species must have critical properties defined to use the high pressure
        transport model. This test uses a YAML file with a specified value of the
        a and b parameters for the Peng-Robinson equation of state, which should
        be parsed by the thermo model and will set the critical properties to
        non-physical values. These non-physical values should trigger an error
        in the high-pressure transport model.
        """
        gas = ct.Solution('methane_co2_noCritProp.yaml')

        with pytest.raises(ct.CanteraError, match="must have critical properties defined"):
            gas.transport_model = 'high-pressure-Chung'

        with pytest.raises(ct.CanteraError, match="must have critical properties defined"):
            gas.transport_model = 'high-pressure'

    def test_high_pressure_transport(self, request, test_data_path):
        """
        This test compares the viscosity and thermal conductivities of a mixture of
        CH4 and CO2 over a range of pressures using the high-pressure transport
        model and compares the results to a reference file. This is a regression test.
        """
        referenceFile = "HighPressureTest.csv"

        phasedef = """
        phases:
        - name: methane_co2
          species:
          - gri30.yaml/species: [CH4, CO2]
          thermo: Peng-Robinson
          transport: mixture-averaged
          state: {T: 300, P: 1 atm}
        """
        gas = ct.Solution(yaml=phasedef)
        gas.transport_model = 'high-pressure'

        pressures = np.linspace(101325, 6e7, 100)
        # Collect viscosities and thermal conductivities
        viscosities = []
        thermal_conductivities = []
        diffusion_coefficients = []
        for pressure in pressures:
            gas.TPX = 350, pressure, 'CH4:0.755, CO2:0.245'
            viscosities.append(gas.viscosity)
            thermal_conductivities.append(gas.thermal_conductivity)
            diffusion_coefficients.append(gas.mix_diff_coeffs)

        data = np.empty((len(viscosities), 3+ len(diffusion_coefficients[0])))
        data[:,0] = pressures
        data[:,1] = viscosities
        data[:,2] = thermal_conductivities
        for i in range(len(diffusion_coefficients[0])):
            data[:,3+i] = [d[i] for d in diffusion_coefficients]

        saveReference = request.config.getoption("--save-reference")
        if saveReference == 'high_pressure_transport':
            np.savetxt(referenceFile, data, '%11.6e', ', ')
        else:
            bad = compareProfiles(test_data_path / referenceFile, data,
                                  rtol=1e-2, atol=1e-8, xtol=1e-2)
            assert not bad, bad

    def test_high_pressure_chung_transport(self, request, test_data_path):
        referenceFile = "HighPressureChungTest.csv"

        phasedef = """
        phases:
        - name: methane_co2
          species:
          - gri30.yaml/species: [CH4, CO2]
          thermo: Peng-Robinson
          transport: mixture-averaged
          state: {T: 300, P: 1 atm}
        """
        gas = ct.Solution(yaml=phasedef)
        gas.transport_model = 'high-pressure-Chung'

        pressures = np.linspace(101325, 6e7, 100)
        # Collect viscosities and thermal conductivities
        viscosities = []
        thermal_conductivities = []
        diffusion_coefficients = []
        for pressure in pressures:
            gas.TPX = 350, pressure, 'CH4:0.755, CO2:0.245'
            viscosities.append(gas.viscosity)
            thermal_conductivities.append(gas.thermal_conductivity)
            diffusion_coefficients.append(gas.mix_diff_coeffs)

        data = np.empty((len(viscosities), 3+gas.n_species))
        data[:,0] = pressures
        data[:,1] = viscosities
        data[:,2] = thermal_conductivities
        for i in range(len(diffusion_coefficients[0])):
            data[:,3+i] = [d[i] for d in diffusion_coefficients]

        saveReference = request.config.getoption("--save-reference")
        if saveReference == 'high_pressure_chung_transport':
            np.savetxt(referenceFile, data, '%11.6e', ', ')
        else:
            bad = compareProfiles(test_data_path / referenceFile, data,
                                  rtol=1e-2, atol=1e-8, xtol=1e-2)
            assert not bad, bad
