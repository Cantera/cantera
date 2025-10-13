import numpy as np
import pytest
from pytest import approx
import re
pd = pytest.importorskip("pandas")

import cantera as ct

from .utilities import yaml
from .utilities import (
    compareProfiles
)


class TestOnedim:
    def test_instantiate(self):
        gas = ct.Solution("h2o2.yaml")
        free = ct.FreeFlow(gas)
        axi = ct.AxisymmetricFlow(gas)

    def test_badInstantiate(self):
        solid = ct.Solution("diamond.yaml", "diamond")
        with pytest.raises(ct.CanteraError,
                           match="An appropriate transport model\nshould be set when instantiating"):
            ct.FreeFlow(solid)

    def test_instantiateSurface(self):
        gas = ct.Solution("diamond.yaml", "gas")
        solid = ct.Solution("diamond.yaml", "diamond")
        interface = ct.Solution("diamond.yaml", "diamond_100", (gas, solid))

        surface = ct.ReactingSurface1D(phase=interface)

    def test_boundaryProperties(self):
        gas1 = ct.Solution("h2o2.yaml")
        gas2 = ct.Solution("h2o2.yaml")
        inlet = ct.Inlet1D(name='something', phase=gas1)
        flame = ct.FreeFlow(gas1)
        sim = ct.Sim1D((inlet, flame))

        assert inlet.name == 'something'

        gas2.TPX = 400, 101325, 'H2:0.3, O2:0.5, AR:0.2'
        Xref = gas2.X
        Yref = gas2.Y
        inlet.Y = Yref

        assert inlet.Y == approx(Yref)
        assert inlet.X == approx(Xref)

        gas2.TPX = 400, 101325, 'H2:0.5, O2:0.2, AR:0.3'
        Xref = gas2.X
        Yref = gas2.Y
        inlet.X = Xref
        assert inlet.X == approx(Xref)
        assert inlet.Y == approx(Yref)

        inlet.X = {'H2':0.3, 'O2':0.5, 'AR':0.2}
        assert inlet.X[gas2.species_index('H2')] == approx(0.3)

    def test_grid_check(self):
        gas = ct.Solution("h2o2.yaml")
        flame = ct.FreeFlow(gas)

        with pytest.raises(ct.CanteraError, match='monotonically'):
            flame.grid = [0, 0.1, 0.1, 0.2]

        with pytest.raises(ct.CanteraError, match='monotonically'):
            flame.grid = [0, 0.1, 0.2, 0.05]

    def test_unpicklable(self):
        import pickle
        gas = ct.Solution("h2o2.yaml")
        flame = ct.FreeFlow(gas)
        with pytest.raises(NotImplementedError):
            pickle.dumps(flame)

    def test_uncopyable(self):
        import copy
        gas = ct.Solution("h2o2.yaml")
        flame = ct.FreeFlow(gas)
        with pytest.raises(NotImplementedError):
            copy.copy(flame)

    def test_exceptions(self):
        with pytest.raises(TypeError, match="Argument 'phase' has incorrect type"):
            ct.Inlet1D(None)
        gas = ct.Solution("h2o2.yaml")
        with pytest.raises(ct.CanteraError, match="incompatible ThermoPhase type"):
            ct.ReactingSurface1D(gas)
        interface = ct.Solution("diamond.yaml", "diamond_100")
        with pytest.raises(TypeError, match="unexpected keyword"):
            ct.ReactingSurface1D(interface, foo="bar")
        surf = ct.ReactingSurface1D(interface)

    def test_invalid_property(self):
        gas1 = ct.Solution("h2o2.yaml")
        inlet = ct.Inlet1D(name='something', phase=gas1)
        flame = ct.FreeFlow(gas1)
        sim = ct.Sim1D((inlet, flame))

        for x in (inlet, flame, sim):
            with pytest.raises(AttributeError):
                x.foobar = 300
            with pytest.raises(AttributeError):
                x.foobar

    def test_tolerances(self):
        gas = ct.Solution("h2o2.yaml")
        left = ct.Inlet1D(gas)
        flame = ct.FreeFlow(gas)
        right = ct.Outlet1D(gas)
        # Some things don't work until the domains have been added to a Sim1D
        sim = ct.Sim1D((left, flame, right))

        with pytest.raises(ct.CanteraError, match="Component 'foobar' not found"):
            flame.set_steady_tolerances(foobar=(3e-4, 3e-6))

        flame.set_steady_tolerances(default=(5e-3, 5e-5),
                                    T=(3e-4, 3e-6),
                                    Y=(7e-7, 7e-9))
        flame.set_transient_tolerances(default=(6e-3, 6e-5),
                                       T=(4e-4, 4e-6),
                                       Y=(2e-7, 2e-9))

        # Flow domain
        atol_ss = set(flame.steady_abstol())
        atol_ts = set(flame.transient_abstol())
        rtol_ss = set(flame.steady_reltol())
        rtol_ts = set(flame.transient_reltol())

        assert atol_ss == set((5e-5, 3e-6, 7e-9))
        assert atol_ts == set((6e-5, 4e-6, 2e-9))
        assert rtol_ss == set((5e-3, 3e-4, 7e-7))
        assert rtol_ts == set((6e-3, 4e-4, 2e-7))

    def test_switch_transport(self):
        gas = ct.Solution('h2o2.yaml')
        gas.set_equivalence_ratio(0.9, 'H2', 'O2:0.21, N2:0.79')
        flame = ct.FreeFlame(gas, width=0.1)
        flame.set_initial_guess()

        assert gas.transport_model == flame.transport_model == 'mixture-averaged'

        flame.transport_model = 'unity-Lewis-number'
        assert gas.transport_model == flame.transport_model == 'unity-Lewis-number'
        Dkm_flame = flame.mix_diff_coeffs
        assert all(Dkm_flame[1,:] == Dkm_flame[2,:])

        gas.transport_model = 'multicomponent'
        assert flame.transport_model == 'multicomponent'

        with pytest.raises(ct.CanteraError, match="Invalid Transport model"):
            flame.transport_model = 'none'

        gas.transport_model = 'unity-Lewis-number'
        with pytest.raises(ct.CanteraError, match="Invalid Transport model"):
            gas.transport_model = 'none'

    def test_width_grid(self):
        gas = ct.Solution('h2o2.yaml')
        for cls in ct.FlameBase.__subclasses__():
            with pytest.raises(ValueError, match="mutually exclusive"):
                sim = cls(gas, grid=[0, 0.1, 0.2], width=0.4)

def check_component_order(fname: str, group: str):
    with fname.open("r", encoding="utf-8") as fid:
        reader = yaml.YAML(typ="safe")
        contents = reader.load(fid)
    components = []
    index = -1
    for key, value in contents[group]["flame"].items():
        if components:
            index += 1
            field = components.pop()
            msg = f"Found mismatch of components at index {index}: expected {field!r}"
            assert key == field, msg
        if key == "components":
            components = value[::-1]


class TestFreeFlame:
    tol_ss = [1.0e-5, 1.0e-14]  # [rtol atol] for steady-state problem
    tol_ts = [1.0e-4, 1.0e-11]  # [rtol atol] for time stepping

    def create_sim(self, p, Tin, reactants, width=0.05, mech="h2o2.yaml"):
        # Solution object used to compute mixture properties
        self.gas = ct.Solution(mech)
        self.gas.TPX = Tin, p, reactants

        # Flame object
        self.sim = ct.FreeFlame(self.gas, width=width)
        self.sim.flame.set_steady_tolerances(default=self.tol_ss)
        self.sim.flame.set_transient_tolerances(default=self.tol_ts)

        # Set properties of the upstream fuel-air mixture
        self.sim.inlet.T = Tin
        self.sim.inlet.X = reactants

    def solve_fixed_T(self, loglevel=0):
        # Solve with the energy equation disabled
        self.sim.energy_enabled = False
        self.sim.solve(loglevel=loglevel, refine_grid=False)

        assert not self.sim.energy_enabled

    def solve_mix(self, ratio=3.0, slope=0.3, curve=0.2, prune=0.0, refine=True,
                  auto=False):
        # Solve with the energy equation enabled
        self.sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)
        self.sim.energy_enabled = True
        self.sim.solve(loglevel=0, refine_grid=refine, auto=auto)

        assert self.sim.energy_enabled
        assert self.sim.transport_model == 'mixture-averaged'

    def solve_multi(self):
        self.sim.transport_model = 'multicomponent'
        self.sim.solve(loglevel=0, refine_grid=True)

        assert self.sim.transport_model == 'multicomponent'

        data = self.test_work_path / f"multicomponent.yaml"
        data.unlink(missing_ok=True)
        group = "multicomponent"
        self.sim.save(data, group)

        arr = ct.SolutionArray(self.sim.gas)
        arr.restore(data, "multicomponent/flame")
        assert arr.transport_model == "multicomponent"

    def test_flow_type(self):
        Tin = 300
        p = ct.one_atm
        reactants = 'H2:0.65, O2:0.5, AR:2'
        self.create_sim(p, Tin, reactants, width=0.0001)

        assert self.sim.flame.domain_type == "free-flow"

    def test_fixed_temperature(self):
        # test setting of fixed temperature
        Tin = 300
        p = ct.one_atm
        reactants = 'H2:0.65, O2:0.5, AR:2'
        self.create_sim(p, Tin, reactants, width=0.0001)
        self.sim.set_initial_guess()
        zvec = self.sim.grid
        tvec = self.sim.T
        tfixed = 900.
        self.sim.fixed_temperature = tfixed
        zfixed = np.interp(tfixed, tvec, zvec)
        assert self.sim.fixed_temperature == approx(tfixed)
        assert self.sim.fixed_temperature_location == approx(zfixed)

        # test some alias names
        assert np.allclose(self.sim.T, self.sim.flame.temperature)
        zeros = self.sim.T * 0
        with pytest.raises(ct.CanteraError, match="is not used by 'free-flow'"):
            self.sim.flame.set_values("spread_rate", zeros)
        with pytest.raises(ct.CanteraError, match="is not used by 'free-flow'"):
            self.sim.flame.set_values("spreadRate", zeros)
        with pytest.raises(ct.CanteraError, match="is not used by 'free-flow'"):
            self.sim.flame.set_values("radial-pressure-gradient", zeros)
        with pytest.raises(ct.CanteraError, match="is not used by 'free-flow'"):
            self.sim.flame.set_values("L", zeros)
        with pytest.raises(ct.CanteraError, match="is not used by 'free-flow'"):
            self.sim.flame.set_values("Lambda", zeros)
        with pytest.raises(ct.CanteraError, match="is not used by 'free-flow'"):
            self.sim.flame.set_values("electric-field", zeros)
        with pytest.raises(ct.CanteraError, match="is not used by 'free-flow'"):
            self.sim.flame.set_values("E", zeros)
        with pytest.raises(ct.CanteraError, match="is not used by 'free-flow'"):
            self.sim.flame.set_values("eField", zeros)
        with pytest.raises(ct.CanteraError, match="is not used by 'free-flow'"):
            self.sim.flame.set_values("oxidizer-velocity", zeros)
        with pytest.raises(ct.CanteraError, match="is not used by 'free-flow'"):
            self.sim.flame.set_values("Uo", zeros)

    @pytest.mark.slow_test
    def test_auto_width(self):
        Tin = 300
        p = ct.one_atm
        reactants = 'H2:0.65, O2:0.5, AR:2'
        self.create_sim(p, Tin, reactants, width=0.0001)
        self.sim.set_refine_criteria(ratio=3, slope=0.3, curve=0.2)
        self.sim.solve(loglevel=0, refine_grid=True, auto=True)

        self.gas.TPX = Tin, p, reactants
        self.gas.equilibrate('HP')
        Tad = self.gas.T
        assert Tad == approx(self.sim.T[-1], rel=2e-2)

        # Re-solving with auto=False should not trigger a DomainTooNarrow
        # exception, and should leave domain width constant
        self.sim.flame.grid *= 0.3
        old_width = self.sim.grid[-1]
        self.sim.solve(loglevel=0, refine_grid=True, auto=False)
        assert self.sim.grid[-1] == approx(old_width)

    def test_auto_width2(self):
        self.create_sim(p=ct.one_atm, Tin=400, reactants='H2:0.8, O2:0.5',
                        width=0.1)

        self.sim.set_refine_criteria(ratio=4, slope=0.8, curve=0.8)
        self.sim.flame.set_steady_tolerances(T=(2e-4, 1e-8))
        self.sim.solve(refine_grid=True, auto=True, loglevel=0)

        assert self.sim.velocity[0] == approx(17.02, rel=1e-1)
        assert self.sim.grid[-1] < 2.0 # grid should not be too wide
        assert self.sim.flame.tolerances("T") == (2e-4, 1e-8)

    @pytest.mark.slow_test
    def test_converge_adiabatic(self):
        # Test that the adiabatic flame temperature and species profiles
        # converge to the correct equilibrium values as the grid is refined

        reactants = 'H2:1.1, O2:1, AR:5'
        p = ct.one_atm
        Tin = 300

        self.create_sim(p, Tin, reactants)
        self.solve_fixed_T()

        self.gas.TPX = Tin, p, reactants
        self.gas.equilibrate('HP')
        Tad = self.gas.T
        Xad = self.gas.X

        self.solve_mix(slope=0.5, curve=0.3)
        T1 = self.sim.T[-1]
        X1 = self.sim.X[:,-1]

        self.solve_mix(slope=0.2, curve=0.1)
        T2 = self.sim.T[-1]
        X2 = self.sim.X[:,-1]

        self.solve_mix(slope=0.1, curve=0.05)
        T3 = self.sim.T[-1]
        X3 = self.sim.X[:,-1]

        assert abs(T2-Tad) < abs(T1-Tad)
        assert abs(T3-Tad) < abs(T2-Tad)

        for k in range(self.gas.n_species):
            if Xad[k] > self.tol_ss[1]:
                assert abs(X2[k]-Xad[k]) <= abs(X1[k]-Xad[k])
                assert abs(X3[k]-Xad[k]) <= abs(X2[k]-Xad[k])

    def run_mix(self, phi, T, width, p, refine):
        reactants = {'H2': phi, 'O2': 0.5, 'AR': 2}
        self.create_sim(p * ct.one_atm, T, reactants, width)
        self.solve_mix(refine=refine)

        rhou = self.sim.inlet.mdot

        # Check continuity
        for rhou_j in self.sim.density * self.sim.velocity:
            assert rhou_j == approx(rhou, rel=1e-4)

    def test_solution_array_output(self):
        self.run_mix(phi=1.0, T=300, width=2.0, p=1.0, refine=False)

        flow = self.sim.to_array(normalize=True)
        assert self.sim.grid == approx(flow.grid)
        assert self.sim.T == approx(flow.T)

        f2 = ct.FreeFlame(self.gas)
        f2.from_array(flow)
        assert self.sim.grid == approx(f2.grid)
        assert self.sim.T == approx(f2.T)

        inlet = self.sim.to_array(self.sim.inlet)
        f2.from_array(inlet, f2.inlet)
        assert self.sim.inlet.T == approx(f2.inlet.T)
        assert self.sim.inlet.mdot == approx(f2.inlet.mdot)
        assert self.sim.inlet.Y == approx(f2.inlet.Y)

    def test_restart_array(self):
        # restart from SolutionArray
        self.run_restart("array")

    def test_restart_csv(self):
        # restart from CSV format
        self.run_restart("csv")

    def test_restart_yaml(self):
        # restart from YAML format
        filename, group = self.run_restart("yaml")
        check_component_order(filename, group)

    @pytest.mark.skipif("native" not in ct.hdf_support(),
                        reason="Cantera compiled without HDF support")
    def test_restart_hdf(self):
        # restart from HDF format
        self.run_restart("h5")

    @pytest.mark.skipif(not pd, reason="Pandas not installed")
    def test_restart_pandas(self):
        # restart from Pandas DataFrame
        self.run_restart("pandas", auto=True)

    def run_restart(self, mode, auto=False):
        self.run_mix(phi=1.0, T=300, width=2.0, p=1.0, refine=False)

        group = "restart"
        if mode == "array":
            data = self.sim.to_array()
        elif mode == "pandas":
            data = self.sim.to_pandas()
        else:
            data = self.test_work_path / f"freeflame_restart.{mode}"
            data.unlink(missing_ok=True)
            if mode == "csv":
                self.sim.save(data, basis="mole")
            else:
                self.sim.save(data, group)

        reactants = {'H2': 0.9, 'O2': 0.5, 'AR': 2}
        self.create_sim(1.1 * ct.one_atm, 500, reactants, 2.0)
        self.sim.set_initial_guess(data=data, group=group)
        self.solve_mix(refine=False, auto=auto)

        # Check continuity
        rhou = self.sim.inlet.mdot
        for rhou_j in self.sim.density * self.sim.velocity:
            assert rhou_j == approx(rhou, rel=1e-4)
        return data, group

    def test_settings(self):
        self.create_sim(p=ct.one_atm, Tin=400, reactants='H2:0.8, O2:0.5', width=0.1)
        self.sim.set_initial_guess()
        sim = self.sim

        # new implementation
        new_keys = {
            "type", "points", "tolerances", "transport-model", "phase",
            "radiation-enabled", "energy-enabled", "Soret-enabled",
            "refine-criteria", "fixed-point"}
        settings = sim.flame.settings
        for k in new_keys:
            assert k in settings

        assert settings["radiation-enabled"] == False

        # Apply settings using individual setters
        sim.flame.boundary_emissivities = 0.12, 0.34
        sim.flame.radiation_enabled = True
        sim.flame.set_steady_tolerances(default=(3e-4, 7e-9))
        sim.transport_model = 'unity-Lewis-number'

        # Check that the aggregated settings reflect the changes
        new_settings = sim.flame.settings

        assert new_settings["radiation-enabled"] == True
        assert new_settings["emissivity-left"] == 0.12
        assert new_settings["emissivity-right"] == 0.34
        assert new_settings["transport-model"] == "unity-Lewis-number"

        tolerances = new_settings["tolerances"]
        assert tolerances["steady-reltol"] == approx(3e-4)
        assert tolerances["steady-abstol"] == approx(7e-9)

        assert "fixed-point" in new_settings
        assert "location" in new_settings["fixed-point"]

    def test_log_solution_bounds(self, capsys):
        self.create_sim(p=ct.one_atm, Tin=300, reactants="H2:1.1, O2:1, AR:5")
        self.solve_fixed_T(loglevel=4)

        out = capsys.readouterr().out
        assert "Undamped Newton step takes solution out of bounds" in out
        m = re.search(r"     ={60,}\n(.*?)\n     ={60,}\n(.*?)\n     ={60,}", out,
                      re.DOTALL)
        assert m is not None

        header = m.group(1).splitlines()
        body = m.group(2).splitlines()
        assert re.match(" +Value +Min +Max +$", header[0])
        assert re.match(" +Domain +Pt. Component +Value +Change +Bound +Bound +$",
                        header[1])
        assert len(body) > 0
        for line in body:
            assert re.match(r" +flame +\d+ +[A-Z0-9]+ +([0-9.e+\-]+ *){4}$", line)

    def test_timestep_jacobian_limits(self, capsys):
        # Use log output to confirm that maximum timestep count and Jacobian age
        # limits are active
        self.create_sim(p=ct.one_atm, Tin=300, reactants="H2:1.1, O2:1, AR:5")
        self.sim.set_time_step(1e-8, [2, 3, 10])
        self.sim.set_max_jac_age(3, 5)
        self.solve_fixed_T(loglevel=3)
        out = capsys.readouterr().out

        assert "Maximum Jacobian age reached (3)" in out
        assert "Maximum Jacobian age reached (5)" in out
        assert "Attempt 2 timesteps" in out
        assert "Attempt 3 timesteps" in out
        assert "Attempt 10 timesteps" in out

    def test_mixture_averaged_case1(self):
        self.run_mix(phi=0.65, T=300, width=0.03, p=1.0, refine=True)

    @pytest.mark.slow_test
    def test_mixture_averaged_case2(self):
        self.run_mix(phi=0.5, T=300, width=2.0, p=1.0, refine=False)

    @pytest.mark.slow_test
    def test_mixture_averaged_case3(self):
        self.run_mix(phi=0.5, T=500, width=0.05, p=1.0, refine=False)

    @pytest.mark.slow_test
    def test_mixture_averaged_case4(self):
        self.run_mix(phi=0.7, T=400, width=2.0, p=5.0, refine=False)

    @pytest.mark.slow_test
    def test_mixture_averaged_case5(self):
        self.run_mix(phi=1.0, T=300, width=2.0, p=5.0, refine=False)

    @pytest.mark.slow_test
    def test_mixture_averaged_case6(self):
        self.run_mix(phi=1.5, T=300, width=0.2, p=1.0, refine=True)

    @pytest.mark.slow_test
    def test_mixture_averaged_case7(self):
        self.run_mix(phi=1.5, T=500, width=2.0, p=0.1, refine=False)

    @pytest.mark.slow_test
    def test_mixture_averaged_case8(self):
        self.run_mix(phi=2.0, T=400, width=2.0, p=5.0, refine=False)

    def test_mixture_averaged_case9(self):
        self.run_mix(phi=0.8, T=180, width=0.05, p=1.0, refine=False)
        assert self.sim.flame.bounds('T')[0] < 190

    def test_adjoint_sensitivities(self):
        self.run_mix(phi=0.5, T=300, width=0.1, p=1.0, refine=True)
        self.sim.flame.set_steady_tolerances(default=(1e-10, 1e-15))
        self.sim.solve(loglevel=0, refine_grid=False)

        # Adjoint sensitivities
        dSdk_adj = self.sim.get_flame_speed_reaction_sensitivities()

        # Forward sensitivities
        dk = 1e-4
        Su0 = self.sim.velocity[0]
        for m in range(self.gas.n_reactions):
            self.gas.set_multiplier(1.0) # reset all multipliers
            self.gas.set_multiplier(1+dk, m) # perturb reaction m
            self.sim.solve(loglevel=0, refine_grid=False)
            Suplus = self.sim.velocity[0]
            self.gas.set_multiplier(1-dk, m) # perturb reaction m
            self.sim.solve(loglevel=0, refine_grid=False)
            Suminus = self.sim.velocity[0]
            fwd = (Suplus-Suminus)/(2*Su0*dk)
            assert fwd == approx(dSdk_adj[m], rel=5e-3, abs=1e-7)

    def test_jacobian_options(self):
        reactants = {'H2': 0.65, 'O2': 0.5, 'AR': 2}
        self.create_sim(p=ct.one_atm, Tin=300, reactants=reactants, width=0.03)
        assert isinstance(self.sim.linear_solver, ct.BandedJacobian)
        self.sim.linear_solver = ct.EigenSparseDirectJacobian()

        self.sim.set_jacobian_perturbation(1e-7, 1e-12, 1e-20)
        self.solve_mix(refine=True)

        # regression value matching test_mixture_averaged_case1
        assert self.sim.velocity[0] == approx(1.693407, rel=1e-4)

    # @utilities.unittest.skip('sometimes slow')
    def test_multicomponent(self):
        reactants = 'H2:1.1, O2:1, AR:5.3'
        p = ct.one_atm
        Tin = 300

        self.create_sim(p, Tin, reactants)
        self.solve_fixed_T()
        self.solve_mix(ratio=5, slope=0.5, curve=0.3)
        Su_mix = self.sim.velocity[0]

        # Multicomponent flame speed should be similar but not identical to
        # the mixture-averaged flame speed.
        self.solve_multi()
        Su_multi = self.sim.velocity[0]
        assert not self.sim.soret_enabled

        assert Su_mix == approx(Su_multi, rel=5e-2)
        assert Su_mix != Su_multi

        # Flame speed with Soret effect should be similar but not identical to
        # the multicomponent flame speed
        self.sim.soret_enabled = True
        self.sim.solve(loglevel=0, refine_grid=True)
        assert self.sim.soret_enabled
        Su_soret = self.sim.velocity[0]

        assert Su_multi == approx(Su_soret, rel=2e-1)
        assert Su_multi != Su_soret

    def test_unity_lewis(self):
        self.create_sim(ct.one_atm, 300, 'H2:1.1, O2:1, AR:5.3')
        self.sim.transport_model = 'unity-Lewis-number'
        self.sim.set_refine_criteria(ratio=3.0, slope=0.08, curve=0.12)
        self.sim.solve(loglevel=0, auto=True)
        dh_unity_lewis = np.ptp(self.sim.enthalpy_mass)

        self.sim.transport_model = 'mixture-averaged'
        self.sim.solve(loglevel=0)
        dh_mix = np.ptp(self.sim.enthalpy_mass)

        # deviation of enthalpy should be much lower for unity Le model (tends
        # towards zero as grid is refined)
        assert dh_unity_lewis < 0.1 * dh_mix

    def test_flux_gradient_mass(self):
        self.create_sim(ct.one_atm, 300, 'H2:1.1, O2:1, AR:5.3')
        self.sim.transport_model = 'mixture-averaged'
        self.sim.set_refine_criteria(ratio=3.0, slope=0.08, curve=0.12)
        assert self.sim.flux_gradient_basis == 'molar'
        self.sim.solve(loglevel=0, auto=True)
        sl_mole = self.sim.velocity[0]
        self.sim.flux_gradient_basis = 'mass'
        assert self.sim.flux_gradient_basis == 'mass'
        self.sim.solve(loglevel=0)
        sl_mass = self.sim.velocity[0]
        # flame speeds should not be exactly equal
        assert sl_mole != sl_mass
        # but they should be close
        assert sl_mole == approx(sl_mass, rel=0.1)

    def test_soret_with_mix(self):
        self.create_sim(ct.one_atm, 300, 'H2:1.1, O2:1, AR:5.3')
        self.sim.transport_model = 'mixture-averaged'
        self.sim.set_refine_criteria(ratio=3.0, slope=0.08, curve=0.12)
        self.sim.solve(loglevel=0, auto=True)
        sl_without_Soret = self.sim.velocity[0]
        self.sim.soret_enabled = True
        self.sim.solve(loglevel=0, auto=True)
        sl_with_Soret = self.sim.velocity[0]
        assert sl_with_Soret < sl_without_Soret
        assert sl_with_Soret > 0
        self.sim.transport_model = 'multicomponent'
        self.sim.solve(loglevel=0, auto=False)
        sl_multi_Soret = self.sim.velocity[0]
        # flame speeds should not be exactly equal
        # but they should be close even on a coarse mesh
        assert abs(sl_multi_Soret-sl_with_Soret) < 0.1

    @pytest.mark.slow_test
    def test_soret_with_auto(self):
        # Test that auto solving with Soret enabled works
        self.create_sim(101325, 300, 'H2:2.0, O2:1.0')
        self.sim.soret_enabled = True
        self.sim.transport_model = 'multicomponent'
        self.sim.solve(loglevel=0, auto=True)

    def test_set_soret_multi_mix(self):
        # Test that the transport model and Soret diffusion
        # can be set in any order without raising errors

        self.create_sim(101325, 300, 'H2:1.0, O2:1.0')
        self.sim.transport_model = 'multicomponent'
        self.sim.soret_enabled = True

        self.sim.transport_model = 'mixture-averaged'
        self.sim.soret_enabled = False

        self.sim.soret_enabled = True
        self.sim.transport_model = 'multicomponent'

    def test_prune(self):
        reactants = 'H2:1.1, O2:1, AR:5'
        p = ct.one_atm
        Tin = 300

        self.create_sim(p, Tin, reactants)
        self.solve_fixed_T()
        self.solve_mix(slope=0.2, curve=0.1, prune=0.0)
        N1 = len(self.sim.grid)

        self.solve_mix(slope=0.5, curve=0.3, prune=0.1)
        N2 = len(self.sim.grid)

        assert N2 < N1

        # TODO: check that the solution is actually correct (that is, that the
        # residual satisfies the error tolerances) on the new grid.

    @pytest.mark.filterwarnings("ignore:.*legacy YAML format.*:UserWarning")
    def test_restore_legacy_yaml(self):
        reactants = 'H2:1.1, O2:1, AR:5'
        p = 5 * ct.one_atm
        Tin = 600
        self.create_sim(p, Tin, reactants)
        meta = self.sim.restore("adiabatic_flame_legacy.yaml", "setup")
        assert meta["generator"] == "Cantera Sim1D"
        assert meta["cantera-version"] == "2.6.0"
        assert self.sim.inlet.T == 300
        assert self.sim.P == approx(ct.one_atm)
        assert len(self.sim.grid) == 9

    def test_fixed_restore_yaml(self):
        # save and restore using YAML format
        filename, group = self.run_fixed_restore("yaml")
        check_component_order(filename, group)

    @pytest.mark.skipif("native" not in ct.hdf_support(),
                        reason="Cantera compiled without HDF support")
    def test_fixed_restore_hdf(self):
        # save and restore using HDF format
        self.run_fixed_restore("h5")

    def run_fixed_restore(self, mode):
        reactants = "H2:1.1, O2:1, AR:5"
        p = 2 * ct.one_atm
        Tin = 400

        self.create_sim(p, Tin, reactants)
        T_rtol = 1.1e-4
        T_atol = 2e-6
        self.sim.flame.set_steady_tolerances(T=(T_rtol, T_atol))

        self.solve_fixed_T()
        filename = self.test_work_path / f"onedim_fixed_T.{mode}"
        filename.unlink(missing_ok=True)

        Y1 = self.sim.Y
        u1 = self.sim.velocity
        P1 = self.sim.P
        T1 = self.sim.T
        self.sim.save(filename, "test")

        # Save a second solution to the same file
        self.sim.radiation_enabled = True
        self.sim.boundary_emissivities = 0.3, 0.8
        self.sim.save(filename, "test2")

        # Create flame object with dummy initial grid
        self.sim = ct.FreeFlame(self.gas)
        self.sim.restore(filename, "test")

        # Sim is initially in "steady-state" mode, so this returns the
        # steady-state tolerances
        rtol, atol = self.sim.flame.tolerances("T")
        assert rtol == approx(T_rtol)
        assert atol == approx(T_atol)
        assert not self.sim.energy_enabled

        P2a = self.sim.P

        assert p == approx(P1)
        assert P1 == approx(P2a)

        Y2 = self.sim.Y
        u2 = self.sim.velocity
        T2 = self.sim.T

        assert Y1 == approx(Y2)
        assert u1 == approx(u2)
        assert T1 == approx(T2)

        self.solve_fixed_T()
        Y3 = self.sim.Y
        u3 = self.sim.velocity

        # TODO: These tolerances seem too loose, but the tests fail on some
        # systems with tighter tolerances.
        assert Y1 == approx(Y3, rel=3e-3)
        assert u1 == approx(u3, rel=1e-3)

        assert not self.sim.radiation_enabled
        assert not self.sim.soret_enabled

        self.sim.restore(filename, "test2")
        assert self.sim.radiation_enabled
        assert self.sim.boundary_emissivities == (0.3, 0.8)

        self.sim.solve(loglevel=0)

        return filename, "test"

    def test_array_properties(self):
        self.create_sim(ct.one_atm, 300, 'H2:1.1, O2:1, AR:5')
        grid_shape = self.sim.grid.shape

        skip = {
            # Skipped because they are constant, irrelevant, or otherwise not desired
            "P", "Te", "atomic_weights", "charges", "electric_potential", "max_temp",
            "min_temp", "molecular_weights", "product_stoich_coeffs",
            "product_stoich_coeffs", "product_stoich_coeffs_reversible",
            "reactant_stoich_coeffs", "reactant_stoich_coeffs", "reference_pressure",
            "state", "u", "v",
            # Skipped because they are 2D (conversion not implemented)
            "binary_diff_coeffs", "creation_rates_ddX", "destruction_rates_ddX",
            "forward_rates_of_progress_ddX", "net_production_rates_ddX",
            "net_rates_of_progress_ddX", "reverse_rates_of_progress_ddX",
            "net_rates_of_progress_ddCi", "forward_rates_of_progress_ddCi",
            "reverse_rates_of_progress_ddCi", "creation_rates_ddCi",
            "destruction_rates_ddCi", "net_production_rates_ddCi"
        }

        for attr in dir(self.gas):
            if attr.startswith("_") or attr in skip:
                continue

            try:
                soln_value = getattr(self.gas, attr)
            except (ct.CanteraError, ct.ThermoModelMethodError, NotImplementedError):
                continue

            if not isinstance(soln_value, (float, np.ndarray)):
                continue

            flame_value = getattr(self.sim, attr)
            assert flame_value.shape == np.asarray(soln_value).shape + grid_shape

    @pytest.mark.slow_test
    def test_save_restore_add_species_yaml(self):
        reactants = "H2:1.1, O2:1, AR:5"
        p = 2 * ct.one_atm
        Tin = 400

        filename = self.test_work_path / "onedim-add-species.yaml"
        filename.unlink(missing_ok=True)
        self.create_sim(p, Tin, reactants, mech="h2o2.yaml")
        gas1 = self.gas
        self.sim.max_grid_points = 234
        self.solve_fixed_T()
        self.solve_mix(ratio=5, slope=0.5, curve=0.3)
        self.sim.transport_model = "multicomponent"
        self.sim.soret_enabled = True
        self.sim.save(filename, "test")
        T1 = self.sim.T
        Y1 = self.sim.Y

        gas2 = ct.Solution("h2o2-plus.yaml", transport_model="multicomponent")
        self.sim = ct.FreeFlame(gas2)
        self.sim.restore(filename, "test")
        T2 = self.sim.T
        Y2 = self.sim.Y

        assert self.sim.soret_enabled
        assert self.sim.max_grid_points == 234
        assert T1 == approx(T2)
        for k1, species in enumerate(gas1.species_names):
            k2 = gas2.species_index(species)
            assert Y1[k1] == approx(Y2[k2])

    @pytest.mark.slow_test
    def test_save_restore_remove_species_yaml(self):
        reactants = "H2:1.1, O2:1, AR:5"
        p = 2 * ct.one_atm
        Tin = 400

        filename = self.test_work_path / "onedim-remove-species.yaml"
        filename.unlink(missing_ok=True)
        self.create_sim(p, Tin, reactants, mech="h2o2-plus.yaml")
        gas1 = self.gas
        self.solve_fixed_T()
        self.solve_mix(ratio=5, slope=0.5, curve=0.3)
        self.sim.save(filename, "test")
        T1 = self.sim.T
        Y1 = self.sim.Y

        gas2 = ct.Solution("h2o2.yaml")
        self.sim = ct.FreeFlame(gas2)
        self.sim.restore(filename, "test")
        T2 = self.sim.T
        Y2 = self.sim.Y

        assert T1 == approx(T2)
        for k2, species in enumerate(gas2.species_names):
            k1 = gas1.species_index(species)
            assert Y1[k1] == approx(Y2[k2])

    def test_write_csv(self):
        filename = self.test_work_path / "onedim-save.csv"
        filename.unlink(missing_ok=True)

        self.create_sim(2e5, 350, 'H2:1.0, O2:2.0', mech="h2o2.yaml")
        self.sim.save(filename, basis="mole")
        data = ct.SolutionArray(self.gas)
        data.read_csv(filename)
        assert data.grid == approx(self.sim.grid)
        assert data.T == approx(self.sim.T)
        k = self.gas.species_index('H2')
        assert data.X[:, k] == approx(self.sim.X[k, :])

    @pytest.mark.usefixtures("allow_deprecated")
    @pytest.mark.filterwarnings("ignore:.*legacy HDF.*:UserWarning")
    @pytest.mark.skipif("native" not in ct.hdf_support(),
                        reason="Cantera compiled without HDF support")
    def test_restore_legacy_hdf(self, test_data_path):
        # Legacy input file was created using the Cantera 2.6 Python test suite:
        # - restore_legacy.h5 -> test_onedim.py::TestFreeFlame::test_write_hdf
        filename = test_data_path / f"freeflame_legacy.h5"

        self.run_mix(phi=1.1, T=350, width=2.0, p=2.0, refine=False)
        desc = 'mixture-averaged simulation'

        f = ct.FreeFlame(self.gas)
        meta = f.restore(filename, "group0")
        assert meta['description'] == desc
        assert meta['cantera_version'] == "2.6.0"

        # check with relaxed tolerances to account for differences between
        # Cantera 2.6 and Cantera 3.1
        self.check_save_restore(f, tol_T=1e-3, tol_X=1e-1)

    @pytest.mark.skipif("native" not in ct.hdf_support(),
                        reason="Cantera compiled without HDF support")
    def test_save_restore_hdf(self):
        # save and restore with native format (HighFive only)
        self.run_save_restore("h5")

    def test_save_restore_yaml(self):
        self.run_save_restore("yaml")

    def run_save_restore(self, mode):
        filename = self.test_work_path / f"freeflame.{mode}"
        filename.unlink(missing_ok=True)

        self.run_mix(phi=1.1, T=350, width=2.0, p=2.0, refine=False)
        desc = 'mixture-averaged simulation'
        self.sim.save(filename, "group0", description=desc)

        f = ct.FreeFlame(self.gas)
        meta = f.restore(filename, "group0")
        assert meta['description'] == desc
        assert meta['cantera-version'] == ct.__version__
        assert meta['git-commit'] == ct.__git_commit__

        self.check_save_restore(f)

    def check_save_restore(self, f, tol_T=None, tol_X=None):
        # approx is used as equality for floats cannot be guaranteed for loaded
        # HDF5 files if they were created on a different OS and/or architecture
        assert list(f.grid) == approx(list(self.sim.grid))
        assert list(f.T) == approx(list(self.sim.T), rel=tol_T)
        k = self.gas.species_index('H2')
        assert list(f.X[k, :]) == approx(list(self.sim.X[k, :]), rel=tol_X)
        assert list(f.inlet.X) == approx(list(self.sim.inlet.X))

        def check_settings(one, two):
            for k, v in one.items():
                assert k in two
                if isinstance(v, dict):
                    for kk, vv in v.items():
                        if isinstance(vv, float):
                            assert two[k][kk] == approx(vv)
                        else:
                            assert two[k][kk] == vv
                elif isinstance(v, float):
                    assert two[k] == approx(v)
                else:
                    assert two[k] == v

        check_settings(self.sim.flame.settings, f.flame.settings)

        f.solve(loglevel=0)

    def test_refine_criteria_boundscheck(self):
        self.create_sim(ct.one_atm, 300.0, 'H2:1.1, O2:1, AR:5')
        good = [3.0, 0.1, 0.2, 0.05]
        bad = [1.2, 1.1, -2, 0.3]

        self.sim.set_refine_criteria(*good)
        for i in range(4):
            with pytest.raises(ct.CanteraError, match='Refiner::setCriteria'):
                vals = list(good)
                vals[i] = bad[i]
                self.sim.set_refine_criteria(*vals)

    def test_refine_criteria(self):
        self.create_sim(ct.one_atm, 300.0, 'H2:1.1, O2:1, AR:5')
        vals = {'ratio': 3.0, 'slope': 0.1, 'curve': 0.2, 'prune': 0.05}
        self.sim.set_refine_criteria(**vals)
        check = self.sim.get_refine_criteria()
        assert vals == check

    def test_replace_grid(self):
        self.create_sim(ct.one_atm, 300.0, 'H2:1.1, O2:1, AR:5')
        self.solve_fixed_T()
        ub = self.sim.velocity[-1]

        # add some points to the grid
        grid = list(self.sim.grid)
        for i in (7,5,4,2):
            grid.insert(i, 0.5 * (grid[i-1] + grid[i]))
        self.sim.flame.grid = grid
        self.sim.set_initial_guess()

        self.solve_fixed_T()
        assert self.sim.velocity[-1] == approx(ub, rel=1e-2)

    def test_exceed_max_grid_points(self):
        self.create_sim(ct.one_atm, 400.0, 'H2:1.1, O2:1, AR:5')
        # Set the maximum grid points to be a low number that should
        # be exceeded by the refinement
        self.sim.max_grid_points = 10
        with pytest.raises(ct.CanteraError, match='max number of grid points'):
            self.sim.solve(loglevel=0, refine_grid=True)

    def test_set_max_grid_points(self):
        self.create_sim(ct.one_atm, 400.0, 'H2:1.1, O2:1, AR:5')
        assert self.sim.max_grid_points == 1000
        self.sim.max_grid_points = 10
        assert self.sim.max_grid_points == 10


class TestDiffusionFlame:
    """
    Note: to re-create the reference file:
    (1) Set PYTHONPATH to build/python.
    (2) Go into test/python directory and run:
        pytest --save-reference=diffusion test_onedim.py::TestDiffusionFlame::test_mixture_averaged
        pytest --save-reference=diffusion test_onedim.py::TestDiffusionFlame::test_auto
        pytest --save-reference=diffusion test_onedim.py::TestDiffusionFlame::test_mixture_averaged_rad
    (3) Compare the reference files created in the current working directory with
        the ones in test/data and replace them if needed.
    """

    def create_sim(self, p, fuel='H2:1.0, AR:1.0', T_fuel=300, mdot_fuel=0.24,
                   oxidizer='O2:0.2, AR:0.8', T_ox=300, mdot_ox=0.72, width=0.02):

        # Solution object used to compute mixture properties
        self.gas = ct.Solution("h2o2.yaml", "ohmech")
        self.gas.TP = T_fuel, p

        # Flame object
        self.sim = ct.CounterflowDiffusionFlame(self.gas, width=width)

        # Set properties of the fuel and oxidizer mixtures
        self.sim.fuel_inlet.mdot = mdot_fuel
        self.sim.fuel_inlet.X = fuel
        self.sim.fuel_inlet.T = T_fuel

        self.sim.oxidizer_inlet.mdot = mdot_ox
        self.sim.oxidizer_inlet.X = oxidizer
        self.sim.oxidizer_inlet.T = T_ox

    def solve_fixed_T(self):
        # Solve with the energy equation disabled
        self.sim.energy_enabled = False
        self.sim.solve(loglevel=0, refine_grid=False)

        assert not self.sim.energy_enabled

    def solve_mix(self, ratio=3.0, slope=0.1, curve=0.12, prune=0.0):
        # Solve with the energy equation enabled

        self.sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)
        self.sim.energy_enabled = True
        self.sim.solve(loglevel=0, refine_grid=True)

        assert self.sim.energy_enabled
        assert self.sim.transport_model == 'mixture-averaged'

        arr = self.sim.to_array()
        radial_lambda = np.unique(arr.Lambda)
        assert len(radial_lambda) == 1
        assert getattr(arr, "lambda")[0] == radial_lambda[0]
        assert arr.radial_pressure_gradient[0] == radial_lambda[0]
        assert arr.radial_pressure_gradient[0] == radial_lambda[0]
        assert getattr(arr, "radial-pressure-gradient")[0] == radial_lambda[0]

    @pytest.mark.slow_test
    def test_mixture_averaged(self, request, test_data_path):
        referenceFile = "DiffusionFlameTest-h2-mix.csv"
        self.create_sim(p=ct.one_atm)
        self.sim.set_initial_guess()

        nPoints = len(self.sim.grid)
        Tfixed = self.sim.T
        self.solve_fixed_T()
        assert nPoints == len(self.sim.grid)
        assert Tfixed == approx(self.sim.T)

        self.solve_mix()
        data = np.empty((self.sim.flame.n_points, self.gas.n_species + 4))
        data[:,0] = self.sim.grid
        data[:,1] = self.sim.velocity
        data[:,2] = self.sim.spread_rate
        data[:,3] = self.sim.T
        data[:,4:] = self.sim.Y.T

        saveReference = request.config.getoption("--save-reference")
        if saveReference == 'diffusion':
            np.savetxt(referenceFile, data, '%11.6e', ', ')
        else:
            bad = compareProfiles(test_data_path / referenceFile, data,
                                  rtol=1e-2, atol=1e-8, xtol=1e-2)
            assert not bad, bad

    def test_auto(self, request, test_data_path):
        referenceFile = "DiffusionFlameTest-h2-auto.csv"
        self.create_sim(p=ct.one_atm, mdot_fuel=2, mdot_ox=3)

        nPoints = []
        timesteps = []

        def steady_func(x):
            nPoints.append(len(self.sim.T))
            return 0

        def time_step_func(dt):
            timesteps.append(dt)
            assert dt > 0
            return 0

        self.sim.set_steady_callback(steady_func)
        self.sim.set_time_step_callback(time_step_func)

        self.sim.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.12, prune=0.0)
        self.sim.solve(loglevel=0, auto=True)

        assert len(nPoints) != 0
        assert len(timesteps) != 0

        data = np.empty((self.sim.flame.n_points, self.gas.n_species + 4))
        data[:,0] = self.sim.grid
        data[:,1] = self.sim.velocity
        data[:,2] = self.sim.spread_rate
        data[:,3] = self.sim.T
        data[:,4:] = self.sim.Y.T

        saveReference = request.config.getoption("--save-reference")
        if saveReference == 'diffusion':
            np.savetxt(referenceFile, data, '%11.6e', ', ')
        else:
            bad = compareProfiles(test_data_path / referenceFile, data,
                                  rtol=1e-2, atol=1e-8, xtol=1e-2)
            assert not bad, bad

    def run_extinction(self, mdot_fuel, mdot_ox, T_ox, width, P):
        self.create_sim(fuel='H2:1.0', oxidizer='O2:1.0', p=ct.one_atm*P,
                        mdot_fuel=mdot_fuel, mdot_ox=mdot_ox, T_ox=T_ox, width=width)
        self.sim.solve(loglevel=0, auto=True)
        assert not self.sim.extinct()

    def test_extinction_case1(self):
        self.run_extinction(mdot_fuel=0.5, mdot_ox=3.0, T_ox=300, width=0.018, P=1.0)

    @pytest.mark.slow_test
    def test_extinction_case2(self):
        self.run_extinction(mdot_fuel=0.5, mdot_ox=1.0, T_ox=300, width=0.01, P=5.0)

    @pytest.mark.slow_test
    def test_extinction_case3(self):
        self.run_extinction(mdot_fuel=1.0, mdot_ox=0.5, T_ox=500, width=0.02, P=5.0)

    @pytest.mark.slow_test
    def test_extinction_case4(self):
        self.run_extinction(mdot_fuel=1.0, mdot_ox=3.0, T_ox=400, width=0.05, P=2.0)

    @pytest.mark.slow_test
    def test_extinction_case5(self):
        self.run_extinction(mdot_fuel=1.0, mdot_ox=3.0, T_ox=300, width=0.1, P=1.0)

    @pytest.mark.slow_test
    def test_extinction_case6(self):
        self.run_extinction(mdot_fuel=0.5, mdot_ox=0.5, T_ox=600, width=0.2, P=0.05)

    @pytest.mark.slow_test
    def test_extinction_case7(self):
        self.run_extinction(mdot_fuel=0.2, mdot_ox=2.0, T_ox=600, width=0.2, P=0.05)

    @pytest.mark.slow_test
    def test_restart(self):
        self.run_extinction(mdot_fuel=0.5, mdot_ox=3.0, T_ox=300, width=0.018, P=1.0)

        arr = self.sim.to_array()

        self.create_sim(mdot_fuel=5.5, mdot_ox=3.3, T_ox=400, width=0.018,
                        p=ct.one_atm*1.1)
        self.sim.set_initial_guess(data=arr)
        self.sim.solve(loglevel=0, auto=True)

        # Check inlet
        mdot = self.sim.density * self.sim.velocity
        assert mdot[0] == approx(self.sim.fuel_inlet.mdot, rel=1e-4)
        assert self.sim.T[0] == approx(self.sim.fuel_inlet.T, rel=1e-4)
        assert mdot[-1] == approx(-self.sim.oxidizer_inlet.mdot, rel=1e-4)

    def test_mixture_averaged_rad(self, request, test_data_path):
        referenceFile = "DiffusionFlameTest-h2-mix-rad.csv"
        self.create_sim(p=ct.one_atm)
        self.sim.set_initial_guess()

        nPoints = len(self.sim.grid)
        Tfixed = self.sim.T
        self.solve_fixed_T()
        assert nPoints == len(self.sim.grid)
        assert Tfixed == approx(self.sim.T)
        assert not self.sim.radiation_enabled
        self.sim.radiation_enabled = True
        assert self.sim.radiation_enabled
        self.sim.flame.boundary_emissivities = 0.25, 0.15

        self.solve_mix()
        data = np.empty((self.sim.flame.n_points, self.gas.n_species + 4))
        data[:,0] = self.sim.grid
        data[:,1] = self.sim.velocity
        data[:,2] = self.sim.spread_rate
        data[:,3] = self.sim.T
        data[:,4:] = self.sim.Y.T

        qdot = self.sim.flame.radiative_heat_loss
        assert len(qdot) == self.sim.flame.n_points

        saveReference = request.config.getoption("--save-reference")
        if saveReference == 'diffusion':
            np.savetxt(referenceFile, data, '%11.6e', ', ')
        else:
            bad = compareProfiles(test_data_path / referenceFile, data,
                                  rtol=1e-2, atol=1e-8, xtol=1e-2)
            assert not bad, bad

        filename = self.test_work_path / "DiffusionFlameTest-h2-mix-rad.csv"
        filename.unlink(missing_ok=True)
        self.sim.save(filename, basis="mole") # check output
        assert filename.is_file()
        csv_data = np.genfromtxt(filename, dtype=float, delimiter=',', names=True)
        assert 'radiativeheatloss' in csv_data.dtype.names

    def test_strain_rate(self):
        # This doesn't test that the values are correct, just that they can be
        # computed without error

        self.create_sim(p=ct.one_atm)
        self.solve_fixed_T()

        a_max = self.sim.strain_rate('max')
        a_mean = self.sim.strain_rate('mean')
        a_pf_fuel = self.sim.strain_rate('potential_flow_fuel')
        a_pf_oxidizer = self.sim.strain_rate('potential_flow_oxidizer')
        a_stoich1 = self.sim.strain_rate('stoichiometric', fuel='H2')
        a_stoich2 = self.sim.strain_rate('stoichiometric', fuel='H2', stoich=0.5)

        assert a_mean <= a_max
        assert a_pf_fuel <= a_max
        assert a_pf_oxidizer <= a_max
        assert a_stoich1 <= a_max
        assert a_stoich1 == a_stoich2

        with pytest.raises(ValueError):
            self.sim.strain_rate('bad_keyword')
        with pytest.raises(KeyError): # missing 'fuel'
            self.sim.strain_rate('stoichiometric')
        with pytest.raises(KeyError): # missing 'stoich'
            self.sim.strain_rate('stoichiometric', fuel='H2', oxidizer='H2O2')

    def test_mixture_fraction(self):
        self.create_sim(p=ct.one_atm)
        self.sim.set_initial_guess()
        Z = self.sim.mixture_fraction('H')
        assert Z[0] == approx(1.0)
        assert Z[-1] == approx(0.0)
        assert all(Z >= 0)
        assert all(Z <= 1.0)
        Z = self.sim.mixture_fraction('Bilger')
        assert Z[0] == approx(1.0)
        assert Z[-1] == approx(0.0)
        assert all(Z >= 0)
        assert all(Z <= 1.0)

    def test_equivalence_ratio(self):
        self.create_sim(p=ct.one_atm)
        self.sim.set_initial_guess()
        phi = self.sim.equivalence_ratio
        assert phi[0] == np.inf
        assert np.isclose(phi[-1], 0.0)

    def test_bad_boundary_conditions(self):
        gas = ct.Solution("h2o2.yaml")
        gas.TP = 300, ct.one_atm
        sim = ct.CounterflowDiffusionFlame(gas, width=0.2)
        sim.fuel_inlet.X = "H2: 1.0"

        sim.oxidizer_inlet.X = "H2: 1.0, O2: 1.0"
        sim.fuel_inlet.X = "O2: 1.0, N2: 3.76"
        with pytest.raises(ct.CanteraError, match="Mass flow for fuel"):
            sim.set_initial_guess()

        sim.fuel_inlet.mdot = 1.0
        sim.oxidizer_inlet.mdot = 1.0
        with pytest.raises(ct.CanteraError, match="too much oxygen"):
            sim.set_initial_guess()

        sim.oxidizer_inlet.X = "O2: 1.0, N2: 3.76"
        sim.fuel_inlet.X = "H2: 1.0"
        sim.set_initial_guess()

    def run_restore_diffusionflame(self, fname):
        gas = ct.Solution("h2o2.yaml")
        sim = ct.CounterflowDiffusionFlame(gas)
        sim.restore(fname)
        radial_lambda = sim.flame.radial_pressure_gradient
        assert np.isclose(np.unique(radial_lambda)[0], -1.355348889404357e+05)
        assert np.allclose(sim.T, sim.flame.temperature)
        assert np.allclose(radial_lambda, sim.flame.L)
        assert np.allclose(radial_lambda, sim.flame.Lambda)
        assert np.allclose(radial_lambda, getattr(sim.flame, "lambda"))
        zeros = sim.T * 0
        with pytest.raises(ct.CanteraError, match="is not used by 'axisymmetric-flow'"):
            sim.flame.set_values("electric-field", zeros)
        with pytest.raises(ct.CanteraError, match="is not used by 'axisymmetric-flow'"):
            sim.flame.set_values("oxidizer-velocity", zeros)

    def run_restore_solutionarray(self, fname):
        gas = ct.Solution("h2o2.yaml")
        arr = ct.SolutionArray(gas)
        arr.restore(fname, "solution/flame")

        radial_lambda = np.unique(arr.radial_pressure_gradient)
        assert len(radial_lambda) == 1
        assert radial_lambda[0] == pytest.approx(-1.355348889404357e+05)
        assert getattr(arr, "lambda")[0] == radial_lambda[0]
        assert radial_lambda[0] == arr.L[0]
        assert radial_lambda[0] == arr.Lambda[0]
        assert getattr(arr, "radial-pressure-gradient")[0] == radial_lambda[0]

    def test_restore_lambda_yaml(self):
        self.run_restore_diffusionflame("flame_lambda.yaml")
        self.run_restore_solutionarray("flame_lambda.yaml")

    @pytest.mark.skipif("native" not in ct.hdf_support(),
                        reason="Cantera compiled without HDF support")
    def test_restore_lambda_hdf(self):
        self.run_restore_diffusionflame("flame_lambda.h5")
        self.run_restore_solutionarray("flame_lambda.h5")

    def test_two_point_control(self):
        """
        Computes a solution using a standard counterflow diffusion flame as a
        starting solution. The tests then do the following:
          1.) Recomputes the solution with two-point control enabled with no
              temperature decrement is used for the two-point control. Solutions
              should be the same.

          2.) Recomputes the solution with two-point control enabled with a fixed
              temperature decrement, then undoes the temperature decrement. Solutions
              should be the same.

          3.) Saves and restores the flame object to ensure that the two-point control
              settings are saved and restored correctly.
        """
        p = ct.one_atm
        fuel = 'H2:1.0, AR:1.0'
        T_fuel = 300
        mdot_fuel = 0.05
        oxidizer = 'O2:0.21, AR:0.78'
        T_ox = 300
        mdot_ox = 0.1
        width = 0.05

        # Solution object used to compute mixture properties
        gas = ct.Solution("h2o2.yaml")
        gas.TP = T_fuel, p

        # Flame object
        sim = ct.CounterflowDiffusionFlame(gas, width=width)

        # Set properties of the fuel and oxidizer mixtures
        sim.fuel_inlet.mdot = mdot_fuel
        sim.fuel_inlet.X = fuel
        sim.fuel_inlet.T = T_fuel

        sim.oxidizer_inlet.mdot = mdot_ox
        sim.oxidizer_inlet.X = oxidizer
        sim.oxidizer_inlet.T = T_ox

        # Test - Check if two-point control gives the same solution as the solver
        #        without two-point control
        sim.flame.set_steady_tolerances(default=[1.0e-5, 1.0e-11])
        sim.solve(loglevel=0)
        temperature_1 = sim.T

        sim.two_point_control_enabled = True
        spacing = 0.95
        control_temperature = np.min(sim.T) + spacing * (np.max(sim.T) - np.min(sim.T))
        sim.set_left_control_point(control_temperature)
        sim.set_right_control_point(control_temperature)
        sim.solve(loglevel=0, refine_grid=False)
        temperature_2 = sim.T

        # Check difference between the two solutions
        assert temperature_1 == approx(temperature_2, rel=1e-4, abs=1e-6)


        # Test - See if the solution is the same when the temperature is decremented and
        # incremented
        temperature_decrement = 5
        sim.right_control_point_temperature -= temperature_decrement
        sim.left_control_point_temperature -= temperature_decrement
        sim.solve(loglevel=0, refine_grid=False)

        sim.right_control_point_temperature += temperature_decrement
        sim.left_control_point_temperature += temperature_decrement
        sim.solve(loglevel=0, refine_grid=False)
        assert sim.flame.oxidizer_velocity == approx(sim.velocity[-1])
        temperature_4 = sim.T

        # Check the difference between the un-perturbed two-point solution and the
        # round-trip solution.
        assert temperature_2 == approx(temperature_4, rel=1e-4, abs=1e-6)

        # Test 4 - Checking property setters/getters
        val_1 = sim.right_control_point_temperature
        increment = 10
        sim.right_control_point_temperature += increment
        assert sim.right_control_point_temperature == val_1 + increment

        val_2 = sim.left_control_point_temperature
        sim.left_control_point_temperature += increment
        assert sim.left_control_point_temperature == val_2 + increment

        # Test - Check if the parameters are saved and restored correctly
        filename = self.test_work_path / "two_point_control.yaml"
        filename.unlink(missing_ok=True)

        original_settings = sim.flame.settings['continuation-method']
        sim.save(filename)

        sim = ct.CounterflowDiffusionFlame(gas, width=width)
        sim.restore(filename)
        filename.unlink()

        assert 'continuation-method' in sim.flame.settings
        restored_settings = sim.flame.settings['continuation-method']

        assert restored_settings['type'] == original_settings['type']

        assert (sim.left_control_point_coordinate
                == approx(original_settings['left-location'], 1e-4))

        assert (sim.right_control_point_coordinate
                == approx(original_settings['right-location'], 1e-4))

        assert (sim.left_control_point_temperature
                == approx(original_settings['left-temperature'], 1e-4))

        assert (sim.right_control_point_temperature
                == approx(original_settings['right-temperature'], 1e-4))

        assert sim.flame.oxidizer_velocity == approx(sim.velocity[-1])

        # Test - Check error conditions
        sim = ct.CounterflowDiffusionFlame(gas, width=width)
        with pytest.raises(ct.CanteraError, match="two-point control is not enabled"):
            sim.right_control_point_temperature = 300

        with pytest.raises(ct.CanteraError, match="two-point control is not enabled"):
            sim.right_control_point_temperature

        with pytest.raises(ct.CanteraError, match="two-point control is not enabled"):
            sim.left_control_point_temperature = 300

        with pytest.raises(ct.CanteraError, match="two-point control is not enabled"):
            sim.left_control_point_temperature

        with pytest.raises(ct.CanteraError,
                           match="No domain with two-point control enabled was found"):
            sim.set_left_control_point(1)

        with pytest.raises(ct.CanteraError,
                           match="No domain with two-point control enabled was found"):
            sim.set_right_control_point(1)

        sim.two_point_control_enabled = True
        with pytest.raises(ct.CanteraError,
                           match="left control point location is not set"):
            sim.left_control_point_coordinate

        with pytest.raises(ct.CanteraError,
                           match="left control point location is not set"):
            sim.left_control_point_temperature

        with pytest.raises(ct.CanteraError,
                           match="right control point location is not set"):
            sim.right_control_point_coordinate

        with pytest.raises(ct.CanteraError,
                           match="right control point location is not set"):
            sim.right_control_point_temperature


class TestCounterflowPremixedFlame:
    """
    Note: to re-create the reference file:
    (1) Set PYTHONPATH to build/python.
    (2) Go into test/python directory and run:
        pytest --save-reference=counterflow_premixed test_onedim.py::TestCounterflowPremixedFlame::test_mixture_averaged
    (3) Compare the reference files created in the current working directory with
        the ones in test/data and replace them if needed.
    """

    def test_mixture_averaged(self, request, test_data_path):
        T_in = 373.0  # inlet temperature
        comp = 'H2:1.6, O2:1, AR:7'  # premixed gas composition

        gas = ct.Solution("h2o2.yaml")
        gas.TPX = T_in, 0.05 * ct.one_atm, comp
        width = 0.2 # m

        sim = ct.CounterflowPremixedFlame(gas=gas, width=width)

        # set the properties at the inlets
        sim.reactants.mdot = 0.12  # kg/m^2/s
        sim.reactants.X = comp
        sim.reactants.T = T_in
        sim.products.mdot = 0.06  # kg/m^2/s

        sim.flame.set_steady_tolerances(default=[1.0e-5, 1.0e-11])
        sim.flame.set_transient_tolerances(default=[1.0e-5, 1.0e-11])
        sim.set_initial_guess()  # assume adiabatic equilibrium products

        sim.energy_enabled = False
        sim.solve(loglevel=0, refine_grid=False)

        sim.set_refine_criteria(ratio=3, slope=0.2, curve=0.4, prune=0.02)
        sim.energy_enabled = True
        assert not sim.radiation_enabled
        sim.solve(loglevel=0, refine_grid=True)

        data = np.empty((sim.flame.n_points, gas.n_species + 4))
        data[:,0] = sim.grid
        data[:,1] = sim.velocity
        data[:,2] = sim.spread_rate
        data[:,3] = sim.T
        data[:,4:] = sim.Y.T

        referenceFile = "CounterflowPremixedFlame-h2-mix.csv"
        saveReference = request.config.getoption("--save-reference")
        if saveReference == 'counterflow_premixed':
            np.savetxt(referenceFile, data, '%11.6e', ', ')
        else:
            bad = compareProfiles(test_data_path / referenceFile, data,
                                  rtol=1e-2, atol=1e-8, xtol=1e-2)
            assert not bad, bad

        filename = self.test_work_path / "CounterflowPremixedFlame-h2-mix.csv"
        filename.unlink(missing_ok=True)

        sim.save(filename) # check output
        assert filename.is_file()
        csv_data = np.genfromtxt(filename, dtype=float, delimiter=',', names=True)
        assert 'qdot' not in csv_data.dtype.names

    @pytest.fixture(scope='function')
    def run_case(self):
        """Helper function to create and run a simulation case."""
        def _run_case(phi, T, width, P):
            gas = ct.Solution("h2o2.yaml")
            gas.TPX = T, P * ct.one_atm, {'H2': phi, 'O2': 0.5, 'AR': 2}
            sim = ct.CounterflowPremixedFlame(gas=gas, width=width)
            sim.reactants.mdot = 10 * gas.density
            sim.products.mdot = 5 * gas.density
            sim.set_refine_criteria(ratio=6, slope=0.7, curve=0.8, prune=0.2)
            sim.solve(loglevel=0, auto=True)
            assert all(sim.T >= T - 1e-3)
            assert all(sim.spread_rate >= -1e-9)
            radial_lambda = sim.flame.radial_pressure_gradient
            assert np.allclose(radial_lambda, radial_lambda[0])
            return sim
        return _run_case

    @pytest.mark.parametrize(
        "phi, T, width, P",
        [
            (0.4, 400, 0.05, 10.0),
            (0.5, 500, 0.03, 2.0),
            (0.7, 300, 0.05, 2.0),
            (1.5, 400, 0.03, 0.02),
            (2.0, 300, 0.2, 0.2),
        ]
    )
    @pytest.mark.slow_test
    def test_solve_case(self, run_case, phi, T, width, P):
        """Parameterized test cases for different simulation setups."""
        run_case(phi, T, width, P)

    @pytest.mark.slow_test
    def test_restart(self, run_case):
        sim = run_case(phi=2.0, T=300, width=0.2, P=0.2)

        arr = sim.to_array()
        sim.reactants.mdot *= 1.1
        sim.products.mdot *= 1.1
        sim.set_initial_guess(data=arr)
        sim.solve(loglevel=0, auto=True)

        # Check inlet / outlet
        mdot = sim.density * sim.velocity
        assert mdot[0] == approx(sim.reactants.mdot, rel=1e-4)
        assert mdot[-1] == approx(-sim.products.mdot, rel=1e-4)

    def test_bad_boundary_conditions(self):
        gas = ct.Solution("h2o2.yaml")
        gas.TPX = 300, 0.05 * ct.one_atm, "H2:1.6, O2:1, AR:7"
        sim = ct.CounterflowPremixedFlame(gas=gas, width=0.2)
        with pytest.raises(ct.CanteraError, match="must be positive"):
            sim.solve(loglevel=0)


class TestCounterflowPremixedFlameNonIdeal:
    """
    Note: to re-create the reference file:
    (1) Set PYTHONPATH to build/python.
    (2) Go into test/python directory and run:
        pytest --save-reference=counterflow_premixed_nonideal test_onedim.py::TestCounterflowPremixedFlameNonIdeal::test_mixture_averaged
    (3) Compare the reference files created in the current working directory with
        the ones in test/data and replace them if needed.
    """

    def test_mixture_averaged(self, request, test_data_path):
        T_in = 373.0  # inlet temperature
        comp = 'H2:1.6, O2:1, AR:7'  # premixed gas composition

        gas = ct.Solution("h2o2.yaml", "ohmech-RK")
        gas.TPX = T_in, 10 * ct.one_atm, comp
        width = 0.005 # m

        sim = ct.CounterflowPremixedFlame(gas=gas, width=width)

        # set the properties at the inlets
        sim.reactants.mdot = 0.12  # kg/m^2/s
        sim.reactants.X = comp
        sim.reactants.T = T_in
        sim.products.mdot = 0.06  # kg/m^2/s

        sim.flame.set_steady_tolerances(default=[1.0e-5, 1.0e-11])
        sim.flame.set_transient_tolerances(default=[1.0e-5, 1.0e-11])
        sim.set_initial_guess()  # assume adiabatic equilibrium products

        sim.energy_enabled = False
        sim.solve(loglevel=0, refine_grid=False)

        sim.set_refine_criteria(ratio=3, slope=0.2, curve=0.4, prune=0.02)
        sim.energy_enabled = True
        assert not sim.radiation_enabled
        sim.solve(loglevel=0, refine_grid=True)

        data = np.empty((sim.flame.n_points, gas.n_species + 4))
        data[:,0] = sim.grid
        data[:,1] = sim.velocity
        data[:,2] = sim.spread_rate
        data[:,3] = sim.T
        data[:,4:] = sim.Y.T

        referenceFile = "CounterflowPremixedFlame-h2-mix-RK.csv"
        saveReference = request.config.getoption("--save-reference")
        if saveReference == 'counterflow_premixed_nonideal':
            np.savetxt(referenceFile, data, '%11.6e', ', ')
        else:
            bad = compareProfiles(test_data_path / referenceFile, data,
                                  rtol=1e-2, atol=1e-8, xtol=1e-2)
            assert not bad, bad

        filename = self.test_work_path / "CounterflowPremixedFlame-h2-mix-RK.csv"
        filename.unlink(missing_ok=True)

        sim.save(filename) # check output
        assert filename.is_file()
        csv_data = np.genfromtxt(filename, dtype=float, delimiter=',', names=True)
        assert 'qdot' not in csv_data.dtype.names

    @pytest.fixture(scope='function')
    def run_case(self):
        """Helper function to create and run a simulation case."""
        def _run_case(phi, T, width, P):
            gas = ct.Solution("h2o2.yaml", "ohmech-RK")
            gas.TPX = T, P * ct.one_atm, {'H2': phi, 'O2': 0.5, 'AR': 3}
            sim = ct.CounterflowPremixedFlame(gas=gas, width=width)
            sim.reactants.mdot = 10 * gas.density
            sim.products.mdot = 5 * gas.density
            sim.set_refine_criteria(ratio=6, slope=0.7, curve=0.8, prune=0.2)
            sim.solve(loglevel=0, auto=True)
            assert all(sim.T >= T - 1e-3)
            assert all(sim.spread_rate >= -1e-9)
            radial_lambda = sim.flame.radial_pressure_gradient
            assert np.allclose(radial_lambda, radial_lambda[0])
            return sim
        return _run_case

    @pytest.mark.parametrize(
        "phi, T, width, P",
        [
            (0.4, 400, 0.05, 10.0),
            (0.5, 500, 0.03, 2.0),
            (0.7, 300, 0.05, 2.0),
            (1.5, 400, 0.03, 0.02),
            (2.0, 300, 0.2, 0.2),
        ]
    )
    @pytest.mark.slow_test
    def test_solve_cases(self, run_case, phi, T, width, P):
        """Parameterized test cases for different simulation setups"""
        run_case(phi, T, width, P)

    @pytest.mark.slow_test
    def test_restart(self, run_case):
        sim = run_case(phi=2.0, T=300, width=0.2, P=0.2)

        arr = sim.to_array()
        sim.reactants.mdot *= 1.1
        sim.products.mdot *= 1.1
        sim.set_initial_guess(data=arr)
        sim.solve(loglevel=0, auto=True)

        # Check inlet / outlet
        mdot = sim.density * sim.velocity
        assert mdot[0] == approx(sim.reactants.mdot, rel=1e-4)
        assert mdot[-1] == approx(-sim.products.mdot, rel=1e-4)


class TestBurnerFlame:

    @pytest.fixture(scope='function')
    def run_case(self):
        """Helper function to create and run a simulation case."""
        def _run_case(phi, T, width, P):
            gas = ct.Solution("h2o2.yaml")
            gas.TPX = T, ct.one_atm*P, {'H2':phi, 'O2':0.5, 'AR':1.5}
            sim = ct.BurnerFlame(gas=gas, width=width)
            sim.burner.mdot = gas.density * 0.15
            sim.solve(loglevel=0, auto=True)
            assert sim.T[1] > T
        return _run_case

    @pytest.mark.parametrize(
        "phi, T, width, P",
        [
            pytest.param(0.5, 500, 2.0, 0.1, id="case1"),
            pytest.param(2.0, 400, 0.05, 1.0, marks=pytest.mark.slow_test, id="case2"),
            pytest.param(1.7, 300, 0.05, 1.0, marks=pytest.mark.slow_test, id="case3"),
            pytest.param(0.5, 300, 1.0, 5.0, marks=pytest.mark.slow_test, id="case4"),
            pytest.param(1.0, 400, 0.2, 0.01, marks=pytest.mark.slow_test, id="case5"),
        ]
    )
    def test_solve_cases(self, run_case, phi, T, width, P):
        """Parameterized test cases for different simulation setups"""
        run_case(phi, T, width, P)

    def test_fixed_temp(self):
        gas = ct.Solution("h2o2.yaml")
        gas.TPX = 400, 2*ct.one_atm, {'H2':0.7, 'O2':0.5, 'AR':1.5}
        sim = ct.BurnerFlame(gas=gas, width=0.05)
        sim.burner.mdot = gas.density * 0.15
        sim.flame.set_fixed_temp_profile([0, 0.1, 0.9, 1],
                                         [400, 1100, 1100, 500])

        sim.energy_enabled = False
        sim.solve(loglevel=0, refine_grid=True)
        assert sim.T[0] == approx(400)
        assert sim.T[-1] == approx(500)
        assert max(sim.T) == approx(1100)

    def test_blowoff(self):
        gas = ct.Solution("h2o2.yaml")
        gas.set_equivalence_ratio(0.4, 'H2', 'O2:1.0, AR:5')
        gas.TP = 300, ct.one_atm
        sim = ct.BurnerFlame(gas=gas, width=0.1)
        sim.burner.mdot = 1.2
        sim.set_refine_criteria(ratio=3, slope=0.3, curve=0.5, prune=0)
        sim.solve(loglevel=0, auto=True)
        # nonreacting solution
        assert sim.T[-1] == approx(sim.T[0], rel=1e-6)
        assert sim.velocity[-1] == approx(sim.velocity[0], rel=1e-6)
        assert sim.Y[:,0] == approx(sim.Y[:,-1], rel=1e-6, abs=1e-6)

    def test_restart(self):
        gas = ct.Solution("h2o2.yaml")
        gas.set_equivalence_ratio(0.8, 'H2', 'O2:1.0, AR:5')
        gas.TP = 300, ct.one_atm
        sim = ct.BurnerFlame(gas=gas, width=0.1)
        sim.burner.mdot = 1.2
        sim.set_refine_criteria(ratio=3, slope=0.3, curve=0.5, prune=0)
        sim.solve(loglevel=0, auto=True)

        arr = sim.to_array()
        sim.burner.mdot = 1.1
        sim.set_initial_guess(data=arr)
        sim.solve(loglevel=0, auto=True)

        # Check continuity
        rhou = sim.burner.mdot
        for rhou_j in sim.density * sim.velocity:
            assert rhou_j == approx(rhou, rel=1e-4)

@pytest.fixture(scope='function')
def setup_stagnation_tests(request):
    request.cls.gas = ct.Solution("h2o2.yaml")

@pytest.mark.usefixtures("setup_stagnation_tests")
class TestStagnationFlame:

    def create_stagnation(self, comp, tsurf, tinlet, mdot, width):
        p = 0.05 * ct.one_atm  # pressure
        self.gas.TPX = tinlet, p, comp

        sim = ct.ImpingingJet(gas=self.gas, width=width)
        sim.inlet.mdot = mdot
        sim.surface.T = tsurf
        return sim

    def run_stagnation(self, xh2, mdot, width):
        # Simplified version of the example 'stagnation_flame.py'
        tburner = 373.0  # burner temperature
        tsurf = 500.0
        comp = {'H2': xh2, 'O2': 1, 'AR': 7}
        sim = self.create_stagnation(comp, tsurf, tburner, mdot, width)

        sim.set_grid_min(1e-4)
        sim.set_refine_criteria(3., 0.1, 0.2, 0.06)
        sim.set_initial_guess(products='equil')  # assume adiabatic equilibrium products

        sim.solve(loglevel=0, auto=True)

        assert sim.T.max() > tburner + tsurf
        radial_lambda = sim.flame.radial_pressure_gradient
        assert np.allclose(radial_lambda, radial_lambda[0])
        self.sim = sim

    def test_stagnation_case1(self):
        self.run_stagnation(xh2=1.8, mdot=0.06, width=0.2)

    @pytest.mark.skipif("native" not in ct.hdf_support(),
                        reason="Cantera compiled without HDF support")
    def test_restore_hdf(self):
        self.run_save_restore("h5")

    def test_restore_yaml(self):
        filename, group = self.run_save_restore("yaml")
        check_component_order(filename, group)

    def run_save_restore(self, mode):
        filename = self.test_work_path / f"stagnation.{mode}"
        filename.unlink(missing_ok=True)

        self.run_stagnation(xh2=1.8, mdot=0.06, width=0.1)
        self.sim.save(filename, "test")

        jet = ct.ImpingingJet(gas=self.gas)
        jet.restore(filename, "test")

        self.check_save_restore(jet)
        return filename, "test"

    def check_save_restore(self, jet):
        # approx is used as equality for floats cannot be guaranteed for loaded
        # HDF5 files if they were created on a different OS and/or architecture
        assert list(jet.grid) == approx(list(self.sim.grid))
        assert list(jet.T) == approx(list(self.sim.T), 1e-3)
        k = self.sim.gas.species_index('H2')
        assert list(jet.X[k, :]) == approx(list(self.sim.X[k, :]), 1e-4)

        settings = self.sim.flame.settings
        for k, v in jet.flame.settings.items():
            assert k in settings
            if k == "fixed_temperature":
                # fixed temperature is NaN
                continue
            if isinstance(v, dict):
                for kk, vv in v.items():
                    if isinstance(vv, float):
                        assert settings[k][kk] == approx(vv)
                    else:
                        assert settings[k][kk] == vv
            if isinstance(k, float):
                assert settings[k] == approx(v)
            else:
                assert settings[k] == v

        jet.solve(loglevel=0)

@pytest.fixture(scope='function')
def setup_impinging_jet_tests(request):
    request.cls.surf_phase = ct.Interface("ptcombust-simple.yaml", "Pt_surf")
    request.cls.gas = request.cls.surf_phase.adjacent["gas"]

@pytest.mark.usefixtures("setup_impinging_jet_tests")
class TestImpingingJet:

    def create_reacting_surface(self, comp, tsurf, tinlet, width):
        self.gas.TPX = tinlet, ct.one_atm, comp
        self.surf_phase.TP = tsurf, ct.one_atm

        # integrate the coverage equations holding the gas composition fixed
        # to generate a good starting estimate for the coverages.
        self.surf_phase.advance_coverages(1.)

        return ct.ImpingingJet(gas=self.gas, width=width, surface=self.surf_phase)

    def run_reacting_surface(self, xch4, tsurf, mdot, width):
        # Simplified version of the example 'catalytic_combustion.py'
        tinlet = 300.0  # inlet temperature
        comp = {'CH4': xch4, 'O2': 0.21, 'N2': 0.79}
        sim = self.create_reacting_surface(comp, tsurf, tinlet, width)
        sim.set_refine_criteria(10.0, 0.3, 0.4, 0.0)

        sim.inlet.mdot = mdot
        sim.inlet.T = tinlet
        sim.inlet.X = comp
        sim.surface.T = tsurf

        sim.solve(loglevel=0, auto=True)

        assert all(np.diff(sim.T) > 0)
        assert all(np.diff(sim.Y[sim.gas.species_index('CH4')]) < 0)
        assert all(np.diff(sim.Y[sim.gas.species_index('CO2')]) > 0)
        self.sim = sim

    def test_reacting_surface_case1(self):
        self.run_reacting_surface(xch4=0.095, tsurf=900.0, mdot=0.06, width=0.1)

    @pytest.mark.slow_test
    def test_reacting_surface_case2(self):
        self.run_reacting_surface(xch4=0.07, tsurf=1200.0, mdot=0.2, width=0.05)

    @pytest.mark.slow_test
    def test_reacting_surface_case3(self):
        self.run_reacting_surface(xch4=0.2, tsurf=800.0, mdot=0.1, width=0.08)

    @pytest.mark.usefixtures("allow_deprecated")
    @pytest.mark.skipif("native" not in ct.hdf_support(),
                        reason="Cantera compiled without HDF support")
    @pytest.mark.filterwarnings("ignore:.*legacy HDF.*:UserWarning")
    def test_restore_legacy_hdf(self, test_data_path):
        # Legacy input file was created using the Cantera 2.6 Python test suite:
        # - restore_legacy.h5 -> test_onedim.py::TestImpingingJet::test_write_hdf
        filename = test_data_path / f"impingingjet_legacy.h5"

        self.run_reacting_surface(xch4=0.095, tsurf=900.0, mdot=0.06, width=0.1)
        jet = ct.ImpingingJet(gas=self.gas, surface=self.surf_phase)
        jet.restore(filename, "group0")

        # check with relaxed tolerances to account for differences between
        # Cantera 2.6 and Cantera 3.1
        self.check_save_restore(jet, tol_T=1e-3, tol_X=1e-1)

    @pytest.mark.skipif("native" not in ct.hdf_support(),
                        reason="Cantera compiled without HDF support")
    def test_restore_hdf(self):
        self.run_save_restore("h5")

    def test_restore_yaml(self):
        filename, group = self.run_save_restore("yaml")
        check_component_order(filename, group)

    def run_save_restore(self, mode):
        filename = self.test_work_path / f"impingingjet.{mode}"
        filename.unlink(missing_ok=True)

        self.run_reacting_surface(xch4=0.095, tsurf=900.0, mdot=0.06, width=0.1)
        self.sim.save(filename, "test")

        comp = {'CH4': 0.095, 'O2':0.21, 'N2':0.79}
        jet = self.create_reacting_surface(comp, 700.0, 500., width=0.2)
        jet.restore(filename, "test")

        self.check_save_restore(jet)
        return filename, "test"

    def check_save_restore(self, jet, tol_T=None, tol_X=None):
        # approx is used as equality for floats cannot be guaranteed for loaded
        # HDF5 files if they were created on a different OS and/or architecture
        assert list(jet.grid) == approx(list(self.sim.grid))
        assert list(jet.T) == approx(list(self.sim.T), tol_T)
        k = self.sim.gas.species_index('H2')
        assert list(jet.X[k, :]) == approx(list(self.sim.X[k, :]), tol_X)

        settings = self.sim.flame.settings
        for k, v in jet.flame.settings.items():
            assert k in settings
            if k == "fixed_temperature":
                # fixed temperature is NaN
                continue
            if isinstance(v, dict):
                for kk, vv in v.items():
                    if isinstance(vv, float):
                        assert settings[k][kk] == approx(vv)
                    else:
                        assert settings[k][kk] == vv
            if isinstance(k, float):
                assert settings[k] == approx(v)
            else:
                assert settings[k] == v

        assert list(jet.surface.surface.X) == approx(list(self.sim.surface.surface.X))
        for k in self.sim.surface.phase.species_names:
            assert self.sim.surface.value(k) == \
                approx(jet.surface.value(k), tol_X)

        jet.solve(loglevel=0)


class TestTwinFlame:
    def solve(self, phi, T, width, P):
        gas = ct.Solution("h2o2.yaml")
        gas.TP = T, ct.one_atm
        gas.set_equivalence_ratio(phi, 'H2', 'O2:1.0, AR:4.0')
        sim = ct.CounterflowTwinPremixedFlame(gas=gas, width=width)
        sim.set_refine_criteria(ratio=5, slope=0.8, curve=1.0, prune=0.2)
        axial_velocity = 2.0
        sim.reactants.mdot = gas.density * axial_velocity
        sim.solve(loglevel=0, auto=True)
        assert sim.T[-1] > T + 100
        radial_lambda = sim.flame.radial_pressure_gradient
        assert np.allclose(radial_lambda, radial_lambda[0])
        return sim

    def test_restart(self):
        sim = self.solve(phi=0.4, T=300, width=0.05, P=0.1)

        arr = sim.to_array()
        axial_velocity = 2.2
        sim.reactants.mdot *= 1.1
        sim.reactants.T = sim.reactants.T + 100
        sim.set_initial_guess(data=arr)
        sim.solve(loglevel=0, auto=True)

        # Check inlet
        mdot = sim.density * sim.velocity
        assert mdot[0] == approx(sim.reactants.mdot, rel=1e-4)
        assert sim.T[0] == approx(sim.reactants.T, rel=1e-4)

    def test_save_restore_yaml(self):
        # save and restore using YAML format
        filename, group = self.run_save_restore("yaml")
        check_component_order(filename, group)

    @pytest.mark.skipif("native" not in ct.hdf_support(),
                        reason="Cantera compiled without HDF support")
    def test_save_restore_hdf(self):
        # save and restore using HDF format
        self.run_save_restore("h5")

    def run_save_restore(self, mode):
        filename = self.test_work_path / f"twinflame.{mode}"
        filename.unlink(missing_ok=True)

        sim = self.solve(phi=0.4, T=300, width=0.05, P=0.1)
        sim.save(filename, compression=7)

        gas = ct.Solution("h2o2.yaml")
        sim2 = ct.CounterflowTwinPremixedFlame(gas=gas)
        sim2.restore(filename)

        assert sim.grid == approx(sim2.grid)
        assert sim.Y == approx(sim2.Y)

        sim2.solve(loglevel=0)
        return filename, "solution"

    def test_bad_boundary_conditions(self):
        gas = ct.Solution("h2o2.yaml")
        gas.TPX = 300, 0.05 * ct.one_atm, "H2:1.6, O2:1, AR:7"
        sim = ct.CounterflowTwinPremixedFlame(gas=gas, width=0.2)
        with pytest.raises(ct.CanteraError, match="must be positive"):
            sim.solve(loglevel=0)

class TestIonFreeFlame:
    def test_ion_profile(self):
        reactants = 'CH4:0.216, O2:2'
        p = ct.one_atm
        Tin = 300
        width = 0.03

        # Solution object used to compute mixture properties
        self.gas = ct.Solution('ch4_ion.yaml')
        self.gas.TPX = Tin, p, reactants
        self.sim = ct.FreeFlame(self.gas, width=width)
        assert self.sim.transport_model == 'ionized-gas'
        self.sim.set_refine_criteria(ratio=4, slope=0.8, curve=1.0)
        # Ionized species may require tighter absolute tolerances
        self.sim.flame.set_steady_tolerances(Y=(1e-4, 1e-12))

        # stage one
        assert self.sim.electric_field_enabled is False
        self.sim.solve(loglevel=0, auto=True)
        assert max(np.abs(self.sim.flame.electric_field)) == approx(0.0, abs=1e-3)

        #stage two
        self.sim.electric_field_enabled = True
        self.sim.solve(loglevel=0)

        # Regression test
        assert max(self.sim.flame.electric_field) == approx(149.6317905667685, rel=1e-3)
        assert np.allclose(self.sim.E, self.sim.flame.electric_field)
        zeros = self.sim.T
        with pytest.raises(ct.CanteraError, match="is not used by 'free-ion-flow'"):
            self.sim.flame.set_values("oxidizer-velocity", zeros)


class TestIonBurnerFlame:
    def test_ion_profile(self):
        reactants = 'CH4:1.0, O2:2.0, N2:7.52'
        p = ct.one_atm
        Tburner = 400
        width = 0.01

        # Solution object used to compute mixture properties
        self.gas = ct.Solution('ch4_ion.yaml')
        self.gas.TPX = Tburner, p, reactants
        self.sim = ct.BurnerFlame(self.gas, width=width)
        assert self.sim.transport_model == 'ionized-gas'
        self.sim.set_refine_criteria(ratio=4, slope=0.1, curve=0.1)
        self.sim.burner.mdot = self.gas.density * 0.15

        self.sim.electric_field_enabled = True
        self.sim.solve(loglevel=0, auto=True)

        # Regression test
        assert max(self.sim.flame.electric_field) == approx(591.76, rel=1e-2)
        assert max(self.sim.X[self.gas.species_index('E')]) == approx(8.024e-9, rel=1e-2)
