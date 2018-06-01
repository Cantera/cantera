import cantera as ct
from . import utilities
import numpy as np
import os
from os.path import join as pjoin


class TestOnedim(utilities.CanteraTest):
    def test_instantiate(self):
        gas = ct.Solution('h2o2.xml')

        flame = ct.FreeFlow(gas)

    def test_badInstantiate(self):
        solid = ct.Solution('diamond.xml', 'diamond')
        with self.assertRaises(TypeError):
            flame = ct.FreeFlow(solid)

    def test_instantiateSurface(self):
        gas = ct.Solution('diamond.xml', 'gas')
        solid = ct.Solution('diamond.xml', 'diamond')
        interface = ct.Solution('diamond.xml', 'diamond_100', (gas, solid))

        surface = ct.ReactingSurface1D(phase=gas)
        surface.set_kinetics(interface)

    def test_boundaryProperties(self):
        gas1 = ct.Solution('h2o2.xml')
        gas2 = ct.Solution('h2o2.xml')
        inlet = ct.Inlet1D(name='something', phase=gas1)
        flame = ct.FreeFlow(gas1)
        sim = ct.Sim1D((inlet, flame))

        self.assertEqual(inlet.name, 'something')

        gas2.TPX = 400, 101325, 'H2:0.3, O2:0.5, AR:0.2'
        Xref = gas2.X
        Yref = gas2.Y
        inlet.Y = Yref

        self.assertArrayNear(inlet.Y, Yref)
        self.assertArrayNear(inlet.X, Xref)

        gas2.TPX = 400, 101325, 'H2:0.5, O2:0.2, AR:0.3'
        Xref = gas2.X
        Yref = gas2.Y
        inlet.X = Xref
        self.assertArrayNear(inlet.X, Xref)
        self.assertArrayNear(inlet.Y, Yref)

        inlet.X = {'H2':0.3, 'O2':0.5, 'AR':0.2}
        self.assertNear(inlet.X[gas2.species_index('H2')], 0.3)

    def test_grid_check(self):
        gas = ct.Solution('h2o2.xml')
        flame = ct.FreeFlow(gas)

        with self.assertRaises(ct.CanteraError):
            flame.grid = [0, 0.1, 0.1, 0.2]

        with self.assertRaises(ct.CanteraError):
            flame.grid = [0, 0.1, 0.2, 0.05]

    def test_unpicklable(self):
        import pickle
        gas = ct.Solution('h2o2.xml')
        flame = ct.FreeFlow(gas)
        with self.assertRaises(NotImplementedError):
            pickle.dumps(flame)

    def test_uncopyable(self):
        import copy
        gas = ct.Solution('h2o2.xml')
        flame = ct.FreeFlow(gas)
        with self.assertRaises(NotImplementedError):
            copy.copy(flame)

    def test_invalid_property(self):
        gas1 = ct.Solution('h2o2.xml')
        inlet = ct.Inlet1D(name='something', phase=gas1)
        flame = ct.FreeFlow(gas1)
        sim = ct.Sim1D((inlet, flame))

        for x in (inlet, flame, sim):
            with self.assertRaises(AttributeError):
                x.foobar = 300
            with self.assertRaises(AttributeError):
                x.foobar

    def test_tolerances(self):
        gas = ct.Solution('h2o2.xml')
        left = ct.Inlet1D(gas)
        flame = ct.FreeFlow(gas)
        right = ct.Inlet1D(gas)
        # Some things don't work until the domains have been added to a Sim1D
        sim = ct.Sim1D((left, flame, right))

        with self.assertRaises(ct.CanteraError):
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

        self.assertEqual(atol_ss, set((5e-5, 3e-6, 7e-9)))
        self.assertEqual(atol_ts, set((6e-5, 4e-6, 2e-9)))
        self.assertEqual(rtol_ss, set((5e-3, 3e-4, 7e-7)))
        self.assertEqual(rtol_ts, set((6e-3, 4e-4, 2e-7)))


class TestFreeFlame(utilities.CanteraTest):
    tol_ss = [1.0e-5, 1.0e-14]  # [rtol atol] for steady-state problem
    tol_ts = [1.0e-4, 1.0e-11]  # [rtol atol] for time stepping

    def create_sim(self, p, Tin, reactants, width=0.05, mech='h2o2.xml'):
        # IdealGasMix object used to compute mixture properties
        self.gas = ct.Solution(mech)
        self.gas.TPX = Tin, p, reactants

        # Flame object
        self.sim = ct.FreeFlame(self.gas, width=width)
        self.sim.flame.set_steady_tolerances(default=self.tol_ss)
        self.sim.flame.set_transient_tolerances(default=self.tol_ts)

        # Set properties of the upstream fuel-air mixture
        self.sim.inlet.T = Tin
        self.sim.inlet.X = reactants

    def solve_fixed_T(self):
        # Solve with the energy equation disabled
        self.sim.energy_enabled = False
        self.sim.solve(loglevel=0, refine_grid=False)

        self.assertFalse(self.sim.energy_enabled)

    def solve_mix(self, ratio=3.0, slope=0.3, curve=0.2, prune=0.0, refine=True):
        # Solve with the energy equation enabled
        self.sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)
        self.sim.energy_enabled = True
        self.sim.solve(loglevel=0, refine_grid=refine)

        self.assertTrue(self.sim.energy_enabled)
        self.assertEqual(self.sim.transport_model, 'Mix')

    def solve_multi(self):
        self.sim.transport_model = 'Multi'
        self.sim.solve(loglevel=0, refine_grid=True)

        self.assertEqual(self.sim.transport_model, 'Multi')

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
        self.assertNear(Tad, self.sim.T[-1], 2e-2)

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

        self.assertLess(abs(T2-Tad), abs(T1-Tad))
        self.assertLess(abs(T3-Tad), abs(T2-Tad))

        for k in range(self.gas.n_species):
            self.assertLess(abs(X2[k]-Xad[k]), abs(X1[k]-Xad[k]))
            self.assertLess(abs(X3[k]-Xad[k]), abs(X2[k]-Xad[k]))

    def run_mix(self, phi, T, width, p, refine):
        reactants = {'H2': phi, 'O2': 0.5, 'AR': 2}
        self.create_sim(p * ct.one_atm, T, reactants, width)
        self.solve_mix(refine=refine)

        rhou = self.sim.inlet.mdot

        # Check continuity
        for rhou_j in self.sim.density * self.sim.u:
            self.assertNear(rhou_j, rhou, 1e-4)

    def test_mixture_averaged_case1(self):
        self.run_mix(phi=0.65, T=300, width=0.03, p=1.0, refine=True)

    def test_mixture_averaged_case2(self):
        self.run_mix(phi=0.5, T=300, width=2.0, p=1.0, refine=False)

    def test_mixture_averaged_case3(self):
        self.run_mix(phi=0.5, T=500, width=0.05, p=1.0, refine=False)

    def test_mixture_averaged_case4(self):
        self.run_mix(phi=0.7, T=400, width=2.0, p=5.0, refine=False)

    def test_mixture_averaged_case5(self):
        self.run_mix(phi=1.0, T=300, width=2.0, p=5.0, refine=False)

    def test_mixture_averaged_case6(self):
        self.run_mix(phi=1.5, T=300, width=0.2, p=1.0, refine=True)

    def test_mixture_averaged_case7(self):
        self.run_mix(phi=1.5, T=500, width=2.0, p=0.1, refine=False)

    def test_mixture_averaged_case8(self):
        self.run_mix(phi=2.0, T=400, width=2.0, p=5.0, refine=False)

    def test_adjoint_sensitivities(self):
        self.run_mix(phi=0.5, T=300, width=0.1, p=1.0, refine=True)
        self.sim.flame.set_steady_tolerances(default=(1e-10, 1e-15))
        self.sim.solve(loglevel=0, refine_grid=False)

        # Adjoint sensitivities
        dSdk_adj = self.sim.get_flame_speed_reaction_sensitivities()

        # Forward sensitivities
        dk = 1e-4
        Su0 = self.sim.u[0]
        for m in range(self.gas.n_reactions):
            self.gas.set_multiplier(1.0) # reset all multipliers
            self.gas.set_multiplier(1+dk, m) # perturb reaction m
            self.sim.solve(loglevel=0, refine_grid=False)
            Suplus = self.sim.u[0]
            self.gas.set_multiplier(1-dk, m) # perturb reaction m
            self.sim.solve(loglevel=0, refine_grid=False)
            Suminus = self.sim.u[0]
            fwd = (Suplus-Suminus)/(2*Su0*dk)
            self.assertNear(fwd, dSdk_adj[m], 5e-3)

    # @utilities.unittest.skip('sometimes slow')
    def test_multicomponent(self):
        reactants = 'H2:1.1, O2:1, AR:5.3'
        p = ct.one_atm
        Tin = 300

        self.create_sim(p, Tin, reactants)
        self.solve_fixed_T()
        self.solve_mix(ratio=5, slope=0.5, curve=0.3)
        Su_mix = self.sim.u[0]

        # Multicomponent flame speed should be similar but not identical to
        # the mixture-averaged flame speed.
        self.solve_multi()
        Su_multi = self.sim.u[0]
        self.assertFalse(self.sim.soret_enabled)

        self.assertNear(Su_mix, Su_multi, 5e-2)
        self.assertNotEqual(Su_mix, Su_multi)

        # Flame speed with Soret effect should be similar but not identical to
        # the multicomponent flame speed
        self.sim.soret_enabled = True
        self.sim.solve(loglevel=0, refine_grid=True)
        self.assertTrue(self.sim.soret_enabled)
        Su_soret = self.sim.u[0]

        self.assertNear(Su_multi, Su_soret, 2e-1)
        self.assertNotEqual(Su_multi, Su_soret)

    def test_unity_lewis(self):
        self.create_sim(ct.one_atm, 300, 'H2:1.1, O2:1, AR:5.3')
        self.sim.transport_model = 'UnityLewis'
        self.sim.set_refine_criteria(ratio=3.0, slope=0.08, curve=0.12)
        self.sim.solve(loglevel=0, auto=True)
        dh_unity_lewis = self.sim.enthalpy_mass.ptp()

        self.sim.transport_model = 'Mix'
        self.sim.solve(loglevel=0)
        dh_mix = self.sim.enthalpy_mass.ptp()

        # deviation of enthalpy should be much lower for unity Le model (tends
        # towards zero as grid is refined)
        self.assertLess(dh_unity_lewis, 0.1 * dh_mix)

    def test_soret_with_mix(self):
        # Test that enabling Soret diffusion without
        # multicomponent transport results in an error

        self.create_sim(101325, 300, 'H2:1.0, O2:1.0')
        self.assertFalse(self.sim.soret_enabled)
        self.assertFalse(self.sim.transport_model == 'Multi')

        with self.assertRaises(ct.CanteraError):
            self.sim.soret_enabled = True
            self.sim.solve(loglevel=0, auto=False)

    def test_soret_with_auto(self):
        # Test that auto solving with Soret enabled works
        self.create_sim(101325, 300, 'H2:2.0, O2:1.0')
        self.sim.soret_enabled = True
        self.sim.transport_model = 'Multi'
        self.sim.solve(loglevel=0, auto=True)

    def test_set_soret_multi_mix(self):
        # Test that the transport model and Soret diffusion
        # can be set in any order without raising errors

        self.create_sim(101325, 300, 'H2:1.0, O2:1.0')
        self.sim.transport_model = 'Multi'
        self.sim.soret_enabled = True

        self.sim.transport_model = 'Mix'
        self.sim.soret_enabled = False

        self.sim.soret_enabled = True
        self.sim.transport_model = 'Multi'

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

        self.assertLess(N2, N1)

        # TODO: check that the solution is actually correct (i.e. that the
        # residual satisfies the error tolerances) on the new grid.

    def test_save_restore(self):
        reactants = 'H2:1.1, O2:1, AR:5'
        p = 2 * ct.one_atm
        Tin = 400

        self.create_sim(p, Tin, reactants)
        self.solve_fixed_T()
        filename = pjoin(self.test_work_dir, 'onedim-fixed-T{0}.xml'.format(utilities.python_version))
        if os.path.exists(filename):
            os.remove(filename)

        Y1 = self.sim.Y
        u1 = self.sim.u
        V1 = self.sim.V
        P1 = self.sim.P

        self.sim.save(filename, 'test', loglevel=0)

        # Save a second solution to the same file
        self.sim.save(filename, 'test2', loglevel=0)

        # Create flame object with dummy initial grid
        self.sim = ct.FreeFlame(self.gas)
        self.sim.restore(filename, 'test', loglevel=0)

        # Sim is initially in "steady-state" mode, so this returns the
        # steady-state tolerances
        rtol, atol = self.sim.flame.tolerances('T')
        self.assertNear(rtol, self.tol_ss[0])
        self.assertNear(atol, self.tol_ss[1])
        self.assertFalse(self.sim.energy_enabled)

        P2a = self.sim.P

        self.assertNear(p, P1)
        self.assertNear(P1, P2a)

        Y2 = self.sim.Y
        u2 = self.sim.u
        V2 = self.sim.V

        self.assertArrayNear(Y1, Y2)
        self.assertArrayNear(u1, u2)
        self.assertArrayNear(V1, V2)

        self.solve_fixed_T()
        Y3 = self.sim.Y
        u3 = self.sim.u
        V3 = self.sim.V

        # TODO: These tolereances seem too loose, but the tests fail on some
        # systems with tighter tolerances.
        self.assertArrayNear(Y1, Y3, 3e-3)
        self.assertArrayNear(u1, u3, 1e-3)
        self.assertArrayNear(V1, V3, 1e-3)

    def test_array_properties(self):
        self.create_sim(ct.one_atm, 300, 'H2:1.1, O2:1, AR:5')

        for attr in ct.FlameBase.__dict__:
            if isinstance(ct.FlameBase.__dict__[attr], property):
                getattr(self.sim, attr)

    def test_save_restore_add_species(self):
        reactants = 'H2:1.1, O2:1, AR:5'
        p = 2 * ct.one_atm
        Tin = 400

        filename = pjoin(self.test_work_dir, 'onedim-add-species{0}.xml'.format(utilities.python_version))
        if os.path.exists(filename):
            os.remove(filename)

        self.create_sim(p, Tin, reactants, mech='h2o2.xml')
        gas1 = self.gas
        self.solve_fixed_T()
        self.solve_mix(ratio=5, slope=0.5, curve=0.3)
        self.sim.save(filename, 'test', loglevel=0)
        T1 = self.sim.T
        Y1 = self.sim.Y

        gas2 = ct.Solution('h2o2-plus.xml')
        self.sim = ct.FreeFlame(gas2)
        self.sim.restore(filename, 'test', loglevel=0)
        T2 = self.sim.T
        Y2 = self.sim.Y

        self.assertArrayNear(T1, T2)
        for k1, species in enumerate(gas1.species_names):
            k2 = gas2.species_index(species)
            self.assertArrayNear(Y1[k1], Y2[k2])

    def test_save_restore_remove_species(self):
        reactants = 'H2:1.1, O2:1, AR:5'
        p = 2 * ct.one_atm
        Tin = 400

        filename = pjoin(self.test_work_dir, 'onedim-add-species{0}.xml'.format(utilities.python_version))
        if os.path.exists(filename):
            os.remove(filename)

        self.create_sim(p, Tin, reactants, mech='h2o2-plus.xml')
        gas1 = self.gas
        self.solve_fixed_T()
        self.solve_mix(ratio=5, slope=0.5, curve=0.3)
        self.sim.save(filename, 'test', loglevel=0)
        T1 = self.sim.T
        Y1 = self.sim.Y

        gas2 = ct.Solution('h2o2.xml')
        self.sim = ct.FreeFlame(gas2)
        self.sim.restore(filename, 'test', loglevel=0)
        T2 = self.sim.T
        Y2 = self.sim.Y

        self.assertArrayNear(T1, T2)
        for k2, species in enumerate(gas2.species_names):
            k1 = gas1.species_index(species)
            self.assertArrayNear(Y1[k1], Y2[k2])

    def test_write_csv(self):
        filename = pjoin(self.test_work_dir, 'onedim-write_csv{0}.csv'.format(utilities.python_version))
        if os.path.exists(filename):
            os.remove(filename)

        self.create_sim(2e5, 350, 'H2:1.0, O2:2.0', mech='h2o2.xml')
        self.sim.write_csv(filename)
        data = np.genfromtxt(filename, delimiter=',', skip_header=1)
        self.assertArrayNear(data[:,0], self.sim.grid)
        self.assertArrayNear(data[:,3], self.sim.T)
        k = self.gas.species_index('H2')
        self.assertArrayNear(data[:,5+k], self.sim.X[k,:])

    def test_refine_criteria_boundscheck(self):
        self.create_sim(ct.one_atm, 300.0, 'H2:1.1, O2:1, AR:5')
        good = [3.0, 0.1, 0.2, 0.05]
        bad = [1.2, 1.1, -2, 0.3]

        self.sim.set_refine_criteria(*good)
        for i in range(4):
            with self.assertRaises(ct.CanteraError):
                vals = list(good)
                vals[i] = bad[i]
                self.sim.set_refine_criteria(*vals)

    def test_refine_criteria(self):
        self.create_sim(ct.one_atm, 300.0, 'H2:1.1, O2:1, AR:5')
        vals = {'ratio': 3.0, 'slope': 0.1, 'curve': 0.2, 'prune': 0.05}
        self.sim.set_refine_criteria(**vals)
        check = self.sim.get_refine_criteria()
        self.assertEqual(vals, check)

    def test_replace_grid(self):
        self.create_sim(ct.one_atm, 300.0, 'H2:1.1, O2:1, AR:5')
        self.solve_fixed_T()
        ub = self.sim.u[-1]

        # add some points to the grid
        grid = list(self.sim.grid)
        for i in (7,5,4,2):
            grid.insert(i, 0.5 * (grid[i-1] + grid[i]))
        self.sim.flame.grid = grid
        self.sim.set_initial_guess()

        self.solve_fixed_T()
        self.assertNear(self.sim.u[-1], ub, 1e-2)

    def test_exceed_max_grid_points(self):
        self.create_sim(ct.one_atm, 400.0, 'H2:1.1, O2:1, AR:5')
        # Set the maximum grid points to be a low number that should
        # be exceeded by the refinement
        self.sim.max_grid_points = 10
        with self.assertRaises(ct.CanteraError):
            self.sim.solve(loglevel=0, refine_grid=True)

    def test_set_max_grid_points(self):
        self.create_sim(ct.one_atm, 400.0, 'H2:1.1, O2:1, AR:5')
        self.assertEqual(self.sim.max_grid_points, 1000)
        self.sim.max_grid_points = 10
        self.assertEqual(self.sim.max_grid_points, 10)


class TestDiffusionFlame(utilities.CanteraTest):
    # Note: to re-create the reference file:
    # (1) set PYTHONPATH to build/python2 or build/python3.
    # (2) Start Python in the test/work directory and run:
    #     >>> import cantera.test
    #     >>> t = cantera.test.test_onedim.TestDiffusionFlame("test_mixture_averaged")
    #     >>> t.test_mixture_averaged(True)

    def create_sim(self, p, fuel='H2:1.0, AR:1.0', T_fuel=300, mdot_fuel=0.24,
                   oxidizer='O2:0.2, AR:0.8', T_ox=300, mdot_ox=0.72, width=0.02):

        # IdealGasMix object used to compute mixture properties
        self.gas = ct.Solution('h2o2.xml', 'ohmech')
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

        self.sim.set_initial_guess()

    def solve_fixed_T(self):
        # Solve with the energy equation disabled
        self.sim.energy_enabled = False
        self.sim.solve(loglevel=0, refine_grid=False)

        self.assertFalse(self.sim.energy_enabled)

    def solve_mix(self, ratio=3.0, slope=0.1, curve=0.12, prune=0.0):
        # Solve with the energy equation enabled

        self.sim.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)
        self.sim.energy_enabled = True
        self.sim.solve(loglevel=0, refine_grid=True)

        self.assertTrue(self.sim.energy_enabled)
        self.assertEqual(self.sim.transport_model, 'Mix')

    def test_mixture_averaged(self, saveReference=False):
        referenceFile = pjoin(self.test_data_dir, 'DiffusionFlameTest-h2-mix.csv')
        self.create_sim(p=ct.one_atm)

        nPoints = len(self.sim.grid)
        Tfixed = self.sim.T
        self.solve_fixed_T()
        self.assertEqual(nPoints, len(self.sim.grid))
        self.assertArrayNear(Tfixed, self.sim.T)

        self.solve_mix()
        data = np.empty((self.sim.flame.n_points, self.gas.n_species + 4))
        data[:,0] = self.sim.grid
        data[:,1] = self.sim.u
        data[:,2] = self.sim.V
        data[:,3] = self.sim.T
        data[:,4:] = self.sim.Y.T

        if saveReference:
            np.savetxt(referenceFile, data, '%11.6e', ', ')
        else:
            bad = utilities.compareProfiles(referenceFile, data,
                                            rtol=1e-2, atol=1e-8, xtol=1e-2)
            self.assertFalse(bad, bad)

    def test_auto(self, saveReference=False):
        referenceFile = pjoin(self.test_data_dir, 'DiffusionFlameTest-h2-auto.csv')
        self.create_sim(p=ct.one_atm, mdot_fuel=2, mdot_ox=3)

        nPoints = []
        timesteps = []

        def steady_func(x):
            nPoints.append(len(self.sim.T))
            return 0

        def time_step_func(dt):
            timesteps.append(dt)
            self.assertGreater(dt, 0)
            return 0

        self.sim.set_steady_callback(steady_func)
        self.sim.set_time_step_callback(time_step_func)

        self.sim.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.12, prune=0.0)
        self.sim.solve(loglevel=0, auto=True)

        self.assertNotEqual(len(nPoints), 0)
        self.assertNotEqual(len(timesteps), 0)

        data = np.empty((self.sim.flame.n_points, self.gas.n_species + 4))
        data[:,0] = self.sim.grid
        data[:,1] = self.sim.u
        data[:,2] = self.sim.V
        data[:,3] = self.sim.T
        data[:,4:] = self.sim.Y.T

        if saveReference:
            np.savetxt(referenceFile, data, '%11.6e', ', ')
        else:
            bad = utilities.compareProfiles(referenceFile, data,
                                            rtol=1e-2, atol=1e-8, xtol=1e-2)
            self.assertFalse(bad, bad)

    def run_extinction(self, mdot_fuel, mdot_ox, T_ox, width, P):
        self.create_sim(fuel='H2:1.0', oxidizer='O2:1.0', p=ct.one_atm*P,
                        mdot_fuel=mdot_fuel, mdot_ox=mdot_ox, width=width)
        self.sim.solve(loglevel=0, auto=True)
        self.assertFalse(self.sim.extinct())

    def test_extinction_case1(self):
        self.run_extinction(mdot_fuel=0.5, mdot_ox=3.0, T_ox=300, width=0.018, P=1.0)

    def test_extinction_case2(self):
        self.run_extinction(mdot_fuel=0.5, mdot_ox=1.0, T_ox=300, width=0.01, P=5.0)

    def test_extinction_case3(self):
        self.run_extinction(mdot_fuel=1.0, mdot_ox=0.5, T_ox=500, width=0.02, P=5.0)

    def test_extinction_case4(self):
        self.run_extinction(mdot_fuel=1.0, mdot_ox=3.0, T_ox=400, width=0.05, P=2.0)

    def test_extinction_case5(self):
        self.run_extinction(mdot_fuel=1.0, mdot_ox=3.0, T_ox=300, width=0.1, P=1.0)

    def test_extinction_case6(self):
        self.run_extinction(mdot_fuel=0.5, mdot_ox=0.5, T_ox=600, width=0.2, P=0.05)

    def test_extinction_case7(self):
        self.run_extinction(mdot_fuel=0.2, mdot_ox=2.0, T_ox=600, width=0.2, P=0.05)

    def test_mixture_averaged_rad(self, saveReference=False):
        referenceFile = pjoin(self.test_data_dir, 'DiffusionFlameTest-h2-mix-rad.csv')
        self.create_sim(p=ct.one_atm)

        nPoints = len(self.sim.grid)
        Tfixed = self.sim.T
        self.solve_fixed_T()
        self.assertEqual(nPoints, len(self.sim.grid))
        self.assertArrayNear(Tfixed, self.sim.T)
        self.assertFalse(self.sim.radiation_enabled)
        self.sim.radiation_enabled = True
        self.assertTrue(self.sim.radiation_enabled)
        self.sim.set_boundary_emissivities(0.25,0.15)

        self.solve_mix()
        data = np.empty((self.sim.flame.n_points, self.gas.n_species + 4))
        data[:,0] = self.sim.grid
        data[:,1] = self.sim.u
        data[:,2] = self.sim.V
        data[:,3] = self.sim.T
        data[:,4:] = self.sim.Y.T

        if saveReference:
            np.savetxt(referenceFile, data, '%11.6e', ', ')
        else:
            bad = utilities.compareProfiles(referenceFile, data,
                                            rtol=1e-2, atol=1e-8, xtol=1e-2)
            self.assertFalse(bad, bad)

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

        self.assertLessEqual(a_mean, a_max)
        self.assertLessEqual(a_pf_fuel, a_max)
        self.assertLessEqual(a_pf_oxidizer, a_max)
        self.assertLessEqual(a_stoich1, a_max)
        self.assertEqual(a_stoich1, a_stoich2)

        with self.assertRaises(ValueError):
            self.sim.strain_rate('bad_keyword')
        with self.assertRaises(KeyError): # missing 'fuel'
            self.sim.strain_rate('stoichiometric')
        with self.assertRaises(KeyError): # missing 'stoich'
            self.sim.strain_rate('stoichiometric', fuel='H2', oxidizer='H2O2')

    def test_mixture_fraction(self):
        self.create_sim(p=ct.one_atm)
        Z = self.sim.mixture_fraction('H')
        self.assertNear(Z[0], 1.0)
        self.assertNear(Z[-1], 0.0)
        self.assertTrue(all(Z >= 0))
        self.assertTrue(all(Z <= 1.0))


class TestCounterflowPremixedFlame(utilities.CanteraTest):
    # Note: to re-create the reference file:
    # (1) set PYTHONPATH to build/python2 or build/python3.
    # (2) Start Python in the test/work directory and run:
    #     >>> import cantera.test
    #     >>> t = cantera.test.test_onedim.TestCounterflowPremixedFlame("test_mixture_averaged")
    #     >>> t.test_mixture_averaged(True)

    def test_mixture_averaged(self, saveReference=False):
        T_in = 373.0  # inlet temperature
        comp = 'H2:1.6, O2:1, AR:7'  # premixed gas composition

        gas = ct.Solution('h2o2.xml')
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
        self.assertFalse(sim.radiation_enabled)
        sim.solve(loglevel=0, refine_grid=True)

        data = np.empty((sim.flame.n_points, gas.n_species + 4))
        data[:,0] = sim.grid
        data[:,1] = sim.u
        data[:,2] = sim.V
        data[:,3] = sim.T
        data[:,4:] = sim.Y.T

        referenceFile = pjoin(self.test_data_dir, 'CounterflowPremixedFlame-h2-mix.csv')
        if saveReference:
            np.savetxt(referenceFile, data, '%11.6e', ', ')
        else:
            bad = utilities.compareProfiles(referenceFile, data,
                                            rtol=1e-2, atol=1e-8, xtol=1e-2)
            self.assertFalse(bad, bad)

    def run_case(self, phi, T, width, P):
        gas = ct.Solution('h2o2.xml')
        gas.TPX = T, P * ct.one_atm, {'H2':phi, 'O2':0.5, 'AR':2}
        sim = ct.CounterflowPremixedFlame(gas=gas, width=width)
        sim.reactants.mdot = 10 * gas.density
        sim.products.mdot = 5 * gas.density
        sim.set_refine_criteria(ratio=6, slope=0.7, curve=0.8, prune=0.4)
        sim.solve(loglevel=0, auto=True)
        self.assertTrue(all(sim.T >= T - 1e-3))
        self.assertTrue(all(sim.V >= -1e-9))
        return sim

    def test_solve_case1(self):
        self.run_case(phi=0.4, T=400, width=0.05, P=10.0)

    def test_solve_case2(self):
        self.run_case(phi=0.5, T=500, width=0.03, P=2.0)

    def test_solve_case3(self):
        self.run_case(phi=0.7, T=300, width=0.05, P=2.0)

    def test_solve_case4(self):
        self.run_case(phi=1.5, T=400, width=0.03, P=0.02)

    def test_solve_case5(self):
        self.run_case(phi=2.0, T=300, width=0.2, P=0.2)


class TestBurnerFlame(utilities.CanteraTest):
    def solve(self, phi, T, width, P):
        gas = ct.Solution('h2o2.xml')
        gas.TPX = T, ct.one_atm*P, {'H2':phi, 'O2':0.5, 'AR':1.5}
        sim = ct.BurnerFlame(gas=gas, width=width)
        sim.burner.mdot = gas.density * 0.15
        sim.solve(loglevel=0, auto=True)
        self.assertGreater(sim.T[1], T)

    def test_case1(self):
        self.solve(phi=0.5, T=500, width=2.0, P=0.1)

    def test_case2(self):
        self.solve(phi=2.0, T=400, width=0.05, P=1.0)

    def test_case3(self):
        self.solve(phi=1.7, T=300, width=0.05, P=1.0)

    def test_case4(self):
        self.solve(phi=0.5, T=300, width=1.0, P=5.0)

    def test_case5(self):
        self.solve(phi=1.0, T=400, width=0.2, P=0.01)

    def test_fixed_temp(self):
        gas = ct.Solution('h2o2.xml')
        gas.TPX = 400, 2*ct.one_atm, {'H2':0.7, 'O2':0.5, 'AR':1.5}
        sim = ct.BurnerFlame(gas=gas, width=0.05)
        sim.burner.mdot = gas.density * 0.15
        sim.flame.set_fixed_temp_profile([0, 0.1, 0.9, 1],
                                         [400, 1100, 1100, 500])

        sim.energy_enabled = False
        sim.solve(loglevel=0, refine_grid=True)
        self.assertNear(sim.T[0], 400)
        self.assertNear(sim.T[-1], 500)
        self.assertNear(max(sim.T), 1100)

    def test_blowoff(self):
        gas = ct.Solution('h2o2.cti')
        gas.set_equivalence_ratio(0.4, 'H2', 'O2:1.0, AR:5')
        gas.TP = 300, ct.one_atm
        sim = ct.BurnerFlame(gas=gas, width=0.1)
        sim.burner.mdot = 1.2
        sim.set_refine_criteria(ratio=3, slope=0.3, curve=0.5, prune=0)
        sim.solve(loglevel=0, auto=True)
        # nonreacting solution
        self.assertNear(sim.T[-1], sim.T[0], 1e-6)
        self.assertNear(sim.u[-1], sim.u[0], 1e-6)
        self.assertArrayNear(sim.Y[:,0], sim.Y[:,-1], 1e-6, atol=1e-6)


class TestImpingingJet(utilities.CanteraTest):
    def run_reacting_surface(self, xch4, tsurf, mdot, width):
        # Simplified version of the example 'catalytic_combustion.py'
        gas = ct.Solution('ptcombust-simple.cti', 'gas')
        surf_phase = ct.Interface('ptcombust-simple.cti',
                                  'Pt_surf', [gas])

        tinlet = 300.0  # inlet temperature
        comp = {'CH4': xch4, 'O2':0.21, 'N2':0.79}
        gas.TPX = tinlet, ct.one_atm, comp
        surf_phase.TP = tsurf, ct.one_atm

        # integrate the coverage equations holding the gas composition fixed
        # to generate a good starting estimate for the coverages.
        surf_phase.advance_coverages(1.0)

        sim = ct.ImpingingJet(gas=gas, width=width, surface=surf_phase)
        sim.set_refine_criteria(10.0, 0.3, 0.4, 0.0)

        sim.inlet.mdot = mdot
        sim.inlet.T = tinlet
        sim.inlet.X = comp
        sim.surface.T = tsurf

        sim.solve(loglevel=0, auto=True)

        self.assertTrue(all(np.diff(sim.T) > 0))
        self.assertTrue(all(np.diff(sim.Y[gas.species_index('CH4')]) < 0))
        self.assertTrue(all(np.diff(sim.Y[gas.species_index('CO2')]) > 0))

    def test_reacting_surface_case1(self):
        self.run_reacting_surface(xch4=0.095, tsurf=900.0, mdot=0.06, width=0.1)

    def test_reacting_surface_case2(self):
        self.run_reacting_surface(xch4=0.07, tsurf=1200.0, mdot=0.2, width=0.05)

    def test_reacting_surface_case3(self):
        self.run_reacting_surface(xch4=0.2, tsurf=800.0, mdot=0.1, width=0.2)


class TestTwinFlame(utilities.CanteraTest):
    def solve(self, phi, T, width, P):
        gas = ct.Solution('h2o2.xml')
        gas.TP = T, ct.one_atm
        gas.set_equivalence_ratio(phi, 'H2', 'O2:1.0, AR:4.0')
        sim = ct.CounterflowTwinPremixedFlame(gas=gas, width=width)
        sim.set_refine_criteria(ratio=5, slope=0.6, curve=0.8, prune=0.1)
        axial_velocity = 2.0
        sim.reactants.mdot = gas.density * axial_velocity
        sim.solve(loglevel=0, auto=True)
        self.assertGreater(sim.T[-1], T + 100)

    def test_case1(self):
        self.solve(phi=0.4, T=300, width=0.05, P=0.1)


class TestIonFlame(utilities.CanteraTest):
    def test_ion_profile(self):
        reactants = 'CH4:0.216, O2:2'
        p = ct.one_atm
        Tin = 300
        width = 0.03

        # IdealGasMix object used to compute mixture properties
        self.gas = ct.Solution('ch4_ion.cti')
        self.gas.TPX = Tin, p, reactants
        self.sim = ct.IonFlame(self.gas, width=width)
        self.sim.set_refine_criteria(ratio=4, slope=0.8, curve=1.0)
        # Ionized species may require tighter absolute tolerances
        self.sim.flame.set_steady_tolerances(Y=(1e-4, 1e-12))
        self.sim.transport_model = 'Ion'

        # stage one
        self.sim.solve(loglevel=0, auto=True)

        #stage two
        self.sim.solve(loglevel=0, stage=2, enable_energy=True)

        # Regression test
        self.assertNear(max(self.sim.E), 131.9956, 1e-3)

class TestOneDimNumCont(utilities.CanteraTest):
    def get_C():
        #start the clock
        timestart = time.clock()
        index_array = np.zeros((len(composition)), dtype=np.int)
        for i,name in enumerate(composition):
            index_array[i] = gas.species_index(name)
        C_array = np.zeros(len(f.grid), dtype=np.double)
        for j in range(0,len(f.grid)-1):
            C=0.0
            for i in index_array:
                C += f.Y[i,j]
            if j != 0 and C < max(C_array) and C > 0.5: 
                break
            C_array[j] = C
        timeend = time.clock()

        return C_array   

    def get_mixture_fraction():
        #start the clock
        timestart = time.clock()

        # create arrays needed
        mixture_fraction_array =  np.zeros((len(f.grid)), dtype=np.double)
        atom_end = np.zeros((len(z_composition), len(f.grid)), dtype=np.double)
        denominator = np.zeros((len(z_composition)), dtype=np.double)

        # convert dictionary to array for faster execution
        denominator_Z = 0.0
        for i, atom in enumerate(z_composition):
            atom_end[i] = f.elemental_mass_fraction(atom) 
            denominator_Z += f.elemental_mass_fraction(atom)[0] - f.elemental_mass_fraction(atom)[-1]   
                 
        # calculate mixture fraction
        for j in range(0,len(f.grid)):              
            numerator_Z = 0.0
            for i in range(len(z_composition)):
                numerator_Z += atom_end[i,j] - atom_end[i, -1]
            mixture_fraction = numerator_Z / (denominator_Z + 10**-30)
            mixture_fraction_array[j] = mixture_fraction
        timeend = time.clock()

        return mixture_fraction_array

    
    def iloc(percent):
        Ttarget = percent * np.max(f.T)
        for i in range(0,len(f.grid),1):
            if f.T[i] > Ttarget:
                z0 = i
                break
        for i in range(len(f.grid)-1,0,-1):
            if f.T[i] > Ttarget:
                z1 = i
                break
        return z0, z1

    def stepControl(ieval, N_opt, dT, MAXdT, MINdT):
        failed = False
        xi = float(N_opt) / (ieval + 1)
        if xi < 0.5:
            xi = 0.5
        elif xi > 2.0:
            xi = 2.0
        dT *= xi
        if dT > MAXdT:  
            dT = MAXdT
            failed = True
        elif dT < MINdT:
            dT = MINdT  
            failed = True
        return failed, dT
    
    def setFuelBoundary():
        reactants = ''
        for i, molecules in enumerate(fuel):
            reactants += str(fuel[i]) + ':' + str(fuel_fraction[i])
            if len(fuel) > 1 and i < len(fuel) - 1:
                reactants += ','
        if mole_fractions:
            f.fuel_inlet.X = reactants
            gas.TPX = Tfuel, target_pressure, reactants
        else:
            f.fuel_inlet.Y = reactants
            gas.TPY = Tfuel, target_pressure, reactants
    
    def setOxidBoundary():
        reactants = ''
        for i, molecules in enumerate(oxidizer):
            reactants += str(oxidizer[i]) + ':' + str(oxidizer_fraction[i])
            if len(oxidizer) > 1 and i < len(oxidizer) - 1:
                reactants += ','
        if mole_fractions:
            f.oxidizer_inlet.X = reactants
            gas.TPX = Toxidizer, target_pressure, reactants
        else:
            f.oxidizer_inlet.Y = reactants
            gas.TPY = Toxidizer, target_pressure, reactants

    def molarFuelAirRatio():
        setFuelBoundary()
        setOxidBoundary()
        x = 0.0
        y = 0.0
        for j in range(0, gas.n_elements):
            if gas.element_names[j] == 'C':
                for i in range(gas.n_species):
                    x += f.fuel_inlet.X[i] * gas.n_atoms(gas.species_names[i], gas.element_names[j])
            if gas.element_names[j] == 'H':
                for i in range(gas.n_species):
                    y += f.fuel_inlet.X[i] * gas.n_atoms(gas.species_names[i], gas.element_names[j])
        stoic_num_oxid_moles = x + y / 4.0
        
        return stoic_num_oxid_moles

    def stoicFuelAirRatio():
        """
        Here we are assuming that fuel is only on the left and oxidizer on the right.
        This formulation is only for diffusion flamelets
        """
        setFuelBoundary()
        setOxidBoundary()

        numerator = np.dot(f.fuel_inlet.X, gas.molecular_weights)           
        denominator = np.dot (f.oxidizer_inlet.X, gas.molecular_weights)
        denominator *= molarFuelAirRatio() * np.sum(f.oxidizer_inlet.X) / f.oxidizer_inlet.X[f.flame.component_index('O2')-5]

        return numerator /  denominator


    def limits():
        file =  open("equilibrium_table.txt", "w")
        Zstoic = stoicFuelAirRatio()
        Zfuel = 1.0
        Zoxid = 0.0
        Z = np.linspace(Zoxid, Zfuel, 101)
        
        modified_fuel_fraction = np.zeros(len(fuel_fraction), dtype= 'float')
        modified_fuel_fraction[:] = fuel_fraction
        for n in range(len(Z)):
            if Z[n] == 0.0:
            #maximum value of progress variable (C) as function of mixture fractiton(Z)
                file.write("%f %f\n" % (Zoxid, 0.0)) 

            elif Z[n] > 0.0 and Z[n] < 1.0:
                phi = np.asscalar(Z[n]) / (1.0 -  np.asscalar(Z[n])) / Zstoic
                
                reactants = ''
                shared_species = []
                # if species is in both the fuel and oxidizer streams Cantera cannot take duplicates
                for j, fuel_species in enumerate(fuel):
                    for i, oxidizer_species in enumerate(oxidizer):
                        if oxidizer_species == fuel_species:
                            shared_species.append(fuel_species)
                            temp = molarFuelAirRatio() * oxidizer_fraction[i] / oxidizer_fraction[0] / phi                             
                            modified_fuel_fraction[j] = fuel_fraction[j] + temp
                for i, molecules in enumerate(fuel):
                    reactants += str(fuel[i]) + ':' + str(modified_fuel_fraction[i])
                    if len(fuel) > 1 and i < len(fuel) - 1:
                        reactants += ','                    
                    if len(shared_species) < len(oxidizer):
                        reactants += ','
                for i, molecules in enumerate(oxidizer):
                    if molecules not in shared_species:
                        reactants += str(oxidizer[i]) + ':' \
                        + str( molarFuelAirRatio() * oxidizer_fraction[i] / \
                        oxidizer_fraction[0]/ phi)
                        if len(oxidizer) > 1 and i < len(oxidizer) - 1:
                            reactants += ','    
                if mole_fractions:
                    gas.TPX = 0.5*(Tfuel+Toxidizer), pressure, reactants
                else:
                    gas.TPY = 0.5*(Tfuel+Toxidizer), pressure, reactants
                gas.equilibrate('HP')
                
                C_limit = 0.0
                for species in composition:
                    C_limit += gas.Y[gas.species_index(species)]
                file.write("%f %f\n" % (Z[n], C_limit)) 

            elif Z[n] == 1.0:
            #mixture fraction(Z) vs. maximum value of progress variable (C)
                file.write("%f %f\n" % (Zfuel, 0.0)) 

        file.close()  
    def test_diffusion_flame_cont(self):

        # start of main, set plot mode to interactive
        init_plot = True
        PrintFlush("Start flamelet generation of FPV type")

        tempdir = tempfile.TemporaryDirectory()
        data_directory = tempdir.name + "\\"

        UnityLewisNumber = True
        gas = ct.Solution('C:\cantera-2.3.0\data\inputs\h2o2.cti')
        reaction_mechanism = 'h2o2'
        pressure = 101325
        Tfuel =  300.
        Toxidizer = 300.
        mole_fractions = True
        fuel = ['H2']
        fuel_fraction = [1.0]
        oxidizer = ['O2']
        oxidizer_fraction = [1.0]
        init_fuel_flux = 0.1
        initial_grid = np.linspace(0, 0.2, 21) 
        tol_ss = [1.0e-10,1.0e-20] 
        tol_ts = [1.0e-10,1.0e-20] 
        ss_age = 10  
        ts_age = 10  
        stepsize = 1e-5 
        nstep = [2,5,10,20] 
        loglevel = 0 
        refine = True  
        temperature_limit_extinction = 500
        ratio = 5.0
        slope = 0.2
        curve = 0.2
        prune = 0.05
        z_composition = ['H']
        composition = ['H2O'] 
        z_points = 21 
        StrainEq = False
        precalculations = 1
        UserStrainFactor = 3.0
        UpperBranchHomotopicIter = 500
        UserMaxTempPercentage = 0.90
        FuelSign = -1
        OxidSign = -1
        OnePointControl = 0
        TwoPointControl = 1
        dT = 20.
        FailedSolutionsToTerminate = 10
        maxNumberOfFlames = 10000
        N_opt = 20
        MAXdT = 20.
        MINdT = 1.0
        OxidizerMassFlux = 0.0
        transport_model = 'Mix'

        # Exponents for the initial solution variation with changes in strain rate
        # Taken from Fiala and Sattelmayer (2014)
        exp_d_a = - 1. / 2.
        exp_u_a = 1. / 2.
        exp_V_a = 1.
        exp_lam_a = 2.
        exp_mdot_a = 1. / 2.

        #Flame object
        f = ct.CounterflowDiffusionFlame(gas, initial_grid)
        f.set_max_grid_points(1, 250)
        f.set_max_jac_age(ss_age,ts_age)  #for some reason this was commented
        f.set_time_step(stepsize,nstep)  # this was not here 
        f.flame.set_steady_tolerances(default=tol_ss)
        f.flame.set_transient_tolerances(default=tol_ts)
        f.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)
        f.set_grid_min(1e-5) # we can keep this input to python script only
        f.transport_model = transport_model

        f.fuel_inlet.T = Tfuel
        f.oxidizer_inlet.T = Toxidizer
        f.fuel_inlet.spread_rate = 0.0
        f.oxidizer_inlet.spread_rate = 0.0
        target_pressure = pressure
        f.fuel_inlet.mdot = init_fuel_flux
        refine_grid = True # There is no need to add this one again
        if OnePointControl == True:
            oxidizer_flux = OxidizerMassFlux
        else:
            oxidizer_flux = f.oxidizer_inlet.mdot # oxidizer flux initialization

        #############################
        file_name = 'restart.xml'  
        restart = False

        #############################################
        middleBranch = False
        lowerBranch = False
        turning = False
        flameControl = False
        # loop trying to solve for current temperatures
        # allows retry if a solution fails
        limits()

        Continuation = True
        stepSizeControl = False # Active/Deactive step size control

        #first time in solve loop
        T_max = []
        a_max = []
        m_dots = []

        for i in range(precalculations):
            PrintFlush(" Precalculations")
            target_pressure = pressure * (i + 1) / precalculations
            f.P = target_pressure

            reactants = ''
            for i, molecules in enumerate(fuel):
                reactants += str(fuel[i]) + ':' + str(fuel_fraction[i])
                if len(fuel) > 1 and i < len(fuel) - 1:
                    reactants += ','
            if mole_fractions:
                f.fuel_inlet.X = reactants
                gas.TPX = Tfuel, target_pressure, reactants
            else:
                f.fuel_inlet.Y = reactants
                gas.TPY = Tfuel, target_pressure, reactants
            fuel_density = gas.density_mass

            reactants = ''
            for i, molecules in enumerate(oxidizer):
                reactants += str(oxidizer[i]) + ':' + str(oxidizer_fraction[i])
                if len(oxidizer) > 1 and i < len(oxidizer) - 1:
                    reactants += ','
            if mole_fractions:
                f.oxidizer_inlet.X = reactants
                gas.TPX = Toxidizer, target_pressure, reactants
            else:
                f.oxidizer_inlet.Y = reactants
                gas.TPY = Toxidizer, target_pressure, reactants
            oxidizer_density = gas.density_mass

            if i > 0:
                f.fuel_inlet.mdot *= 1.0 # pow((i + 1)/i , 1.25)
            f.oxidizer_inlet.mdot = f.fuel_inlet.mdot * np.sqrt(oxidizer_density/fuel_density)

            ###########################################################################
            f.set_flame_control(1, StrainEq, UnityLewisNumber, False, False, -2000, -20, -2000, -20, False)
            if i == 0:
                f.solve(loglevel=loglevel, refine_grid=refine_grid, auto=True)
            else: 
                f.solve(loglevel=loglevel, refine_grid=refine_grid, auto=False)
            PrintFlush("  Pressure      Tmax Fuel_mdot")
            PrintFlush("%10.1f %10.1f %10.2f" % (f.P, max(f.T), f.fuel_inlet.mdot))
        
        f.save(data_directory + file_name, name='solution',
        description='Cantera version ' + ct.__version__ +
        ', reaction mechanism ' + reaction_mechanism, loglevel=loglevel)  
        
        mixture_fraction_array = get_mixture_fraction()    

        T_max.append(max(f.T))
        a_max.append(f.strain_rate('max'))
        m_dots.append(f.fuel_inlet.mdot)

        C_array = get_C()
        max_C = np.max(C_array)
        if len(T_max) == 1:
            global_max_C = max_C
            iFlamesInUpperBranch = 1
            iFlamesInMiddleBranch = 0
            iFlamesInLowerBranch = 0
            n = 1

        mixture_fraction_array = get_mixture_fraction()

        iFailedSolutions = 0

        iflag = 0
        f.clear_stats()

        while Continuation:
            stepFailed = False
            Tfuel_j, Toxid_j = iloc(UserMaxTempPercentage)
            Tfuel_fixed = f.T[Tfuel_j] + FuelSign * dT
            Toxid_fixed = f.T[Toxid_j] + OxidSign * dT
        
            if restart == False and iFlamesInUpperBranch > UpperBranchHomotopicIter:
                flameControl = True
        
            if n > 50 and turning == False:
                if np.sign((a_max[-1] - a_max[-2]) * (a_max[-2] - a_max[-3])) < 0:
                    middleBranch = True

            if middleBranch == True and iFlamesInMiddleBranch > 2000:
                strain_factor = a_max[-1] / a_max[-2]
                f.flame.grid *= strain_factor ** exp_d_a   
                                    
            if flameControl == False and middleBranch == False:
                strain_factor = (n + UserStrainFactor ) / n
                # Create an initial guess based on the previous solution
                # Update grid
                f.flame.grid *= strain_factor ** exp_d_a
                normalized_grid = f.grid / (f.grid[-1] - f.grid[0])
                # Update mass fluxes
                f.fuel_inlet.mdot *= strain_factor ** exp_mdot_a
                f.oxidizer_inlet.mdot *= strain_factor ** exp_mdot_a
                # Update velocities
                f.set_profile('u', normalized_grid, f.u * strain_factor ** exp_u_a)
                f.set_profile('V', normalized_grid, f.V * strain_factor ** exp_V_a)
                # Update pressure curvature
                f.set_profile('lambda', normalized_grid, f.L * strain_factor ** exp_lam_a)

            if flameControl == True:
                oxidizer_flux = f.oxidizer_inlet.mdot
                if OnePointControl == True:
                    f.oxidizer_inlet.mdot = oxidizer_flux
                f.set_flame_control(1, StrainEq, UnityLewisNumber, OnePointControl, TwoPointControl, \
                                    Tfuel_fixed, Tfuel_j, Toxid_fixed, Toxid_j, True)

            try:
                success = True
                f.solve(loglevel=loglevel, refine_grid=refine_grid, auto=False)  
        
            except Exception:
                success = False

            if success == True:
                if f.extinct() == False:
                    if n % 1 == 0:
                        T_max.append(max(f.T))
                        a_max.append(f.strain_rate('max'))
                        m_dots.append(f.fuel_inlet.mdot)
                              
                    if n % 1 == 0 and T_max[-1] < T_max[-2]:
                        f.save(data_directory + file_name, name='solution',
                        description='Cantera version ' + ct.__version__ +
                        ', reaction mechanism ' + reaction_mechanism, loglevel=0)

                    C_array = get_C()
                    max_C = np.max(C_array)
                    mixture_fraction_array = get_mixture_fraction()
                   
                    ieval = f.eval_count_stats[-1]
                    if middleBranch == False and lowerBranch == False:
                        branch = 'upper branch'
                    elif middleBranch == True and lowerBranch == False:
                        branch = 'middle branch'
                    elif middleBranch == False and lowerBranch == True:
                        branch = 'lower branch'
                    if n ==1 or n % 10 == 0:
                        PrintFlush("     #n      Tmax       Amax         dT     branch   FlameControl?   # gridpoints")
                    PrintFlush("%6d %10.1f %10.1f %10.2f %10s %10s %12d" \
                            % (n, max(f.T), f.strain_rate('max'), dT, branch.split(' ')[0], str(flameControl),  f.flame.n_points))
                    if middleBranch == False:
                        iFlamesInUpperBranch += 1
                    if middleBranch == True:
                        iFlamesInMiddleBranch += 1     
                        if T_max[-1] > T_max[-2]:
                            PrintFlush("***********Maximum temperature higher than that of previous solution******")
                            T_max.pop(-1)
                            a_max.pop(-1)
                            m_dots.pop(-1)
                            Continuation = False
                        else:
                             continue

                    if iFailedSolutions > FailedSolutionsToTerminate:
                        Continuation = False
                else:          
                    if turning == False:
                        iflag = 1
                        f.restore(data_directory + file_name, name='solution', loglevel=loglevel)
                        flameControl = True
            else:
                iflag = 1
                iFailedSolutions += 1
                stepSizeControl = True     
                flameControl = True
                stepFailed, dT = stepControl(ieval, N_opt, dT, MAXdT, MINdT)
                PrintFlush('Calculation failed; the new value of dT is %f' % dT)
                f.restore(data_directory + file_name, name='solution', loglevel=loglevel)
                if stepFailed == False:
                    continue
                else:
                    Continuation = False

            if Continuation == True:
                # check for maximum number of flames
                n += 1
                if n > maxNumberOfFlames:
                    continuation = False
                    retrySolve = False;
                    break   
            
                # update step control 
                if stepSizeControl == True and flameControl == True:
                    stepFailed, dT = stepControl(ieval, N_opt, dT, MAXdT, MINdT)

            else:
                restart = False
                retrySolve = False
                break

        # Calculating C = 0.0
        f.fuel_inlet.mdot *= 100
        f.oxidizer_inlet.mdot *=100
        f.set_max_grid_points(1, 101)
        for j in range(len(f.grid)):
            Tinter = min(Tfuel,Toxidizer)
            f.set_value(1, 'T', j, Tinter)    
        f.set_flame_control(1, StrainEq, UnityLewisNumber, False, False, -2000, -20, -2000, -20, False)
        f.solve(loglevel=loglevel, refine_grid=False, auto=False)
        T_max.append(max(f.T))
        a_max.append(f.strain_rate('max'))
        m_dots.append(f.fuel_inlet.mdot)
        C_array = get_C()
        max_C = np.abs(np.max(C_array))
        mixture_fraction_array = get_mixture_fraction()
        PrintFlush("%6d %10.1f %10.1f %10.2f %10s %10s %12d" \
                % (n, max(f.T), f.strain_rate('max'), dT, "lower", "False", f.flame.n_points))    
>>>>>>> one dimensional continuation python test
