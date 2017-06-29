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

    def test_soret_flag(self):
        self.create_sim(101325, 300, 'H2:1.0, O2:1.0')
        self.assertFalse(self.sim.soret_enabled)
        with self.assertRaises(ct.CanteraError):
            self.sim.soret_enabled = True
        self.sim.transport_model = 'Multi'
        self.sim.soret_enabled = True

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
