import cantera as ct
from . import utilities
import os


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

        surface = ct.ReactingSurface1D()
        surface.setKinetics(interface)


class TestFreeFlame(utilities.CanteraTest):
    def create_sim(self, p, Tin, reactants):

        initial_grid = [0.0, 0.001, 0.01, 0.02, 0.029, 0.03]  # m
        tol_ss = [1.0e-5, 1.0e-13]  # [rtol atol] for steady-state problem
        tol_ts = [1.0e-4, 1.0e-10]  # [rtol atol] for time stepping

        # IdealGasMix object used to compute mixture properties
        self.gas = ct.Solution('h2o2.xml')
        self.gas.TPX = Tin, p, reactants

        # Flame object
        self.sim = ct.FreeFlame(self.gas, initial_grid)
        self.sim.flame.setSteadyTolerances(default=tol_ss)
        self.sim.flame.setTransientTolerances(default=tol_ts)

        # Set properties of the upstream fuel-air mixture
        self.sim.inlet.T = Tin
        self.sim.inlet.X = reactants

    def solve_fixed_T(self):
        # Solve with the energy equation disabled
        self.sim.energyEnabled = False
        self.sim.setMaxJacAge(50, 50)
        self.sim.setTimeStep(1e-5, [2, 5, 10, 20])
        self.sim.solve(loglevel=0, refine_grid=False)

        self.assertFalse(self.sim.energyEnabled)

    def solve_mix(self, ratio=3.0, slope=0.3, curve=0.2, prune=0.0):
        # Solve with the energy equation enabled

        self.sim.setRefineCriteria(ratio=ratio, slope=slope, curve=curve, prune=prune)
        self.sim.energyEnabled = True
        self.sim.solve(loglevel=0, refine_grid=True)

        self.assertTrue(self.sim.energyEnabled)
        self.assertEqual(self.sim.transportModel, 'Mix')

    def solve_multi(self):
        self.sim.transportModel = 'Multi'
        self.sim.solve(loglevel=0, refine_grid=True)

        self.assertEqual(self.sim.transportModel, 'Multi')

    def test_converge_adiabatic(self):
        # Test that the adiabatic flame temperature and species profiles
        # converge to the correct equilibrium values as the grid is refined

        reactants= 'H2:1.1, O2:1, AR:5'
        p = ct.OneAtm
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

        for k in range(self.gas.nSpecies):
            self.assertLess(abs(X2[k]-Xad[k]), abs(X1[k]-Xad[k]))
            self.assertLess(abs(X3[k]-Xad[k]), abs(X2[k]-Xad[k]))

    def test_mixture_averaged(self):
        reactants= 'H2:1.1, O2:1, AR:5'
        p = ct.OneAtm
        Tin = 300

        self.create_sim(p, Tin, reactants)
        self.solve_fixed_T()
        self.solve_mix()

        self.gas.TPX = Tin, p, reactants
        self.gas.equilibrate('HP')

        rhou = self.sim.density[0] * self.sim.u[0]

        # Check continuity
        for rhou_j in self.sim.density * self.sim.u:
            self.assertNear(rhou_j, rhou, 1e-4)

    def test_multicomponent(self):
        reactants= 'H2:1.1, O2:1, AR:5.3'
        p = ct.OneAtm
        Tin = 300

        self.create_sim(p, Tin, reactants)
        self.solve_fixed_T()
        self.solve_mix(ratio=5, slope=0.5, curve=0.3)
        Su_mix = self.sim.u[0]

        # Multicomponent flame speed should be similar but not identical to
        # the mixture-averaged flame speed.
        self.solve_multi()
        Su_multi = self.sim.u[0]
        self.assertFalse(self.sim.soretEnabled)

        self.assertNear(Su_mix, Su_multi, 5e-2)
        self.assertNotEqual(Su_mix, Su_multi)

        # Flame speed with Soret effect should be similar but not identical to
        # the multicomponent flame speed
        self.sim.soretEnabled = True
        self.sim.solve(loglevel=0, refine_grid=True)
        self.assertTrue(self.sim.soretEnabled)
        Su_soret = self.sim.u[0]

        self.assertNear(Su_multi, Su_soret, 2e-1)
        self.assertNotEqual(Su_multi, Su_soret)

    def test_prune(self):
        reactants= 'H2:1.1, O2:1, AR:5'
        p = ct.OneAtm
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
        reactants= 'H2:1.1, O2:1, AR:5'
        p = 2 * ct.OneAtm
        Tin = 400

        self.create_sim(p, Tin, reactants)
        self.solve_fixed_T()
        filename = 'onedim-fixed-T.xml'
        if os.path.exists(filename):
            os.remove(filename)

        Y1 = self.sim.Y
        u1 = self.sim.u
        V1 = self.sim.V
        P1 = self.sim.P

        self.sim.save(filename, 'test', loglevel=0)

        self.create_sim(ct.OneAtm, Tin, reactants)
        self.sim.energyEnabled = False
        self.sim.restore(filename, 'test', loglevel=0)

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

        self.assertArrayNear(Y1, Y3, 1e-3)
        self.assertArrayNear(u1, u3, 1e-3)
        self.assertArrayNear(V1, V3, 1e-3)

    def test_array_properties(self):
        self.create_sim(ct.OneAtm, 300, 'H2:1.1, O2:1, AR:5')

        for attr in ct.FlameBase.__dict__:
            if isinstance(ct.FlameBase.__dict__[attr], property):
                getattr(self.sim, attr)
