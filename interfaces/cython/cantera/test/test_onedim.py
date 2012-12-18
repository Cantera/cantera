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

    def solve_mix(self, ratio=3.0, slope=0.2, curve=0.1, prune=0.0):
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

        self.solve_mix(slope=0.2, curve=0.1)
        T1 = self.sim.T[-1]
        X1 = self.sim.X[:,-1]

        self.solve_mix(slope=0.1, curve=0.05)
        T2 = self.sim.T[-1]
        X2 = self.sim.X[:,-1]

        self.solve_mix(slope=0.05, curve=0.025)
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

        self.sim.setGasState(0)
        rhou = self.gas.density * self.sim.u[0]

        for j in range(self.sim.flame.nPoints):
            self.sim.setGasState(j)

            # Check continuity
            self.assertNear(self.gas.density * self.sim.u[j], rhou, 1e-4)

    def test_multicomponent(self):
        reactants= 'H2:1.1, O2:1, AR:5.5'
        p = ct.OneAtm
        Tin = 300

        self.create_sim(p, Tin, reactants)
        self.solve_fixed_T()
        self.solve_mix(ratio=4, slope=0.4, curve=0.2)
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
        self.solve_mix(slope=0.1, curve=0.05, prune=0.0)
        N1 = len(self.sim.grid)

        self.solve_mix(slope=0.4, curve=0.2, prune=0.05)
        N2 = len(self.sim.grid)

        self.assertLess(N2, N1)

        # TODO: check that the solution is actually correct (i.e. that the
        # residual satisfies the error tolerances) on the new grid.

    def test_save_restore(self):
        reactants= 'H2:1.1, O2:1, AR:5'
        p = ct.OneAtm
        Tin = 300

        self.create_sim(p, Tin, reactants)
        self.solve_fixed_T()
        filename = 'onedim-fixed-T.xml'
        if os.path.exists(filename):
            os.remove(filename)

        Y1 = self.sim.Y
        u1 = self.sim.u
        V1 = self.sim.V

        self.sim.save(filename, 'test')

        self.create_sim(p, Tin, reactants)
        self.sim.restore(filename, 'test')
        Y2 = self.sim.Y
        u2 = self.sim.u
        V2 = self.sim.V

        self.assertArrayNear(Y1, Y2)
        self.assertArrayNear(u1, u2)
        self.assertArrayNear(V1, V2)
