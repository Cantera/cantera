import cantera as ct
from . import utilities

class TestPureFluid(utilities.CanteraTest):
    def setUp(self):
        self.water = ct.Water()

    def test_critical_properties(self):
        self.assertNear(self.water.critical_pressure, 22.089e6)
        self.assertNear(self.water.critical_temperature, 647.286)
        self.assertNear(self.water.critical_density, 317.0)

    def test_set_state(self):
        self.water.PX = 101325, 0.5
        self.assertNear(self.water.P, 101325)
        self.assertNear(self.water.X, 0.5)

        self.water.TX = 500, 0.8
        self.assertNear(self.water.T, 500)
        self.assertNear(self.water.X, 0.8)

    def test_set_minmax(self):
        self.water.TP = self.water.min_temp, 101325
        self.assertNear(self.water.T, self.water.min_temp)

        self.water.TP = self.water.max_temp, 101325
        self.assertNear(self.water.T, self.water.max_temp)

    def check_fd_properties(self, T1, P1, T2, P2, tol):
        # Properties which are computed as finite differences
        self.water.TP = T1, P1
        cp1 = self.water.cp_mass
        cv1 = self.water.cv_mass
        k1 = self.water.isothermal_compressibility
        alpha1 = self.water.thermal_expansion_coeff

        self.water.TP = T2, P2
        cp2 = self.water.cp_mass
        cv2 = self.water.cv_mass
        k2 = self.water.isothermal_compressibility
        alpha2 = self.water.thermal_expansion_coeff

        self.assertNear(cp1, cp2, tol)
        self.assertNear(cv1, cv2, tol)
        self.assertNear(k1, k2, tol)
        self.assertNear(alpha1, alpha2, tol)

    def test_properties_near_min(self):
        self.check_fd_properties(self.water.min_temp*(1+1e-5), 101325,
                              self.water.min_temp*(1+1e-4), 101325, 1e-2)

    def test_properties_near_max(self):
        self.check_fd_properties(self.water.max_temp*(1-1e-5), 101325,
                              self.water.max_temp*(1-1e-4), 101325, 1e-2)
