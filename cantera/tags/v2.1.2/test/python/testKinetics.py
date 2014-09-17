import utilities
import Cantera as ct

class ExplicitForwardOrderTest(utilities.CanteraTest):
    def setUp(self):
        self.gas = ct.IdealGasMix('../data/explicit-forward-order.xml')
        self.gas.set(T=800, P=101325, X=[0.01, 0.90, 0.02, 0.03, 0.04])

    def test_irreversibility(self):
        # Reactions are irreversible
        Rr = self.gas.revRateConstants()
        self.assertEqual(Rr[0], 0.0)
        self.assertEqual(Rr[1], 0.0)

    def test_rateConstants(self):
        # species order: [H, AR, R1A, R1B, P1]

        C = self.gas.moleFractions() * self.gas.molarDensity()
        Rf = self.gas.fwdRatesOfProgress()
        kf = self.gas.fwdRateConstants()
        self.assertNear(Rf[0], kf[0] * C[2]**1.5 * C[3]**0.5)
        self.assertNear(Rf[1], kf[1] * C[0]**1.0 * C[4]**0.2)

    def test_ratio1(self):
        rop1 = self.gas.fwdRatesOfProgress()
        # Double concentration of H and R1A
        self.gas.set(T=800, P=101325, X=[0.02, 0.87, 0.04, 0.03, 0.04])
        rop2 = self.gas.fwdRatesOfProgress()
        ratio = rop2/rop1
        self.assertNear(ratio[0], 2**1.5) # order of R1A is 1.5
        self.assertNear(ratio[1], 2**1.0) # order of H is 1.0

    def test_ratio2(self):
        rop1 = self.gas.fwdRatesOfProgress()
        # Double concentration of P1 and R1B
        self.gas.set(T=800, P=101325, X=[0.01, 0.83, 0.02, 0.06, 0.08])
        rop2 = self.gas.fwdRatesOfProgress()
        ratio = rop2/rop1
        self.assertNear(ratio[0], 2**0.5) # order of R1B is 0.5
        self.assertNear(ratio[1], 2**0.2) # order of P1 is 1.0
