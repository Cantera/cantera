import unittest
import numpy as np

import cantera as ct
from . import utilities

class TestKinetics(utilities.CanteraTest):
    def setUp(self):
        self.phase = ct.Solution('h2o2.xml')
        self.phase.X = [0.1, 1e-4, 1e-5, 0.2, 2e-4, 0.3, 1e-6, 5e-5, 0.4]
        self.phase.TP = 800, 2*ct.OneAtm

    def test_counts(self):
        self.assertEqual(self.phase.nReactions, 27)
        self.assertEqual(self.phase.nTotalSpecies, 9)
        self.assertEqual(self.phase.nPhases, 1)
        self.assertEqual(self.phase.reactionPhaseIndex, 0)

    def test_isReversible(self):
        for i in range(self.phase.nReactions):
            self.assertTrue(self.phase.isReversible(i))

    def test_multiplier(self):
        fwd_rates0 = self.phase.fwdRatesOfProgress
        rev_rates0 = self.phase.revRatesOfProgress

        self.phase.setMultiplier(2.0, 0)
        self.phase.setMultiplier(0.1, 6)

        fwd_rates1 = self.phase.fwdRatesOfProgress
        rev_rates1 = self.phase.revRatesOfProgress

        self.assertNear(2 * fwd_rates0[0], fwd_rates1[0])
        self.assertNear(0.1 * fwd_rates0[6], fwd_rates1[6])
        self.assertNear(2 * rev_rates0[0], rev_rates1[0])
        self.assertNear(0.1 * rev_rates0[6], rev_rates1[6])
        for i in range(self.phase.nReactions):
            if i not in (0,6):
                self.assertNear(fwd_rates0[i], fwd_rates1[i])
                self.assertNear(rev_rates0[i], rev_rates1[i])

        self.phase.setMultiplier(0.5)
        fwd_rates2 = self.phase.fwdRatesOfProgress
        rev_rates2 = self.phase.revRatesOfProgress
        self.assertArrayNear(0.5 * fwd_rates0, fwd_rates2)
        self.assertArrayNear(0.5 * rev_rates0, rev_rates2)

    def test_reactionType(self):
        self.assertNear(self.phase.reactionType(0), 2) # 3rd body
        self.assertNear(self.phase.reactionType(2), 1) # elementary
        self.assertNear(self.phase.reactionType(19), 4) # falloff

        self.assertRaises(ValueError, self.phase.reactionType, 33)
        self.assertRaises(ValueError, self.phase.reactionType, -2)

    def test_reactionEquations(self):
        self.assertEqual(self.phase.nReactions,
                         len(self.phase.reactionEquations()))
        self.assertEqual(self.phase.reactionEquation(16),
                         'H + H2O2 <=> HO2 + H2')

    def test_stoichCoeffs(self):
        nu_r = self.phase.reactantStoichCoeffs()
        nu_p = self.phase.productStoichCoeffs()

        def check_reactant(k, i, value):
            self.assertEqual(self.phase.reactantStoichCoeff(k,i), value)
            self.assertEqual(nu_r[k,i], value)

        def check_product(k, i, value):
            self.assertEqual(self.phase.productStoichCoeff(k,i), value)
            self.assertEqual(nu_p[k,i], value)

        # H + H2O2 <=> HO2 + H2
        check_reactant(1, 16, 1)
        check_reactant(7, 16, 1)
        check_reactant(6, 16, 0)
        check_reactant(0, 16, 0)

        check_product(1, 16, 0)
        check_product(7, 16, 0)
        check_product(6, 16, 1)
        check_product(0, 16, 1)

        # 2 O + M <=> O2 + M
        check_reactant(2, 0, 2)
        check_reactant(3, 0, 0)
        check_product(2, 0, 0)
        check_product(3, 0, 1)

    def test_ratesOfProgress(self):
        self.assertEqual(len(self.phase.netRatesOfProgress),
                         self.phase.nReactions)
        self.assertArrayNear(self.phase.fwdRatesOfProgress - self.phase.revRatesOfProgress,
                             self.phase.netRatesOfProgress)

    def test_rateConstants(self):
        self.assertEqual(len(self.phase.fwdRateConstants), self.phase.nReactions)
        self.assertArrayNear(self.phase.fwdRateConstants / self.phase.revRateConstants,
                             self.phase.equilibriumConstants)

    def test_species_rates(self):
        nu_p = self.phase.productStoichCoeffs()
        nu_r = self.phase.reactantStoichCoeffs()
        creation = (np.dot(nu_p, self.phase.fwdRatesOfProgress) +
                    np.dot(nu_r, self.phase.revRatesOfProgress))
        destruction = (np.dot(nu_r, self.phase.fwdRatesOfProgress) +
                       np.dot(nu_p, self.phase.revRatesOfProgress))

        self.assertArrayNear(self.phase.creationRates, creation)
        self.assertArrayNear(self.phase.destructionRates, destruction)
        self.assertArrayNear(self.phase.netProductionRates,
                             creation - destruction)

    def test_reaction_deltas(self):
        self.assertArrayNear(self.phase.deltaEnthalpy -
                             self.phase.deltaEntropy * self.phase.T,
                             self.phase.deltaGibbs)
        self.assertArrayNear(self.phase.deltaStandardEnthalpy -
                             self.phase.deltaStandardEntropy * self.phase.T,
                             self.phase.deltaStandardGibbs)
