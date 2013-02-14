import unittest
import numpy as np

import cantera as ct
from . import utilities

class TestKinetics(utilities.CanteraTest):
    def setUp(self):
        self.phase = ct.Solution('h2o2.xml')
        self.phase.X = [0.1, 1e-4, 1e-5, 0.2, 2e-4, 0.3, 1e-6, 5e-5, 0.4]
        self.phase.TP = 800, 2*ct.one_atm

    def test_counts(self):
        self.assertEqual(self.phase.n_reactions, 27)
        self.assertEqual(self.phase.n_total_species, 9)
        self.assertEqual(self.phase.n_phases, 1)
        self.assertEqual(self.phase.reaction_phase_index, 0)

    def test_is_reversible(self):
        for i in range(self.phase.n_reactions):
            self.assertTrue(self.phase.is_reversible(i))

    def test_multiplier(self):
        fwd_rates0 = self.phase.forward_rates_of_progress
        rev_rates0 = self.phase.reverse_rates_of_progress

        self.phase.set_multiplier(2.0, 0)
        self.phase.set_multiplier(0.1, 6)

        fwd_rates1 = self.phase.forward_rates_of_progress
        rev_rates1 = self.phase.reverse_rates_of_progress

        self.assertNear(2 * fwd_rates0[0], fwd_rates1[0])
        self.assertNear(0.1 * fwd_rates0[6], fwd_rates1[6])
        self.assertNear(2 * rev_rates0[0], rev_rates1[0])
        self.assertNear(0.1 * rev_rates0[6], rev_rates1[6])
        for i in range(self.phase.n_reactions):
            if i not in (0,6):
                self.assertNear(fwd_rates0[i], fwd_rates1[i])
                self.assertNear(rev_rates0[i], rev_rates1[i])

        self.phase.set_multiplier(0.5)
        fwd_rates2 = self.phase.forward_rates_of_progress
        rev_rates2 = self.phase.reverse_rates_of_progress
        self.assertArrayNear(0.5 * fwd_rates0, fwd_rates2)
        self.assertArrayNear(0.5 * rev_rates0, rev_rates2)

    def test_reaction_type(self):
        self.assertNear(self.phase.reaction_type(0), 2) # 3rd body
        self.assertNear(self.phase.reaction_type(2), 1) # elementary
        self.assertNear(self.phase.reaction_type(19), 4) # falloff

        self.assertRaises(ValueError, self.phase.reaction_type, 33)
        self.assertRaises(ValueError, self.phase.reaction_type, -2)

    def test_reaction_equations(self):
        self.assertEqual(self.phase.n_reactions,
                         len(self.phase.reaction_equations()))
        self.assertEqual(self.phase.reaction_equation(16),
                         'H + H2O2 <=> HO2 + H2')

    def test_stoich_coeffs(self):
        nu_r = self.phase.reactant_stoich_coeffs()
        nu_p = self.phase.product_stoich_coeffs()

        def check_reactant(k, i, value):
            self.assertEqual(self.phase.reactant_stoich_coeff(k,i), value)
            self.assertEqual(nu_r[k,i], value)

        def check_product(k, i, value):
            self.assertEqual(self.phase.product_stoich_coeff(k,i), value)
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

    def test_rates_of_progress(self):
        self.assertEqual(len(self.phase.net_rates_of_progress),
                         self.phase.n_reactions)
        self.assertArrayNear(self.phase.forward_rates_of_progress - self.phase.reverse_rates_of_progress,
                             self.phase.net_rates_of_progress)

    def test_rate_constants(self):
        self.assertEqual(len(self.phase.forward_rate_constants), self.phase.n_reactions)
        self.assertArrayNear(self.phase.forward_rate_constants / self.phase.reverse_rate_constants,
                             self.phase.equilibrium_constants)

    def test_species_rates(self):
        nu_p = self.phase.product_stoich_coeffs()
        nu_r = self.phase.reactant_stoich_coeffs()
        creation = (np.dot(nu_p, self.phase.forward_rates_of_progress) +
                    np.dot(nu_r, self.phase.reverse_rates_of_progress))
        destruction = (np.dot(nu_r, self.phase.forward_rates_of_progress) +
                       np.dot(nu_p, self.phase.reverse_rates_of_progress))

        self.assertArrayNear(self.phase.creation_rates, creation)
        self.assertArrayNear(self.phase.destruction_rates, destruction)
        self.assertArrayNear(self.phase.net_production_rates,
                             creation - destruction)

    def test_reaction_deltas(self):
        self.assertArrayNear(self.phase.delta_enthalpy -
                             self.phase.delta_entropy * self.phase.T,
                             self.phase.delta_gibbs)
        self.assertArrayNear(self.phase.delta_standard_enthalpy -
                             self.phase.delta_standard_entropy * self.phase.T,
                             self.phase.delta_standard_gibbs)


class TestEmptyKinetics(utilities.CanteraTest):
    def test_empty(self):
        gas = ct.Solution('air-no-reactions.xml')

        self.assertEqual(gas.n_reactions, 0)
        self.assertArrayNear(gas.creation_rates, np.zeros(gas.n_species))
        self.assertArrayNear(gas.destruction_rates, np.zeros(gas.n_species))
        self.assertArrayNear(gas.net_production_rates, np.zeros(gas.n_species))
