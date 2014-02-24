import unittest
import numpy as np
import re

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


class KineticsRepeatability(utilities.CanteraTest):
    """
    Tests to make sure that lazily evaluated of terms in the rate expression
    are always updated correctly.
    """
    T0 = 1200
    T1 = 1300
    rho0 = 2.4
    rho1 = 3.1
    P0 = 1.4e5
    P1 = 3.7e6

    def setup_gas(self, mech):
        gas = ct.Solution(mech)
        self.X0 = 1 + np.sin(range(1, gas.n_species+1))
        self.X1 = 1 + np.sin(range(2, gas.n_species+2))
        return gas

    def check_rates_composition(self, mech):
        gas = self.setup_gas(mech)
        gas.TRX = self.T0, self.rho0, self.X0
        w1 = gas.net_production_rates

        # change everything to guarantee recomputation of rates
        gas.TRX = self.T1, self.rho1, self.X1
        w2 = gas.net_production_rates

        gas.TRX = self.T0, self.rho0, self.X1
        w3 = gas.net_production_rates

        # change only composition, and make sure the rates match
        gas.TRX = self.T0, self.rho0, self.X0
        w4 = gas.net_production_rates

        self.assertArrayNear(w1, w4)

    def check_rates_temperature1(self, mech):
        gas = self.setup_gas(mech)
        gas.TRX = self.T0, self.rho0, self.X0
        w1 = gas.net_production_rates

        # change everything to guarantee recomputation of rates
        gas.TRX = self.T1, self.rho1, self.X1
        w2 = gas.net_production_rates

        gas.TRX = self.T1, self.rho0, self.X0
        w3 = gas.net_production_rates

        # change only temperature, and make sure the rates match
        gas.TRX = self.T0, self.rho0, self.X0
        w4 = gas.net_production_rates

        self.assertArrayNear(w1, w4)

    def check_rates_temperature2(self, mech):
        gas = self.setup_gas(mech)
        gas.TPX = self.T0, self.P0, self.X0
        w1 = gas.net_production_rates

        # change everything to guarantee recomputation of rates
        gas.TPX = self.T1, self.P1, self.X1
        w2 = gas.net_production_rates

        gas.TPX = self.T1, self.P0, self.X0
        w3 = gas.net_production_rates

        # change only temperature, and make sure the rates match
        gas.TPX = self.T0, self.P0, self.X0
        w4 = gas.net_production_rates

        self.assertArrayNear(w1, w4)

    def check_rates_pressure(self, mech):
        gas = self.setup_gas(mech)
        gas.TPX = self.T0, self.P0, self.X0
        w1 = gas.net_production_rates

        # change everything to guarantee recomputation of rates
        gas.TPX = self.T1, self.P1, self.X1
        w2 = gas.net_production_rates

        gas.TPX = self.T0, self.P1, self.X0
        w3 = gas.net_production_rates

        # change only pressure, and make sure the rates match
        gas.TPX = self.T0, self.P0, self.X0
        w4 = gas.net_production_rates

        self.assertArrayNear(w1, w4)

    def test_gri30_composition(self):
        self.check_rates_composition('gri30.xml')

    def test_gri30_temperature(self):
        self.check_rates_temperature1('gri30.xml')
        self.check_rates_temperature2('gri30.xml')

    def test_gri30_pressure(self):
        self.check_rates_pressure('gri30.xml')

    def test_h2o2_composition(self):
        self.check_rates_composition('h2o2.xml')

    def test_h2o2_temperature(self):
        self.check_rates_temperature1('h2o2.xml')
        self.check_rates_temperature2('h2o2.xml')

    def test_h2o2_pressure(self):
        self.check_rates_pressure('h2o2.xml')

    def test_pdep_composition(self):
        self.check_rates_composition('pdep-test.xml')

    def test_pdep_temperature(self):
        self.check_rates_temperature1('pdep-test.xml')
        self.check_rates_temperature2('pdep-test.xml')

    def test_pdep_pressure(self):
        self.check_rates_pressure('pdep-test.xml')


class TestEmptyKinetics(utilities.CanteraTest):
    def test_empty(self):
        gas = ct.Solution('air-no-reactions.xml')

        self.assertEqual(gas.n_reactions, 0)
        self.assertArrayNear(gas.creation_rates, np.zeros(gas.n_species))
        self.assertArrayNear(gas.destruction_rates, np.zeros(gas.n_species))
        self.assertArrayNear(gas.net_production_rates, np.zeros(gas.n_species))


class TestReactionPath(utilities.CanteraTest):
    def test_dot_output(self):
        gas = ct.Solution('gri30.xml')
        gas.TPX = 1300.0, ct.one_atm, 'CH4:0.4, O2:1, N2:3.76'
        r = ct.IdealGasReactor(gas)
        net = ct.ReactorNet([r])
        T = r.T
        while T < 1900:
            net.step(1.0)
            T = r.T

        for element in ['N','C','H','O']:
            diagram = ct.ReactionPathDiagram(gas, element)
            diagram.label_threshold = 0.01

            dot = diagram.get_dot()
            dot = dot.replace('\n', ' ')
            nodes1 = set()
            nodes2 = set()
            species = set()
            for line in dot.split(';'):
                m = re.match(r'(.*)\[(.*)\]', line)
                if not m:
                    continue
                A, B = m.groups()
                if '->' in A:
                    # edges
                    nodes1.update(s.strip() for s in A.split('->'))
                else:
                    # nodes
                    nodes2.add(A.strip())
                    spec = re.search('label="(.*?)"', B).group(1)
                    self.assertTrue(spec not in species)
                    species.add(spec)

            # Make sure that the output was actually parsable and that we
            # found some nodes
            self.assertTrue(nodes1)
            self.assertTrue(species)

            # All nodes should be connected to some edge (this does not
            # require the graph to be connected)
            self.assertEqual(nodes1, nodes2)

            # All of the species in the graph should contain the element whose
            # flux we're looking at
            for spec in species:
                self.assertTrue(gas.n_atoms(spec, element) > 0)


class TestChemicallyActivated(utilities.CanteraTest):
    def test_rate_evaluation(self):
        gas = ct.Solution('chemically-activated-reaction.xml')
        P = [2026.5, 202650.0, 10132500.0] # pressure

        # forward rate of progress, computed using Chemkin
        Rf = [2.851022e-04, 2.775924e+00, 2.481792e+03]

        for i in range(len(P)):
            gas.TPX = 900.0, P[i], [0.01, 0.01, 0.04, 0.10, 0.84]
            self.assertNear(gas.forward_rates_of_progress[0], Rf[i], 2e-5)
