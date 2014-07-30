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

        with self.assertRaises(ValueError):
            self.phase.reaction_type(33)
        with self.assertRaises(ValueError):
            self.phase.reaction_type(-2)

    def test_reaction_equations(self):
        self.assertEqual(self.phase.n_reactions,
                         len(self.phase.reaction_equations()))
        self.assertEqual(self.phase.reaction_equation(16),
                         'H + H2O2 <=> HO2 + H2')

    def test_reactants_products(self):
        for i in range(self.phase.n_reactions):
            R = self.phase.reactants(i)
            P = self.phase.products(i)
            self.assertTrue(self.phase.reaction_equation(i).startswith(R))
            self.assertTrue(self.phase.reaction_equation(i).endswith(P))
            for k in range(self.phase.n_species):
                if self.phase.reactant_stoich_coeff(k,i) != 0:
                    self.assertIn(self.phase.species_name(k), R)
                if self.phase.product_stoich_coeff(k,i) != 0:
                    self.assertIn(self.phase.species_name(k), P)

    def test_stoich_coeffs(self):
        nu_r = self.phase.reactant_stoich_coeffs()
        nu_p = self.phase.product_stoich_coeffs()

        def check_reactant(s, i, value):
            k = self.phase.kinetics_species_index(s)
            self.assertEqual(self.phase.reactant_stoich_coeff(s,i), value)
            self.assertEqual(self.phase.reactant_stoich_coeff(k,i), value)
            self.assertEqual(nu_r[k,i], value)

        def check_product(s, i, value):
            k = self.phase.kinetics_species_index(s)
            self.assertEqual(self.phase.product_stoich_coeff(k,i), value)
            self.assertEqual(self.phase.product_stoich_coeff(s,i), value)
            self.assertEqual(nu_p[k,i], value)

        # H + H2O2 <=> HO2 + H2
        check_reactant('H', 16, 1)
        check_reactant('H2O2', 16, 1)
        check_reactant('HO2', 16, 0)
        check_reactant('H2', 16, 0)

        check_product('H', 16, 0)
        check_product('H2O2', 16, 0)
        check_product('HO2', 16, 1)
        check_product('H2', 16, 1)

        # 2 O + M <=> O2 + M
        check_reactant('O', 0, 2)
        check_reactant('O2', 0, 0)
        check_product('O', 0, 0)
        check_product('O2', 0, 1)

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
        gas.TDX = self.T0, self.rho0, self.X0
        w1 = gas.net_production_rates

        # change everything to guarantee recomputation of rates
        gas.TDX = self.T1, self.rho1, self.X1
        w2 = gas.net_production_rates

        gas.TDX = self.T0, self.rho0, self.X1
        w3 = gas.net_production_rates

        # change only composition, and make sure the rates match
        gas.TDX = self.T0, self.rho0, self.X0
        w4 = gas.net_production_rates

        self.assertArrayNear(w1, w4)

    def check_rates_temperature1(self, mech):
        gas = self.setup_gas(mech)
        gas.TDX = self.T0, self.rho0, self.X0
        w1 = gas.net_production_rates

        # change everything to guarantee recomputation of rates
        gas.TDX = self.T1, self.rho1, self.X1
        w2 = gas.net_production_rates

        gas.TDX = self.T1, self.rho0, self.X0
        w3 = gas.net_production_rates

        # change only temperature, and make sure the rates match
        gas.TDX = self.T0, self.rho0, self.X0
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
                    self.assertNotIn(spec, species)
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


class ExplicitForwardOrderTest(utilities.CanteraTest):
    def setUp(self):
        self.gas = ct.Solution('explicit-forward-order.xml')
        self.gas.TPX = 800, 101325, [0.01, 0.90, 0.02, 0.03, 0.04]

    def test_irreversibility(self):
        # Reactions are irreversible
        Rr = self.gas.reverse_rate_constants
        self.assertEqual(Rr[0], 0.0)
        self.assertEqual(Rr[1], 0.0)

    def test_rateConstants(self):
        # species order: [H, AR, R1A, R1B, P1]
        C = self.gas.concentrations
        Rf = self.gas.forward_rates_of_progress
        kf = self.gas.forward_rate_constants
        self.assertNear(Rf[0], kf[0] * C[2]**1.5 * C[3]**0.5)
        self.assertNear(Rf[1], kf[1] * C[0]**1.0 * C[4]**0.2)

    def test_ratio1(self):
        rop1 = self.gas.forward_rates_of_progress
        # Double concentration of H and R1A
        self.gas.TPX = None, None, [0.02, 0.87, 0.04, 0.03, 0.04]
        rop2 = self.gas.forward_rates_of_progress
        ratio = rop2/rop1
        self.assertNear(ratio[0], 2**1.5) # order of R1A is 1.5
        self.assertNear(ratio[1], 2**1.0) # order of H is 1.0

    def test_ratio2(self):
        rop1 = self.gas.forward_rates_of_progress
        # Double concentration of P1 and R1B
        self.gas.TPX = None, None, [0.01, 0.83, 0.02, 0.06, 0.08]
        rop2 = self.gas.forward_rates_of_progress
        ratio = rop2/rop1
        self.assertNear(ratio[0], 2**0.5) # order of R1B is 0.5
        self.assertNear(ratio[1], 2**0.2) # order of P1 is 1.0


class TestSofcKinetics(utilities.CanteraTest):
    """ Test based on sofc.py """
    def test_sofc(self):
        mech = 'sofc-test.xml'
        T = 1073.15  # T in K
        P = ct.one_atm
        TPB_length_per_area = 1.0e7  # TPB length per unit area [1/m]

        def newton_solve(f, xstart, C=0.0):
            """ Solve f(x) = C by Newton iteration. """
            x0 = xstart
            dx = 1.0e-6
            while True:
                f0 = f(x0) - C
                x0 -= f0/(f(x0 + dx) - C - f0)*dx
                if abs(f0) < 0.00001:
                    return x0

        # Anode-side phases
        gas_a, anode_bulk, oxide_a = ct.import_phases(mech,
                                                      ['gas', 'metal', 'oxide_bulk',])
        anode_surf = ct.Interface(mech, 'metal_surface', [gas_a])
        oxide_surf_a = ct.Interface(mech, 'oxide_surface', [gas_a, oxide_a])
        tpb_a = ct.Interface(mech, 'tpb', [anode_bulk, anode_surf, oxide_surf_a])

        # Cathode-side phases
        gas_c, cathode_bulk, oxide_c = ct.import_phases(mech,
                                                        ['gas', 'metal', 'oxide_bulk'])
        cathode_surf = ct.Interface(mech, 'metal_surface', [gas_c])
        oxide_surf_c = ct.Interface(mech, 'oxide_surface', [gas_c, oxide_c])
        tpb_c = ct.Interface(mech, 'tpb', [cathode_bulk, cathode_surf,
                                                 oxide_surf_c])

        def anode_curr(E):
            anode_bulk.electric_potential = E
            w = tpb_a.net_production_rates
            return ct.faraday * w[0] * TPB_length_per_area

        def cathode_curr(E):
            cathode_bulk.electric_potential = E + oxide_c.electric_potential
            w = tpb_c.net_production_rates
            return -ct.faraday * w[0] * TPB_length_per_area

        # initialization
        gas_a.TPX = T, P, 'H2:0.97, H2O:0.03'
        gas_a.equilibrate('TP')
        gas_c.TPX = T, P, 'O2:1.0, H2O:0.001'
        gas_c.equilibrate('TP')

        for p in [anode_bulk, anode_surf, oxide_surf_a, oxide_a, cathode_bulk,
                  cathode_surf, oxide_surf_c, oxide_c, tpb_a, tpb_c]:
            p.TP = T, P

        for s in [anode_surf, oxide_surf_a, cathode_surf, oxide_surf_c]:
            s.advance_coverages(50.0)

        # These values are just a regression test with no theoretical basis
        self.assertArrayNear(anode_surf.coverages,
                             [6.18736755e-01, 3.81123779e-01, 8.63037850e-05,
                              2.59274708e-06, 5.05702339e-05])
        self.assertArrayNear(oxide_surf_a.coverages,
                             [4.99435780e-02, 9.48927983e-01, 1.12840577e-03,
                              3.35936530e-08])
        self.assertArrayNear(cathode_surf.coverages,
                             [1.48180380e-07, 7.57234727e-14, 9.99999827e-01,
                              2.49235513e-08, 4.03296469e-13])
        self.assertArrayNear(oxide_surf_c.coverages,
                             [4.99896947e-02, 9.49804199e-01, 2.06104969e-04,
                              1.11970271e-09])

        Ea0 = newton_solve(anode_curr, xstart=-0.51)
        Ec0 = newton_solve(cathode_curr, xstart=0.51)

        data = []

        # vary the anode overpotential, from cathodic to anodic polarization
        for Ea in np.linspace(Ea0 - 0.25, Ea0 + 0.25, 20):
            anode_bulk.electric_potential = Ea
            curr = anode_curr(Ea)
            delta_V = curr * 5.0e-5 / 2.0
            phi_oxide_c = -delta_V
            oxide_c.electric_potential = phi_oxide_c
            oxide_surf_c.electric_potential = phi_oxide_c
            Ec = newton_solve(cathode_curr, xstart=Ec0+0.1, C=curr)
            cathode_bulk.electric_potential = phi_oxide_c + Ec
            data.append([Ea - Ea0, 0.1*curr, Ec - Ec0, delta_V,
                             cathode_bulk.electric_potential -
                             anode_bulk.electric_potential])

        self.compare(data, '../data/sofc-test.csv')
