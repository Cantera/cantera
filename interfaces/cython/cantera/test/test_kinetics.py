import numpy as np
import re
import itertools
from os.path import join as pjoin
import os

import cantera as ct
from . import utilities


class TestKinetics(utilities.CanteraTest):
    def setUp(self):
        self.phase = ct.Solution('h2o2.xml')
        self.phase.X = [0.1, 1e-4, 1e-5, 0.2, 2e-4, 0.3, 1e-6, 5e-5, 0.4]
        self.phase.TP = 800, 2*ct.one_atm

    def test_counts(self):
        self.assertEqual(self.phase.n_reactions, 28)
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
        self.assertNear(self.phase.reaction_type(20), 4) # falloff

        with self.assertRaisesRegex(ValueError, 'out of range'):
            self.phase.reaction_type(33)
        with self.assertRaisesRegex(ValueError, 'out of range'):
            self.phase.reaction_type(-2)

    def test_reaction_equations(self):
        self.assertEqual(self.phase.n_reactions,
                         len(self.phase.reaction_equations()))
        r,p = [x.split() for x in self.phase.reaction_equation(17).split('<=>')]
        self.assertIn('H', r)
        self.assertIn('H2O2', r)
        self.assertIn('HO2', p)
        self.assertIn('H2', p)

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
        check_reactant('H', 17, 1)
        check_reactant('H2O2', 17, 1)
        check_reactant('HO2', 17, 0)
        check_reactant('H2', 17, 0)

        check_product('H', 17, 0)
        check_product('H2O2', 17, 0)
        check_product('HO2', 17, 1)
        check_product('H2', 17, 1)

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

    def test_heat_release(self):
        hrr = - self.phase.partial_molar_enthalpies.dot(self.phase.net_production_rates)
        self.assertNear(hrr, self.phase.heat_release_rate)
        self.assertNear(hrr, sum(self.phase.heat_production_rates))

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


class KineticsFromReactions(utilities.CanteraTest):
    """
    Test for Kinetics objects which are constructed directly from Reaction
    objects instead of from input files.
    """
    def test_idealgas(self):
        gas1 = ct.Solution('h2o2.xml')

        S = ct.Species.listFromFile('h2o2.xml')
        R = ct.Reaction.listFromFile('h2o2.xml')
        gas2 = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                           species=S, reactions=R)

        self.assertEqual(gas1.n_reactions, gas2.n_reactions)
        gas1.TPY = 800, 2*ct.one_atm, 'H2:0.3, O2:0.7, OH:2e-4, O:1e-3, H:5e-5'
        gas2.TPY = gas1.TPY

        self.assertTrue((gas1.reactant_stoich_coeffs() ==
                         gas2.reactant_stoich_coeffs()).all())
        self.assertTrue((gas1.product_stoich_coeffs() ==
                         gas2.product_stoich_coeffs()).all())

        self.assertArrayNear(gas1.delta_gibbs,
                             gas2.delta_gibbs)
        self.assertArrayNear(gas1.reverse_rate_constants,
                             gas2.reverse_rate_constants)
        self.assertArrayNear(gas1.net_production_rates,
                             gas2.net_production_rates)

    def test_surface(self):
        gas_species = ct.Species.listFromFile('gri30.xml')
        surf_species = ct.Species.listFromFile('ptcombust.xml')
        reactions = ct.Reaction.listFromFile('ptcombust.xml')

        gas = ct.Solution('ptcombust.xml', 'gas')
        surf1 = ct.Interface('ptcombust.xml', 'Pt_surf', [gas])

        surf2 = ct.Interface(thermo='Surface', kinetics='interface',
                             species=surf_species, reactions=reactions,
                             adjacent=[gas])
        surf1.site_density = surf2.site_density = 5e-9
        gas.TP = surf2.TP = surf1.TP = 900, 2*ct.one_atm
        surf2.concentrations = surf1.concentrations

        self.assertEqual(surf1.n_reactions, surf2.n_reactions)

        for k,i in itertools.product(['PT(S)','H2','OH','OH(S)'],
                                     range(surf1.n_species)):
            self.assertEqual(surf1.reactant_stoich_coeff(k,i),
                             surf2.reactant_stoich_coeff(k,i))
            self.assertEqual(surf1.product_stoich_coeff(k,i),
                             surf2.product_stoich_coeff(k,i))

        for i in range(surf1.n_reactions):
            r1 = surf1.reaction(i)
            r2 = surf2.reaction(i)
            self.assertEqual(r1.reactants, r2.reactants)
            self.assertEqual(r1.products, r2.products)
            self.assertEqual(r1.rate.pre_exponential_factor,
                             r2.rate.pre_exponential_factor)
            self.assertEqual(r1.rate.temperature_exponent,
                             r2.rate.temperature_exponent)
            self.assertEqual(r1.rate.activation_energy,
                             r2.rate.activation_energy)

        self.assertArrayNear(surf1.delta_enthalpy,
                             surf2.delta_enthalpy)
        self.assertArrayNear(surf1.forward_rate_constants,
                             surf2.forward_rate_constants)
        self.assertArrayNear(surf1.reverse_rate_constants,
                             surf2.reverse_rate_constants)

        rop1 = surf1.net_production_rates
        rop2 = surf2.net_production_rates
        for k in gas.species_names + surf1.species_names:
            k1 = surf1.kinetics_species_index(k)
            k2 = surf2.kinetics_species_index(k)
            self.assertNear(rop1[k1], rop2[k2])

    def test_add_reaction(self):
        gas1 = ct.Solution('h2o2.xml')

        S = ct.Species.listFromFile('h2o2.xml')
        R = ct.Reaction.listFromFile('h2o2.xml')
        gas2 = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                           species=S, reactions=R[:5])

        gas1.TPY = 800, 2*ct.one_atm, 'H2:0.3, O2:0.7, OH:2e-4, O:1e-3, H:5e-5'
        gas2.TPY = gas1.TPY

        for r in R[5:]:
            gas2.add_reaction(r)

        self.assertEqual(gas1.n_reactions, gas2.n_reactions)

        self.assertTrue((gas1.reactant_stoich_coeffs() ==
                         gas2.reactant_stoich_coeffs()).all())
        self.assertTrue((gas1.product_stoich_coeffs() ==
                         gas2.product_stoich_coeffs()).all())

        self.assertArrayNear(gas1.delta_gibbs,
                             gas2.delta_gibbs)
        self.assertArrayNear(gas1.reverse_rate_constants,
                             gas2.reverse_rate_constants)
        self.assertArrayNear(gas1.net_production_rates,
                             gas2.net_production_rates)


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

    def test_modify_thermo(self):
        # Make sure that thermo modifications propagate through to Kinetics

        # Set a gas state that is near enough to equilibrium that changes in the
        # reverse rate always show up in the net rate
        gas = self.setup_gas('gri30.xml')
        gas.TPX = self.T0, self.P0, self.X0
        gas.equilibrate('TP')
        gas.TP = gas.T + 20, None

        S = {sp.name: sp for sp in ct.Species.listFromFile('gri30.xml')}
        w1 = gas.net_rates_of_progress

        OH = gas.species('OH')
        OH.thermo = S['CO2'].thermo
        gas.modify_species(gas.species_index('OH'), OH)
        w2 = gas.net_rates_of_progress

        for i,R in enumerate(gas.reactions()):
            if ('OH' in R.reactants or 'OH' in R.products) and R.reversible:
                # Rate should be different if reaction involves OH
                self.assertNotAlmostEqual(w2[i] / w1[i], 1.0)
            else:
                # Rate should be the same if reaction does not involve OH
                self.assertAlmostEqual(w2[i] / w1[i], 1.0)


class TestEmptyKinetics(utilities.CanteraTest):
    def test_empty(self):
        gas = ct.Solution('air-no-reactions.xml')

        self.assertEqual(gas.n_reactions, 0)
        self.assertArrayNear(gas.creation_rates, np.zeros(gas.n_species))
        self.assertArrayNear(gas.destruction_rates, np.zeros(gas.n_species))
        self.assertArrayNear(gas.net_production_rates, np.zeros(gas.n_species))


class TestReactionPath(utilities.CanteraTest):
    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        cls.gas = ct.Solution('gri30.xml')
        cls.gas.TPX = 1300.0, ct.one_atm, 'CH4:0.4, O2:1, N2:3.76'
        r = ct.IdealGasReactor(cls.gas)
        net = ct.ReactorNet([r])
        T = r.T
        while T < 1900:
            net.step()
            T = r.T

    def check_dot(self, diagram, element):
        diagram.label_threshold = 0
        diagram.threshold = 0
        dot = diagram.get_dot()
        dot = dot.replace('\n', '')
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
            self.assertTrue(self.gas.n_atoms(spec, element) > 0)

        # return fluxes from the dot file for further tests
        return [float(re.search('label *= *"(.*?)"', line).group(1))
                for line in dot.split(';')
                if line.startswith('s') and 'arrowsize' in line]

    def get_fluxes(self, diagram):
        directional = {}
        net = {}
        for line in diagram.get_data().strip().split('\n')[1:]:
            s = line.split()
            fwd = float(s[2])
            rev = float(s[3])
            directional[s[0], s[1]] = fwd
            directional[s[1], s[0]] = -rev
            if fwd + rev > 0:
                net[s[0], s[1]] = fwd + rev
            else:
                net[s[1], s[0]] = - fwd - rev
        return directional, net

    def test_dot_net_autoscaled(self):
        for element in ['N', 'C', 'H', 'O']:
            diagram = ct.ReactionPathDiagram(self.gas, element)
            dot_fluxes = self.check_dot(diagram, element)
            self.assertEqual(max(dot_fluxes), 1.0)

    def test_dot_net_unscaled(self):
        for element in ['N', 'C', 'H', 'O']:
            diagram = ct.ReactionPathDiagram(self.gas, element)
            diagram.scale = 1.0
            dot_fluxes = sorted(self.check_dot(diagram, element))
            _, fluxes = self.get_fluxes(diagram)
            fluxes = sorted(fluxes.values())

            for i in range(1, 20):
                self.assertNear(dot_fluxes[-i], fluxes[-i], 1e-2)

    def test_dot_oneway_autoscaled(self):
        for element in ['N', 'C', 'H', 'O']:
            diagram = ct.ReactionPathDiagram(self.gas, element)
            diagram.flow_type = 'OneWayFlow'
            dot_fluxes = self.check_dot(diagram, element)
            self.assertEqual(max(dot_fluxes), 1.0)

    def test_dot_oneway_unscaled(self):
        for element in ['N', 'C', 'H', 'O']:
            diagram = ct.ReactionPathDiagram(self.gas, element)
            diagram.scale = 1.0
            diagram.flow_type = 'OneWayFlow'
            dot_fluxes = sorted(self.check_dot(diagram, element))
            fluxes, _ = self.get_fluxes(diagram)
            fluxes = sorted(fluxes.values())

            for i in range(1, 20):
                self.assertNear(dot_fluxes[-i], fluxes[-i], 1e-2)

    def test_fluxes(self):
        gas = ct.Solution('h2o2.cti')
        gas.TPX = 1100, 10*ct.one_atm, 'H:0.1, HO2:0.1, AR:6'
        diagram = ct.ReactionPathDiagram(gas, 'H')
        ropf = gas.forward_rates_of_progress
        ropr = gas.reverse_rates_of_progress
        fluxes, _ = self.get_fluxes(diagram)

        self.assertNear(fluxes['HO2','H'], ropr[5] + ropr[8], 1e-5)
        self.assertNear(fluxes['H', 'H2'], 2*ropf[10] + ropf[15], 1e-5)
        self.assertNear(fluxes['H', 'H2O'], ropf[14], 1e-5)
        self.assertNear(fluxes['H', 'OH'], ropf[16], 1e-5)
        self.assertNear(fluxes['HO2','H2'], ropf[15], 1e-5)
        self.assertNear(fluxes['HO2', 'H2O'], ropf[14], 1e-5)
        self.assertNear(fluxes['HO2', 'OH'], ropf[16], 1e-5)
        self.assertNear(fluxes['HO2', 'H2O2'], 2*ropf[25] + 2*ropf[26], 1e-5)


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
        for i in range(3):
            self.assertEqual(Rr[i], 0.0)

    def test_rateConstants(self):
        # species order: [H, AR, R1A, R1B, P1]
        C = self.gas.concentrations
        Rf = self.gas.forward_rates_of_progress
        kf = self.gas.forward_rate_constants
        self.assertNear(Rf[0], kf[0] * C[2]**1.5 * C[3]**0.5)
        self.assertNear(Rf[1], kf[1] * C[0]**1.0 * C[4]**0.2)
        self.assertNear(Rf[2], kf[2] * C[2]**3.0)

    def test_ratio1(self):
        rop1 = self.gas.forward_rates_of_progress
        # Double concentration of H and R1A
        self.gas.TPX = None, None, [0.02, 0.87, 0.04, 0.03, 0.04]
        rop2 = self.gas.forward_rates_of_progress
        ratio = rop2/rop1
        self.assertNear(ratio[0], 2**1.5) # order of R1A is 1.5
        self.assertNear(ratio[1], 2**1.0) # order of H is 1.0
        self.assertNear(ratio[2], 2**3) # order of R1A is 3

    def test_ratio2(self):
        rop1 = self.gas.forward_rates_of_progress
        # Double concentration of P1 and R1B
        self.gas.TPX = None, None, [0.01, 0.83, 0.02, 0.06, 0.08]
        rop2 = self.gas.forward_rates_of_progress
        ratio = rop2/rop1
        self.assertNear(ratio[0], 2**0.5) # order of R1B is 0.5
        self.assertNear(ratio[1], 2**0.2) # order of P1 is 1.0
        self.assertNear(ratio[2], 2**0.0) # order of R1B is 0


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

            n = 0
            while True:
                n += 1
                f0 = f(x0) - C
                x0 -= f0/(f(x0 + dx) - C - f0)*dx
                if n > 1000:
                    raise Exception('No convergence in Newton solve')
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
                             [6.18736878e-01, 3.81123655e-01, 8.6303646e-05,
                              2.59274203e-06, 5.05700981e-05], 1e-7)
        self.assertArrayNear(oxide_surf_a.coverages,
                             [4.99435780e-02, 9.48927983e-01, 1.12840418e-03,
                              3.35936530e-08], 1e-7)
        self.assertArrayNear(cathode_surf.coverages,
                             [1.48180380e-07, 7.57234727e-14, 9.99999827e-01,
                              2.49235513e-08, 4.03296469e-13], 1e-7)
        self.assertArrayNear(oxide_surf_c.coverages,
                             [4.99896947e-02, 9.49804199e-01, 2.06104679e-04,
                              1.11970271e-09], 1e-7)

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

        self.compare(data, pjoin(self.test_data_dir, 'sofc-test.csv'), rtol=1e-7)


class TestDuplicateReactions(utilities.CanteraTest):
    infile = 'duplicate-reactions.cti'

    def check(self, name):
        with self.assertRaisesRegex(ct.CanteraError, 'duplicate reaction'):
            ct.Solution(self.infile, name)

    def test_forward_multiple(self):
        self.check('A')

    def test_opposite_direction1(self):
        self.check('B')

    def test_opposite_direction2(self):
        self.check('C')

    def test_opposite_direction3(self):
        self.check('D')

    def test_opposite_direction4(self):
        gas = ct.Solution(self.infile, 'E')
        self.assertEqual(gas.n_reactions, 2)

    def test_common_efficiencies(self):
        self.check('F')

    def test_disjoint_efficiencies(self):
        gas = ct.Solution(self.infile, 'G')
        self.assertEqual(gas.n_reactions, 2)

    def test_different_type(self):
        gas = ct.Solution(self.infile, 'H')
        self.assertEqual(gas.n_reactions, 2)

    def test_declared_duplicate(self):
        gas = ct.Solution(self.infile, 'I')
        self.assertEqual(gas.n_reactions, 2)

    def test_unmatched_duplicate(self):
        self.check('J')

    def test_nonreacting_species(self):
        gas = ct.Solution(self.infile, 'K')
        self.assertEqual(gas.n_reactions, 3)


class TestReaction(utilities.CanteraTest):
    @classmethod
    def setUpClass(self):
        utilities.CanteraTest.setUpClass()
        self.gas = ct.Solution('h2o2.xml')
        self.gas.X = 'H2:0.1, H2O:0.2, O2:0.7, O:1e-4, OH:1e-5, H:2e-5'
        self.gas.TP = 900, 2*ct.one_atm
        self.species = ct.Species.listFromFile('h2o2.xml')

    def test_fromCti(self):
        r = ct.Reaction.fromCti('''three_body_reaction('2 O + M <=> O2 + M',
            [1.200000e+11, -1.0, 0.0], efficiencies='AR:0.83 H2:2.4 H2O:15.4')''')

        self.assertTrue(isinstance(r, ct.ThreeBodyReaction))
        self.assertEqual(r.reactants['O'], 2)
        self.assertEqual(r.products['O2'], 1)
        self.assertEqual(r.efficiencies['H2O'], 15.4)
        self.assertEqual(r.rate.temperature_exponent, -1.0)
        self.assertIn('O', r)
        self.assertIn('O2', r)
        self.assertNotIn('H2O', r)

    def test_fromXml(self):
        import xml.etree.ElementTree as ET
        p = os.path.dirname(__file__)
        root = ET.parse(pjoin(p, '..', 'data', 'h2o2.xml')).getroot()
        rxn_node = root.find('.//reaction[@id="0001"]')
        r = ct.Reaction.fromXml(ET.tostring(rxn_node))

        self.assertTrue(isinstance(r, ct.ThreeBodyReaction))
        self.assertEqual(r.reactants['O'], 2)
        self.assertEqual(r.products['O2'], 1)
        self.assertEqual(r.efficiencies['H2O'], 15.4)
        self.assertEqual(r.rate.temperature_exponent, -1.0)

    def test_fromYaml(self):
        r = ct.Reaction.fromYaml(
                "{equation: 2 O + M <=> O2 + M,"
                " type: three-body,"
                " rate-constant: {A: 1.2e+11, b: -1.0, Ea: 0.0},"
                " efficiencies: {H2: 2.4, H2O: 15.4, AR: 0.83}}",
                self.gas)

        self.assertTrue(isinstance(r, ct.ThreeBodyReaction))
        self.assertEqual(r.reactants['O'], 2)
        self.assertEqual(r.products['O2'], 1)
        self.assertEqual(r.efficiencies['H2O'], 15.4)
        self.assertEqual(r.rate.temperature_exponent, -1.0)
        self.assertIn('O', r)
        self.assertIn('O2', r)
        self.assertNotIn('H2O', r)

    def test_listFromFile(self):
        R = ct.Reaction.listFromFile('h2o2.xml')
        eq1 = [r.equation for r in R]
        eq2 = [r.equation for r in self.gas.reactions()]
        self.assertEqual(eq1, eq2)

    def test_listFromCti(self):
        p = os.path.dirname(__file__)
        with open(pjoin(p, '..', 'data', 'h2o2.cti')) as f:
            R = ct.Reaction.listFromCti(f.read())
        eq1 = [r.equation for r in R]
        eq2 = [r.equation for r in self.gas.reactions()]
        self.assertEqual(eq1, eq2)

    def test_listFromXml(self):
        p = os.path.dirname(__file__)
        with open(pjoin(p, '..', 'data', 'h2o2.xml')) as f:
            R = ct.Reaction.listFromCti(f.read())
        eq1 = [r.equation for r in R]
        eq2 = [r.equation for r in self.gas.reactions()]
        self.assertEqual(eq1, eq2)

    def test_listFromYaml(self):
        yaml = """
            - equation: O + H2 <=> H + OH  # Reaction 3
              rate-constant: {A: 3.87e+04, b: 2.7, Ea: 6260.0}
            - equation: O + HO2 <=> OH + O2  # Reaction 4
              rate-constant: {A: 2.0e+13, b: 0.0, Ea: 0.0}
            - equation: O + H2O2 <=> OH + HO2  # Reaction 5
              rate-constant: {A: 9.63e+06, b: 2.0, Ea: 4000.0}
        """
        R = ct.Reaction.listFromYaml(yaml, self.gas)
        self.assertEqual(len(R), 3)
        self.assertIn('HO2', R[2].products)
        self.assertEqual(R[0].rate.temperature_exponent, 2.7)

    def test_elementary(self):
        r = ct.ElementaryReaction({'O':1, 'H2':1}, {'H':1, 'OH':1})
        r.rate = ct.Arrhenius(3.87e1, 2.7, 6260*1000*4.184)

        gas2 = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                           species=self.species, reactions=[r])
        gas2.TPX = self.gas.TPX

        self.assertNear(gas2.forward_rate_constants[0],
                        self.gas.forward_rate_constants[2])
        self.assertNear(gas2.net_rates_of_progress[0],
                        self.gas.net_rates_of_progress[2])

    def test_arrhenius_rate(self):
        R = self.gas.reaction(2)
        self.assertNear(R.rate(self.gas.T), self.gas.forward_rate_constants[2])

    def test_negative_A(self):
        species = ct.Species.listFromFile('gri30.cti')
        r = ct.ElementaryReaction('NH:1, NO:1', 'N2O:1, H:1')
        r.rate = ct.Arrhenius(-2.16e13, -0.23, 0)

        self.assertFalse(r.allow_negative_pre_exponential_factor)

        with self.assertRaisesRegex(ct.CanteraError, 'negative pre-exponential'):
            gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                              species=species, reactions=[r])

        r.allow_negative_pre_exponential_factor = True
        gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                          species=species, reactions=[r])

    def test_negative_A_falloff(self):
        species = ct.Species.listFromFile('gri30.yaml')
        r = ct.FalloffReaction('NH:1, NO:1', 'N2O:1, H:1')
        r.low_rate = ct.Arrhenius(2.16e13, -0.23, 0)
        r.high_rate = ct.Arrhenius(-8.16e12, -0.5, 0)
        self.assertFalse(r.allow_negative_pre_exponential_factor)

        with self.assertRaisesRegex(ct.CanteraError, 'pre-exponential'):
            gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                              species=species, reactions=[r])

        r.allow_negative_pre_exponential_factor = True
        # Should still fail because of mixed positive and negative A factors
        with self.assertRaisesRegex(ct.CanteraError, 'pre-exponential'):
            gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                              species=species, reactions=[r])

        r.low_rate = ct.Arrhenius(-2.16e13, -0.23, 0)
        gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                          species=species, reactions=[r])
        self.assertLess(gas.forward_rate_constants, 0)

    def test_threebody(self):
        r = ct.ThreeBodyReaction()
        r.reactants = {'O':1, 'H':1}
        r.products = {'OH':1}
        r.rate = ct.Arrhenius(5e11, -1.0, 0.0)
        r.efficiencies = {'AR':0.7, 'H2':2.0, 'H2O':6.0}

        gas2 = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                           species=self.species, reactions=[r])
        gas2.TPX = self.gas.TPX

        self.assertNear(gas2.forward_rate_constants[0],
                        self.gas.forward_rate_constants[1])
        self.assertNear(gas2.net_rates_of_progress[0],
                        self.gas.net_rates_of_progress[1])

    def test_falloff(self):
        r = ct.FalloffReaction('OH:2', 'H2O2:1')
        r.high_rate = ct.Arrhenius(7.4e10, -0.37, 0.0)
        r.low_rate = ct.Arrhenius(2.3e12, -0.9, -1700*1000*4.184)
        r.falloff = ct.TroeFalloff((0.7346, 94, 1756, 5182))
        r.efficiencies = {'AR':0.7, 'H2':2.0, 'H2O':6.0}
        self.assertEqual(r.falloff.type, "Troe")

        gas2 = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                           species=self.species, reactions=[r])
        gas2.TPX = self.gas.TPX

        self.assertNear(gas2.forward_rate_constants[0],
                        self.gas.forward_rate_constants[20])
        self.assertNear(gas2.net_rates_of_progress[0],
                        self.gas.net_rates_of_progress[20])

    def test_plog(self):
        gas1 = ct.Solution('pdep-test.cti')
        species = ct.Species.listFromFile('pdep-test.cti')

        r = ct.PlogReaction()
        r.reactants = {'R1A':1, 'R1B':1}
        r.products = {'P1':1, 'H':1}
        r.rates = [
            (0.01*ct.one_atm, ct.Arrhenius(1.2124e13, -0.5779, 10872.7*4184)),
            (1.0*ct.one_atm, ct.Arrhenius(4.9108e28, -4.8507, 24772.8*4184)),
            (10.0*ct.one_atm, ct.Arrhenius(1.2866e44, -9.0246, 39796.5*4184)),
            (100.0*ct.one_atm, ct.Arrhenius(5.9632e53, -11.529, 52599.6*4184))
        ]

        gas2 = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                           species=species, reactions=[r])

        gas2.X = gas1.X = 'R1A:0.3, R1B:0.6, P1:0.1'

        for P in [0.001, 0.01, 0.2, 1.0, 1.1, 9.0, 10.0, 99.0, 103.0]:
            gas1.TP = gas2.TP = 900, P * ct.one_atm
            self.assertNear(gas2.forward_rate_constants[0],
                            gas1.forward_rate_constants[0])
            self.assertNear(gas2.net_rates_of_progress[0],
                            gas1.net_rates_of_progress[0])

    def test_plog_rate(self):
        gas1 = ct.Solution('pdep-test.cti')
        gas1.TP = 800, 2*ct.one_atm
        for i in range(4):
            self.assertNear(gas1.reaction(i)(gas1.T, gas1.P),
                            gas1.forward_rate_constants[i])

    def test_chebyshev(self):
        gas1 = ct.Solution('pdep-test.cti')
        species = ct.Species.listFromFile('pdep-test.cti')

        r = ct.ChebyshevReaction()
        r.reactants = 'R5:1, H:1'
        r.products = 'P5A:1, P5B:1'
        r.set_parameters(Tmin=300.0, Tmax=2000.0, Pmin=1000, Pmax=10000000,
            coeffs=[[ 5.28830e+00, -1.13970e+00, -1.20590e-01,  1.60340e-02],
                    [ 1.97640e+00,  1.00370e+00,  7.28650e-03, -3.04320e-02],
                    [ 3.17700e-01,  2.68890e-01,  9.48060e-02, -7.63850e-03],
                    [-3.12850e-02, -3.94120e-02,  4.43750e-02,  1.44580e-02]])

        gas2 = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                           species=species, reactions=[r])

        gas2.X = gas1.X = 'R5:0.3, P5A:0.6, H:0.1'

        for T,P in itertools.product([300, 500, 1500], [1e4, 4e5, 3e6]):
            gas1.TP = gas2.TP = T, P
            self.assertNear(gas2.forward_rate_constants[0],
                            gas1.forward_rate_constants[4])
            self.assertNear(gas2.net_rates_of_progress[0],
                            gas1.net_rates_of_progress[4])

    def test_chebyshev_single_P(self):
        species = ct.Species.listFromFile('pdep-test.cti')
        r = ct.ChebyshevReaction()
        r.reactants = 'R5:1, H:1'
        r.products = 'P5A:1, P5B:1'
        r.set_parameters(Tmin=300.0, Tmax=2000.0, Pmin=1000, Pmax=10000000,
            coeffs=[[ 5.28830e+00],
                    [ 1.97640e+00],
                    [ 3.17700e-01],
                    [-3.12850e-02]])

        gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                          species=species, reactions=[r])

        # rate constant should be pressure independent
        for T in [300, 500, 1500]:
            gas.TP = T, 1e4
            k1 = gas.forward_rate_constants[0]
            gas.TP = T, 1e6
            k2 = gas.forward_rate_constants[0]
            self.assertNear(k1, k2)

    def test_chebyshev_single_T(self):
        species = ct.Species.listFromFile('pdep-test.cti')
        r = ct.ChebyshevReaction()
        r.reactants = 'R5:1, H:1'
        r.products = 'P5A:1, P5B:1'
        r.set_parameters(Tmin=300.0, Tmax=2000.0, Pmin=1000, Pmax=10000000,
            coeffs=[[ 5.28830e+00, -1.13970e+00, -1.20590e-01,  1.60340e-02]])

        gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                          species=species, reactions=[r])

        # rate constant should be temperature independent
        for P in [1e4, 2e5, 8e6]:
            gas.TP = 400, P
            k1 = gas.forward_rate_constants[0]
            gas.TP = 1700, P
            k2 = gas.forward_rate_constants[0]
            self.assertNear(k1, k2)

    def test_chebyshev_rate(self):
        gas1 = ct.Solution('pdep-test.cti')
        gas1.TP = 800, 2*ct.one_atm
        for i in range(4,6):
            self.assertNear(gas1.reaction(i)(gas1.T, gas1.P),
                            gas1.forward_rate_constants[i])

    def test_chebyshev_bad_shape_cti(self):
        with self.assertRaisesRegex(ct.CanteraError, "same number"):
            r = ct.Reaction.fromCti('''
                chebyshev_reaction('R5 + H (+ M) <=> P5A + P5B (+M)',
                   Tmin=300.0, Tmax=2000.0,
                   Pmin=(0.00986, 'atm'), Pmax=(98.6, 'atm'),
                   coeffs=[[ 8.28830e+00, -1.13970e+00, -1.20590e-01],
                           [ 1.97640e+00,  1.00370e+00,  7.28650e-03, -3.04320e-02],
                           [ 3.17700e-01,  2.68890e-01,  9.48060e-02, -7.63850e-03],
                           [-3.12850e-02, -3.94120e-02,  4.43750e-02,  1.44580e-02]])''')

    def test_chebyshev_bad_shape_yaml(self):
        species = ct.Species.listFromFile('pdep-test.xml')
        gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                          species=species, reactions=[])

        with self.assertRaisesRegex(ct.CanteraError, "Inconsistent"):
            r = ct.Reaction.fromYaml('''
                equation: R5 + H (+ M) <=> P5A + P5B (+M)
                type: Chebyshev
                temperature-range: [300.0, 2000.0]
                pressure-range: [9.86e-03 atm, 98.6 atm]
                data:
                - [8.2883, -1.1397, -0.12059, 0.016034]
                - [1.9764, 1.0037, 7.2865e-03]
                - [0.3177, 0.26889, 0.094806, -7.6385e-03]
                - [-0.031285, -0.039412, 0.044375, 0.014458]''', gas)

    def test_interface(self):
        surf_species = ct.Species.listFromFile('ptcombust.xml')
        gas = ct.Solution('ptcombust.xml', 'gas')
        surf1 = ct.Interface('ptcombust.xml', 'Pt_surf', [gas])
        r1 = ct.InterfaceReaction()
        r1.reactants = 'H(S):2'
        r1.products = 'H2:1, PT(S):2'
        r1.rate = ct.Arrhenius(3.7e20, 0, 67.4e6)
        r1.coverage_deps = {'H(S)': (0, 0, -6e6)}

        self.assertNear(r1.coverage_deps['H(S)'][2], -6e6)

        surf2 = ct.Interface(thermo='Surface', species=surf_species,
                             kinetics='interface', reactions=[r1], adjacent=[gas])

        surf2.site_density = surf1.site_density
        surf1.coverages = surf2.coverages = 'PT(S):0.7, H(S):0.3'
        gas.TP = surf2.TP = surf1.TP

        for T in [300, 500, 1500]:
            gas.TP = surf1.TP = surf2.TP = T, 5*ct.one_atm
            self.assertNear(surf1.forward_rate_constants[1],
                            surf2.forward_rate_constants[0])
            self.assertNear(surf1.net_rates_of_progress[1],
                            surf2.net_rates_of_progress[0])

    def test_modify_invalid(self):
        # different reaction type
        tbr = self.gas.reaction(0)
        R2 = ct.ElementaryReaction(tbr.reactants, tbr.products)
        R2.rate = tbr.rate
        with self.assertRaisesRegex(ct.CanteraError, 'types are different'):
            self.gas.modify_reaction(0, R2)

        # different reactants
        R = self.gas.reaction(7)
        with self.assertRaisesRegex(ct.CanteraError, 'Reactants are different'):
            self.gas.modify_reaction(23, R)

        # different products
        R = self.gas.reaction(14)
        with self.assertRaisesRegex(ct.CanteraError, 'Products are different'):
            self.gas.modify_reaction(15, R)

    def test_modify_elementary(self):
        gas = ct.Solution('h2o2.xml')
        gas.TPX = self.gas.TPX
        R = self.gas.reaction(2)
        A1 = R.rate.pre_exponential_factor
        b1 = R.rate.temperature_exponent
        Ta1 = R.rate.activation_energy / ct.gas_constant
        T = gas.T
        self.assertNear(A1*T**b1*np.exp(-Ta1/T), gas.forward_rate_constants[2])

        A2 = 1.5 * A1
        b2 = b1 + 0.1
        Ta2 = Ta1 * 1.2
        R.rate = ct.Arrhenius(A2, b2, Ta2 * ct.gas_constant)
        gas.modify_reaction(2, R)
        self.assertNear(A2*T**b2*np.exp(-Ta2/T), gas.forward_rate_constants[2])

    def test_modify_third_body(self):
        gas = ct.Solution('h2o2.xml')
        gas.TPX = self.gas.TPX
        R = self.gas.reaction(5)
        A1 = R.rate.pre_exponential_factor
        b1 = R.rate.temperature_exponent
        T = gas.T
        kf1 = gas.forward_rate_constants[5]

        A2 = 1.7 * A1
        b2 = b1 - 0.1
        R.rate = ct.Arrhenius(A2, b2, 0.0)
        gas.modify_reaction(5, R)
        kf2 = gas.forward_rate_constants[5]
        self.assertNear((A2*T**b2) / (A1*T**b1), kf2/kf1)

    def test_modify_falloff(self):
        gas = ct.Solution('gri30.xml')
        gas.TPX = 1100, 3 * ct.one_atm, 'CH4:1.0, O2:0.4, CO2:0.1, H2O:0.05'
        r0 = gas.reaction(11)
        self.assertEqual(r0.falloff.type, "Lindemann")
        # these two reactions happen to have the same third-body efficiencies
        r1 = gas.reaction(49)
        r2 = gas.reaction(53)
        self.assertEqual(r2.falloff.type, "Troe")
        self.assertEqual(r1.efficiencies, r2.efficiencies)
        r2.high_rate = r1.high_rate
        r2.low_rate = r1.low_rate
        r2.falloff = r1.falloff

        gas.modify_reaction(53, r2)
        kf = gas.forward_rate_constants
        self.assertNear(kf[49], kf[53])

    def test_modify_plog(self):
        gas = ct.Solution('pdep-test.cti')
        gas.TPX = 1010, 0.12 * ct.one_atm, 'R1A:0.3, R1B:0.2, H:0.1, R2:0.4'

        r0 = gas.reaction(0)
        r1 = gas.reaction(1)
        r0.rates = r1.rates
        gas.modify_reaction(0, r0)
        kf = gas.forward_rate_constants
        self.assertNear(kf[0], kf[1])

        # Removing the high-pressure rates should have no effect at low P...
        r1.rates = r1.rates[:-4]
        gas.modify_reaction(1, r1)
        self.assertNear(kf[1], gas.forward_rate_constants[1])

        # ... but should change the rate at higher pressures
        gas.TP = 1010, 12.0 * ct.one_atm
        kf = gas.forward_rates_of_progress
        self.assertNotAlmostEqual(kf[0], kf[1])

    def test_modify_chebyshev(self):
        gas = ct.Solution('pdep-test.cti')
        gas.TPX = 1010, 0.34 * ct.one_atm, 'R1A:0.3, R1B:0.2, H:0.1, R2:0.4'

        r1 = gas.reaction(4)
        r2 = gas.reaction(5)
        r1.set_parameters(r2.Tmin, r2.Tmax, r2.Pmin, r2.Pmax, r2.coeffs)

        # rates should be different before calling 'modify_reaction'
        kf = gas.forward_rate_constants
        self.assertNotAlmostEqual(kf[4], kf[5])

        gas.modify_reaction(4, r1)
        kf = gas.forward_rate_constants
        self.assertNear(kf[4], kf[5])

    def test_modify_interface(self):
        gas = ct.Solution('ptcombust.xml', 'gas')
        surf = ct.Interface('ptcombust.xml', 'Pt_surf', [gas])
        surf.coverages = 'O(S):0.1, PT(S):0.5, H(S):0.4'
        gas.TP = surf.TP

        R = surf.reaction(1)
        R.coverage_deps = {'O(S)': (0.0, 0.0, -3e6)}
        surf.modify_reaction(1, R)

        # Rate constant should now be independent of H(S) coverage, but
        # dependent on O(S) coverage
        k1 = surf.forward_rate_constants[1]
        surf.coverages = 'O(S):0.2, PT(S):0.4, H(S):0.4'
        k2 = surf.forward_rate_constants[1]
        surf.coverages = 'O(S):0.2, PT(S):0.6, H(S):0.2'
        k3 = surf.forward_rate_constants[1]
        self.assertNotAlmostEqual(k1, k2)
        self.assertNear(k2, k3)

    def test_modify_sticking(self):
        gas = ct.Solution('ptcombust.xml', 'gas')
        surf = ct.Interface('ptcombust.xml', 'Pt_surf', [gas])
        surf.coverages = 'O(S):0.1, PT(S):0.5, H(S):0.4'
        gas.TP = surf.TP

        R = surf.reaction(2)
        R.rate = ct.Arrhenius(0.25, 0, 0) # original sticking coefficient = 1.0

        k1 = surf.forward_rate_constants[2]
        surf.modify_reaction(2, R)
        k2 = surf.forward_rate_constants[2]
        self.assertNear(k1, 4*k2)

    def test_motz_wise(self):
        # Motz & Wise off for all reactions
        gas1 = ct.Solution('ptcombust.xml', 'gas')
        surf1 = ct.Interface('ptcombust.xml', 'Pt_surf', [gas1])
        surf1.coverages = 'O(S):0.1, PT(S):0.5, H(S):0.4'
        gas1.TP = surf1.TP

        # Motz & Wise correction on for some reactions
        gas2 = ct.Solution('ptcombust-motzwise.cti', 'gas')
        surf2 = ct.Interface('ptcombust-motzwise.cti', 'Pt_surf', [gas2])
        surf2.TPY = surf1.TPY

        k1 = surf1.forward_rate_constants
        k2 = surf2.forward_rate_constants

        # M&W toggled on (globally) for reactions 2 and 7
        self.assertNear(2.0 * k1[2], k2[2]) # sticking coefficient = 1.0
        self.assertNear(1.6 * k1[7], k2[7]) # sticking coefficient = 0.75

        # M&W toggled off (locally) for reaction 4
        self.assertNear(k1[4], k2[4])

        # M&W toggled on (locally) for reaction 9
        self.assertNear(2.0 * k1[9], k2[9]) # sticking coefficient = 1.0
