"""
An example demonstrating how to use Species and Reaction objects to
programmatically extract a reaction submechanism. In this example, the CO/H2
oxidation reactions are extracted from the GRI 3.0 mechanism.

To test the submechanism, a premixed CO/H2 flame is simulated using the original
mechanism and the submechanism, which demonstrates that the submechanism
contains all of the important species and reactions.

Requires: cantera >= 2.5.0, matplotlib >= 2.0
"""

from timeit import default_timer
import cantera as ct
import matplotlib.pyplot as plt

all_species = ct.Species.listFromFile('gri30.yaml')
species = []

# Filter species
for S in all_species:
    comp = S.composition
    if 'C' in comp and 'H' in comp:
        # Exclude all hydrocarbon species
        continue
    if 'N' in comp and comp != {'N': 2}:
        # Exclude all nitrogen compounds except for N2
        continue
    if 'Ar' in comp:
        # Exclude Argon
        continue

    species.append(S)

species_names = {S.name for S in species}
print('Species: {0}'.format(', '.join(S.name for S in species)))

# Filter reactions, keeping only those that only involve the selected species
all_reactions = ct.Reaction.listFromFile('gri30.yaml')
reactions = []

print('\nReactions:')
for R in all_reactions:
    if not all(reactant in species_names for reactant in R.reactants):
        continue

    if not all(product in species_names for product in R.products):
        continue

    reactions.append(R)
    print(R.equation)
print('\n')

gas1 = ct.Solution('gri30.yaml')
gas2 = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                   species=species, reactions=reactions)


def solve_flame(gas):
    gas.TPX = 373, 0.05*ct.one_atm, 'H2:0.4, CO:0.6, O2:1, N2:3.76'

    # Create the flame simulation object
    sim = ct.CounterflowPremixedFlame(gas=gas, width=0.2)

    sim.reactants.mdot = 0.12  # kg/m^2/s
    sim.products.mdot = 0.06  # kg/m^2/s

    sim.set_refine_criteria(ratio=3, slope=0.1, curve=0.2)
    sim.solve(0, auto=True)
    return sim


t1 = default_timer()
sim1 = solve_flame(gas1)
t2 = default_timer()
print('Solved with GRI 3.0 in {0:.2f} seconds'.format(t2-t1))
sim2 = solve_flame(gas2)
t3 = default_timer()
print('Solved with CO/H2 submechanism in {0:.2f} seconds'.format(t3-t2))

plt.plot(sim1.grid, sim1.heat_release_rate,
         lw=2, label='GRI 3.0')

plt.plot(sim2.grid, sim2.heat_release_rate,
         'r--', lw=2, label='CO/H2 submechanism')

plt.ylabel('Heat release rate ($W/m^3$)')
plt.xlabel('position (m)')
plt.legend(loc='upper right')
plt.tight_layout()
plt.show()
