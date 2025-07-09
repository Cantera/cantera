"""
Viewing a reaction path diagram
===============================

This script uses Graphviz to generate an image. You must have Graphviz installed
and the program 'dot' must be on your path for this example to work.
Graphviz can be obtained from https://www.graphviz.org/ or (possibly) installed
using your operating system's package manager.

Requires: cantera >= 2.5.0, graphviz

.. tags:: Python, kinetics, reaction path analysis, pollutant formation
"""

from subprocess import run
from pathlib import Path
import graphviz

import cantera as ct

# %%
# Generate a Solution with a state where reactions are occurring
# --------------------------------------------------------------
#
# these lines can be replaced by any commands that generate
# an object of a class derived from class Kinetics in some state.
gas = ct.Solution('gri30.yaml')
gas.TPX = 1300.0, ct.one_atm, 'CH4:0.4, O2:1, N2:3.76'
r = ct.IdealGasReactor(gas, clone=False)
net = ct.ReactorNet([r])
T = r.T
while T < 1900:
    net.step()
    T = r.T

# %%
# Create a reaction path diagram following nitrogen
# -------------------------------------------------
element = 'N'

diagram = ct.ReactionPathDiagram(gas, element)
diagram.title = 'Reaction path diagram following {0}'.format(element)
diagram.label_threshold = 0.01

# %%
# Save the reaction path diagram
# ------------------------------
dot_file = 'rxnpath.dot'
img_file = 'rxnpath.png'
img_path = Path.cwd().joinpath(img_file)

diagram.write_dot(dot_file)
print(diagram.get_data())

print(f"Wrote graphviz input file to '{Path.cwd().joinpath(dot_file)}'.")

run(f"dot {dot_file} -Tpng -o{img_file} -Gdpi=200".split())
print(f"Wrote graphviz output file to '{img_path}'.")

# %%
# View the reaction path diagram
# ------------------------------
graphviz.Source(diagram.get_dot())
