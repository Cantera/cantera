"""
Viewing a reaction path diagram.

This script uses Graphviz to generate an image. You must have Graphviz installed
and the program 'dot' must be on your path for this example to work.
Graphviz can be obtained from https://www.graphviz.org/ or (possibly) installed
using your operating system's package manager.

Requires: cantera >= 2.5.0
"""

import os
import sys

import cantera as ct

# these lines can be replaced by any commands that generate
# an object of a class derived from class Kinetics in some state.
gas = ct.Solution('gri30.yaml')
gas.TPX = 1300.0, ct.one_atm, 'CH4:0.4, O2:1, N2:3.76'
r = ct.IdealGasReactor(gas)
net = ct.ReactorNet([r])
T = r.T
while T < 1900:
    net.step()
    T = r.T

element = 'N'

diagram = ct.ReactionPathDiagram(gas, element)
diagram.title = 'Reaction path diagram following {0}'.format(element)
diagram.label_threshold = 0.01

dot_file = 'rxnpath.dot'
img_file = 'rxnpath.png'
img_path = os.path.join(os.getcwd(), img_file)

diagram.write_dot(dot_file)
print(diagram.get_data())

print("Wrote graphviz input file to '{0}'.".format(os.path.join(os.getcwd(), dot_file)))

os.system('dot {0} -Tpng -o{1} -Gdpi=200'.format(dot_file, img_file))
print("Wrote graphviz output file to '{0}'.".format(img_path))

if "-view" in sys.argv:
    import webbrowser
    webbrowser.open('file:///' + img_path)
