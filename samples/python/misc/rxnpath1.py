"""
Viewing a reaction path diagram.

This script uses Graphviz to generate an image. You must have Graphviz installed
and the program 'dot' must be on your path for this example to work.
Graphviz can be obtained from http://www.graphviz.org/ or (possibly) installed
using your operating system's package manager.
"""

import os
import sys

from Cantera import *
from Cantera import rxnpath

#-----------------------------------------------------------------------
# these lines can be replaced by any commands that generate
# an object of a class derived from class Kinetics (such as IdealGasMix)
# in some state.
gas = GRI30()
gas.setState_TPX(2500.0, OneAtm, 'CH4:0.4, O2:1, N2:3.76')
gas.equilibrate('TP')
gas.setTemperature(500.0)

#------------------------------------------------------------------------

d = rxnpath.PathDiagram(title = 'reaction path diagram following N',
                        bold_color = 'green')

element = 'N'
dot_file = 'rxnpath2.dot'
img_file = 'rxnpath2.png'
img_path = os.path.join(os.getcwd(), img_file)

rxnpath.write(gas, element, dot_file, d)
print "Wrote graphviz input file to '%s'." % os.path.join(os.getcwd(), dot_file)

os.system('dot %s -Tpng -o%s -Gdpi=200' % (dot_file, img_file))
print "Wrote graphviz output file to '%s'." % img_path

if len(sys.argv) > 1 and sys.argv[1] == "-view":
    rxnpath.view('file:///' + img_path)
