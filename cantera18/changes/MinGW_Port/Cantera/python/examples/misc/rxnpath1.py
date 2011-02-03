"""

Viewing a reaction path diagram in a web browser.

This script uses the public webdot server at webdot.graphviz.org to
perform the rendering of the graph. You don't need to get the 'dot'
program to use this script, but you do need to be able to put your
output file where it can be served by your web server. You can either
run a local server, or mount the remote directory where your web files
are located. You could also modify this script to write the file
locally, and then upload it to your web server.

"""

from os.path import join
from Cantera import *
from Cantera import rxnpath

import sys
opts = sys.argv

#------------ site-specific configuration -----------------------------

# to view the diagram in a browser, set 'output_dir' to a directory that
# is within the document tree of your web server, and set output_url to
# the URL that accesses this directory.

output_dir = 'd:/www/docs'
output_urldir = 'http://your.http.server/'

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
output_file = 'rxnpath2.dot'
url = output_urldir+'/'+output_file

rxnpath.write(gas, element, join(output_dir, output_file), d)

if len(opts) > 1 and opts[1] == "-view":
    # graphics format. Must be one of png, svg, gif, or jpg
    fmt = 'svg'
    rxnpath.view(url, fmt)


