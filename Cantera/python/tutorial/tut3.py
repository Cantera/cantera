######################################################
print """

  Tutorial 3: Getting Help

"""
######################################################

# Suppose you have created a Cantera object and want to know what
# methods are available for it, and get help on using the methods. 
from Cantera import *
g = GRI30()

# The first thing you need to know is the Python class that object g
# belongs to. In Python, the class an object belongs to is stored in
# data member __class__:

print g.__class__

# To get help on this class, type
help(g.__class__)


# You can also use the Python module browser to view this same
# information in a web browser.
# Under Windows, on the Start menu
# select Programs/Python2.x/Module Docs. On unix or Mac OSX, type
# 'pydoc -g' at a shell prompt, A small pop-up window will
# appear. Click on 'open browser', then navigate to the Cantera module, and then select what you want documentation about.

# The module browser can also be started from within a Python script
# as follows:
import pydoc
pydoc.gui()

