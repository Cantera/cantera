######################################################
print """

  Tutorial 3: Getting Help

"""
######################################################

# Python has a built-in help facility. To get help on any class or
# function (in Cantera or not), import it and the call 'help' with the
# class or function as the argument

from Cantera import Solution
help(Solution)

from Cantera import Reactor
help(Reactor)

# You can also use the Python module browser to view this same
# information in a web browser. Under Windows, goto
# Programs/Python2.x/Module Docs on the Start menu. On unix or Mac
# OSX, type 'pydoc -g' at a shell prompt, A small pop-up window will
# appear. Click on 'open browser', then navigate to the Cantera module
