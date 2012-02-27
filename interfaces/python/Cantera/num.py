import _cantera
nummodule = None

try:
    if _cantera.nummod == 'numpy':
        import numpy
        nummodule = numpy
    elif _cantera.nummod == 'numarray':
        import numarray
        nummodule = numarray
    else:
        import Numeric
        nummodule = Numeric
except:
    print """

        ERROR: """+_cantera.nummod+""" not found!

        Cantera uses a set of numerical extensions to Python, but these do
        not appear to be present on your system. To install the required
        package, go to http://sourceforge.net/projects/numpy, and install
        either the """+_cantera.nummod+""" package for your system. If you are
        using a Windows system, use the binary installer to install the
        selected package for you automatically.

        """
    raise "could not import "+_cantera.nummod

zeros = nummodule.zeros
array = nummodule.array
asarray = nummodule.asarray
transpose = nummodule.transpose
ravel = nummodule.ravel
shape = nummodule.shape
ones = nummodule.ones
log10 = nummodule.log10
