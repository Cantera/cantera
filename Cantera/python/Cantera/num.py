
nummodule = None
try:
    import numarray
    nummodule = numarray
except:
    try:
        import Numeric
        nummodule = Numeric
    except:
        print """
        
        ERROR: numarray / numeric not found!
        
        Cantera uses a set of numerical extensions to Python, but these do
        not appear to be present on your system. To install the required
        package, go to http://sourceforge.net/projects/numpy, and install
        either the numarray or Numeric package for your system. If you are
        using a Windows system, use the binary installer to install the
        selected package for you automatically.
        
        """
        raise "could not import numarray or Numeric"

zeros = nummodule.zeros
array = nummodule.array
asarray = nummodule.asarray
transpose = nummodule.transpose
ravel = nummodule.ravel

