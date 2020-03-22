# This is a sample cython code that contains a sample wrapper for the function "setState_PX" in class "ThermoPhase" in namespace "Cantera".
# I have made this wrapper assuming we would pass a third parameter (datatype = string) in the function containing which units we want, which by default would be SI unit(Pascal(Pa) in this case).
# Right now I have added support for only for one Method/Function. I will extend to other Methods, Classes and Units once this sample is approved

from pint import UnitRegistry       # make sure you have pint library installed before running this code

ureg = UnitRegistry()

cdef extern from "ThermoPhase.h" namespace "Cantera":
    cdef cppclass ThermoPhase:
        ThermoPhase()
        virtual void setState_PX(doublereal p, doublereal* x);

cdef class PyThermoPhase:
    cdef ThermoPhase *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.thisptr = new ThermoPhase()
    def __dealloc__(self):
        del self.thisptr
    def setState_PX(self, doublereal p, doublereal* x, string units)
        return self.thisptr.setState_PX(p.to(ureg.Pa),x)    #converts the input from any unit to SI unit for pressure(Pascal(Pa)) to make it compatible with current Cantera library design!
