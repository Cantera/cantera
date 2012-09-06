import numpy as np
cimport numpy as np

from cython.operator cimport dereference as deref

from utils cimport *

cdef class _SolutionBase:
    def __cinit__(self, infile, phaseid=''):
        rootNode = getCtmlTree(stringify(infile))

        # Get XML data
        cdef XML_Node* phaseNode
        if phaseid:
            phaseNode = rootNode.findID(stringify(phaseid))
        else:
            phaseNode = rootNode.findByName(stringify('phase'))
        if phaseNode is NULL:
            raise ValueError("Couldn't read phase node from XML file")

        # Thermo
        if isinstance(self, ThermoPhase):
            self.thermo = newPhase(deref(phaseNode))
        else:
            self.thermo = NULL

        # Kinetics
        cdef vector[CxxThermoPhase*] v

        if isinstance(self, Kinetics):
            v.push_back(self.thermo)
            self.kinetics = newKineticsMgr(deref(phaseNode), v)
        else:
            self.kinetics = NULL

        # Transport
        if isinstance(self, Transport):
            self.transport = newDefaultTransportMgr(self.thermo)
        else:
            self.transport = NULL

    def __dealloc__(self):
        del self.thermo
        del self.kinetics
        del self.transport


cdef class ThermoPhase(_SolutionBase):
    property nSpecies:
        def __get__(self):
            return self.thermo.nSpecies()

    property pressure:
        def __get__(self):
            return self.thermo.pressure()

    property temperature:
        def __get__(self):
            return self.thermo.temperature()

    def setMoleFractions(self, X):
        if len(X) != self.nSpecies:
            raise ValueError("Mole fraction array has incorrect length")
        cdef np.ndarray[np.double_t, ndim=1] X_c = np.ascontiguousarray(X, dtype=np.double)
        self.thermo.setMoleFractions(&X_c[0])

    property massFractions:
        def __get__(self):
            cdef np.ndarray[np.double_t, ndim=1] X_c = np.empty(self.nSpecies)
            self.thermo.getMassFractions(&X_c[0])
            return X_c


cdef class Kinetics(_SolutionBase):
    property nReactions:
        def __get__(self):
            return self.kinetics.nReactions()


cdef class Transport(_SolutionBase):
    property viscosity:
        def __get__(self):
            return self.transport.viscosity()


class Solution(ThermoPhase, Kinetics, Transport):
    def __init__(self, *args, **kwars):
        pass
