cdef enum ThermoBasis:
    mass = 0
    molar = 1

ctypedef void (*thermoMethod1d)(CxxThermoPhase*, double*) except +

cdef class ThermoPhase(_SolutionBase):
    def report(self, show_thermo=True):
        return pystr(self.thermo.report(bool(show_thermo)))

    property nElements:
        def __get__(self):
            return self.thermo.nElements()

    def elementIndex(self, name):
        return self.thermo.elementIndex(stringify(name))

    def elementName(self, m):
        return pystr(self.thermo.elementName(m))

    property nSpecies:
        def __get__(self):
            return self.thermo.nSpecies()

    def speciesName(self, k):
        return pystr(self.thermo.speciesName(k))

    def speciesIndex(self, name):
        return self.thermo.speciesIndex(stringify(name))

    property P:
        def __get__(self):
            return self.thermo.pressure()

    property T:
        def __get__(self):
            return self.thermo.temperature()

    property rho:
        def __get__(self):
            return self.thermo.density()

    cdef np.ndarray _getArray1(self, thermoMethod1d method):
        cdef np.ndarray[np.double_t, ndim=1] data = np.empty(self.nSpecies)
        method(self.thermo, &data[0])
        return data

    cdef void _setArray1(self, thermoMethod1d method, values) except *:
        if len(values) != self.nSpecies:
            raise ValueError("Array has incorrect length")

        cdef np.ndarray[np.double_t, ndim=1] data = \
            np.ascontiguousarray(values, dtype=np.double)
        method(self.thermo, &data[0])

    property molecularWeights:
        def __get__(self):
            return self._getArray1(thermo_getMolecularWeights)

    property Y:
        def __get__(self):
            return self._getArray1(thermo_getMassFractions)
        def __set__(self, Y):
            if isinstance(Y, str):
                self.thermo.setMassFractionsByName(stringify(Y))
            else:
                self._setArray1(thermo_setMassFractions, Y)

    def massFraction(self, int k):
        return self.thermo.massFraction(k)

    property X:
        def __get__(self):
            return self._getArray1(thermo_getMoleFractions)
        def __set__(self, X):
            if isinstance(X, str):
                self.thermo.setMoleFractionsByName(stringify(X))
            else:
                self._setArray1(thermo_setMoleFractions, X)

    def moleFraction(self, int k):
        return self.thermo.moleFraction(k)

    property concentrations:
        def __get__(self):
            return self._getArray1(thermo_getConcentrations)
        def __set__(self, C):
            self._setArray1(thermo_setConcentrations, C)


cdef class InterfacePhase(ThermoPhase):
    cdef CxxSurfPhase* surf
    def __cinit__(self, *args, **kwargs):
        if self.thermo.eosType() not in (thermo_type_surf, thermo_type_edge):
            raise TypeError('Underlying ThermoPhase object is of the wrong type.')
        self.surf = <CxxSurfPhase*>(self.thermo)

    property siteDensity:
        def __get__(self):
            return self.surf.siteDensity()
        def __set__(self, double value):
            self.surf.setSiteDensity(value)


cdef class PureFluid(ThermoPhase):
    property critTemperature:
        def __get__(self):
            return self.thermo.critTemperature()

    property critPressure:
        def __get__(self):
            return self.thermo.critPressure()

    property critDensity:
        def __get__(self):
            return self.thermo.critDensity()

    property vaporFraction:
        def __get__(self):
            return self.thermo.vaporFraction()

    def _setState_Psat(self, double P, double x):
        self.thermo.setState_Psat(P, x)

    def _setState_Tsat(self, double T, double x):
        self.thermo.setState_Tsat(T, x)
