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
