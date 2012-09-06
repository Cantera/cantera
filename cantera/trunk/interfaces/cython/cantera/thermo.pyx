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
