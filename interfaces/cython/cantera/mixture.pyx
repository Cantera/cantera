from mixture cimport *

cdef class Mixture:
    def __cinit__(self, phases):
        self.mix = new CxxMultiPhase()
        self._phases = []

        cdef _SolutionBase phase
        for phase,moles in phases:
            self.mix.addPhase(phase.thermo, moles)
            self._phases.append(phase)

        self.mix.init()
        if self._phases:
            self.pressure = self._phases[0].pressure
            self.temperature = self._phases[0].temperature

    def __dealloc__(self):
        del self.mix

    def phase(self, n):
        return self._phases[n]

    property nSpecies:
        def __get__(self):
            return self.mix.nSpecies()

    property temperature:
        def __get__(self):
            return self.mix.temperature()
        def __set__(self, T):
            self.mix.setTemperature(T)

    property pressure:
        def __get__(self):
            return self.mix.pressure()
        def __set__(self, P):
            self.mix.setPressure(P)
