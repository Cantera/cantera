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
            self.pressure = self._phases[0].P
            self.temperature = self._phases[0].T

    def __dealloc__(self):
        del self.mix

    def phase(self, n):
        return self._phases[n]

    def equilibrate(self, XY, solver=1, int estimateEquil=0, double err=1e-9,
                    int maxsteps=1000,
                    int maxiter=200, int loglevel=0, printlevel=0):
        XY = XY.upper()
        vcs_equilibrate(deref(self.mix), stringify(XY).c_str(),
                        estimateEquil, printlevel, solver, err,
                        maxsteps, maxiter, loglevel)

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
