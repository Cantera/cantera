cdef class Transport(_SolutionBase):
    def __init__(self, *args, **kwargs):
        if self.transport == NULL:
            self.transport = newDefaultTransportMgr(self.thermo)
        super().__init__(*args, **kwargs)

    property viscosity:
        def __get__(self):
            return self.transport.viscosity()

cdef class DustyGasTransport(Transport):
    def __init__(self, *args, **kwargs):
        self.transport = newTransportMgr(stringify("DustyGas"), self.thermo)
        super().__init__(*args, **kwargs)

    def setPorosity(self, value):
        (<CxxDustyGasTransport*>self.transport).setPorosity(value)
