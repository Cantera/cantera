ctypedef void (*transportMethod1d)(CxxTransport*, double*) except +
ctypedef void (*transportMethod2d)(CxxTransport*, size_t, double*) except +

# NOTE: These cdef functions cannot be members of Transport because they would
# cause "layout conflicts" when creating derived classes with multiple bases,
# e.g. class Solution. [Cython 0.16]
cdef np.ndarray get_transport_1d(Transport tran, transportMethod1d method):
    cdef np.ndarray[np.double_t, ndim=1] data = np.empty(tran.thermo.nSpecies())
    method(tran.transport, &data[0])
    if tran._selectedSpecies.size:
        return data[tran._selectedSpecies]
    else:
        return data

cdef np.ndarray get_transport_2d(Transport tran, transportMethod2d method):
    cdef size_t kk = tran.thermo.nSpecies()
    cdef np.ndarray[np.double_t, ndim=2] data = np.empty((kk, kk))
    method(tran.transport, kk, &data[0,0])
    return data


cdef class Transport(_SolutionBase):
    def __init__(self, *args, **kwargs):
        if self.transport == NULL:
            self.transport = newDefaultTransportMgr(self.thermo)
        super().__init__(*args, **kwargs)

    property transportModel:
        def __get__(self):
            return pystr(transportModelName(self.transport.model()))

        def __set__(self, model):
            cdef CxxTransport* old = self.transport
            self.transport = newTransportMgr(stringify(model), self.thermo)
            del old # only if the new transport manager was successfully created

    property viscosity:
        def __get__(self):
            return self.transport.viscosity()

    property thermalConductivity:
        def __get__(self):
            return self.transport.thermalConductivity()

    property mixDiffCoeffs:
        def __get__(self):
            return get_transport_1d(self, tran_getMixDiffCoeffs)

    property mixDiffCoeffsMass:
        def __get__(self):
            return get_transport_1d(self, tran_getMixDiffCoeffsMass)

    property mixDiffCoeffsMole:
        def __get__(self):
            return get_transport_1d(self, tran_getMixDiffCoeffsMole)

    property thermalDiffCoeffs:
        def __get__(self):
            return get_transport_1d(self, tran_getThermalDiffCoeffs)

    property multiDiffCoeffs:
        def __get__(self):
            return get_transport_2d(self, tran_getMultiDiffCoeffs)

    property binaryDiffCoeffs:
        def __get__(self):
            return get_transport_2d(self, tran_getBinaryDiffCoeffs)


cdef class DustyGasTransport(Transport):
    def __init__(self, *args, **kwargs):
        self.transport = newTransportMgr(stringify("DustyGas"), self.thermo)
        super().__init__(*args, **kwargs)

    property porosity:
        def __set__(self, value):
            (<CxxDustyGasTransport*>self.transport).setPorosity(value)

    property tortuosity:
        def __set__(self, value):
            (<CxxDustyGasTransport*>self.transport).setTortuosity(value)

    property meanPoreRadius:
        def __set__(self, value):
            (<CxxDustyGasTransport*>self.transport).setMeanPoreRadius(value)

    property meanParticleDiameter:
        def __set__(self, value):
            (<CxxDustyGasTransport*>self.transport).setMeanParticleDiameter(value)

    property permeability:
        def __set__(self, value):
            (<CxxDustyGasTransport*>self.transport).setPermeability(value)
