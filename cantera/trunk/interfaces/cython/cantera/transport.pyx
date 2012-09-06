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
    """
    This class is used to compute transport properties for a phase of matter.

    Not all transport properties are implemented in all transport models.
    """
    def __init__(self, *args, **kwargs):
        if self.transport == NULL:
            self.transport = newDefaultTransportMgr(self.thermo)
        super().__init__(*args, **kwargs)

    property transportModel:
        """
        Get/Set the transport model associated with this transport model.

        Setting a new transport model deletes the underlying C++ Transport
        object and replaces it with a new one implementing the specified model.
        """
        def __get__(self):
            return pystr(transportModelName(self.transport.model()))

        def __set__(self, model):
            cdef CxxTransport* old = self.transport
            self.transport = newTransportMgr(stringify(model), self.thermo)
            del old # only if the new transport manager was successfully created

    property viscosity:
        """Viscosity [Pa-s]"""
        def __get__(self):
            return self.transport.viscosity()

    property thermalConductivity:
        """Thermal conductivity. [W/m/K]"""
        def __get__(self):
            return self.transport.thermalConductivity()

    property mixDiffCoeffs:
        """
        Mixture-averaged diffusion coefficients [m^2/s] relating the
        mass-averaged diffusive fluxes (with respect to the mass averaged
        velocity) to gradients in the species mole fractions.
        """
        def __get__(self):
            return get_transport_1d(self, tran_getMixDiffCoeffs)

    property mixDiffCoeffsMass:
        """
        Mixture-averaged diffusion coefficients [m^2/s] relating the
        diffusive mass fluxes to gradients in the species mass fractions.
        """
        def __get__(self):
            return get_transport_1d(self, tran_getMixDiffCoeffsMass)

    property mixDiffCoeffsMole:
        """
        Mixture-averaged diffusion coefficients [m^2/s] relating the
        molar diffusive fluxes to gradients in the species mole fractions.
        """
        def __get__(self):
            return get_transport_1d(self, tran_getMixDiffCoeffsMole)

    property thermalDiffCoeffs:
        """
        Return a one-dimensional array of the species thermal diffusion
        coefficients [kg/m/s].
        """
        def __get__(self):
            return get_transport_1d(self, tran_getThermalDiffCoeffs)

    property multiDiffCoeffs:
        """Multicomponent diffusion coefficients [m^2/s]"""
        def __get__(self):
            return get_transport_2d(self, tran_getMultiDiffCoeffs)

    property binaryDiffCoeffs:
        """Binary diffusion coefficients [m^2/s]"""
        def __get__(self):
            return get_transport_2d(self, tran_getBinaryDiffCoeffs)


cdef class DustyGasTransport(Transport):
    """
    Implements the "dusty gas" model for transport in porous media.

    As implemented here, only species transport (`~Transport.multiDiffCoeffs`)
    is handled. The viscosity, thermal conductivity, and thermal diffusion
    coefficients are not implemented.
    """
    def __init__(self, *args, **kwargs):
        self.transport = newTransportMgr(stringify("DustyGas"), self.thermo)
        super().__init__(*args, **kwargs)

    property porosity:
        """Porosity of the porous medium [dimensionless]"""
        def __set__(self, value):
            (<CxxDustyGasTransport*>self.transport).setPorosity(value)

    property tortuosity:
        """Tortuosity of the porous medium [dimensionless]"""
        def __set__(self, value):
            (<CxxDustyGasTransport*>self.transport).setTortuosity(value)

    property meanPoreRadius:
        """Mean pore radius of the porous medium [m]"""
        def __set__(self, value):
            (<CxxDustyGasTransport*>self.transport).setMeanPoreRadius(value)

    property meanParticleDiameter:
        """Mean particle diameter of the porous medium [m]"""
        def __set__(self, value):
            (<CxxDustyGasTransport*>self.transport).setMeanParticleDiameter(value)

    property permeability:
        """Permeability of the porous medium [m^2]"""
        def __set__(self, value):
            (<CxxDustyGasTransport*>self.transport).setPermeability(value)
