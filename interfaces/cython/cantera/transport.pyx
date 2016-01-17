# NOTE: These cdef functions cannot be members of Transport because they would
# cause "layout conflicts" when creating derived classes with multiple bases,
# e.g. class Solution. [Cython 0.16]
cdef np.ndarray get_transport_1d(Transport tran, transportMethod1d method):
    cdef np.ndarray[np.double_t, ndim=1] data = np.empty(tran.thermo.nSpecies())
    method(tran.transport, &data[0])
    if tran._selected_species.size:
        return data[tran._selected_species]
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

    property transport_model:
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
        """Viscosity [Pa-s]."""
        def __get__(self):
            return self.transport.viscosity()

    property electrical_conductivity:
        """Electrical conductivity. [S/m]."""
        def __get__(self):
            return self.transport.electricalConductivity()

    property thermal_conductivity:
        """Thermal conductivity. [W/m/K]."""
        def __get__(self):
            return self.transport.thermalConductivity()

    property mix_diff_coeffs:
        """
        Mixture-averaged diffusion coefficients [m^2/s] relating the
        mass-averaged diffusive fluxes (with respect to the mass averaged
        velocity) to gradients in the species mole fractions.
        """
        def __get__(self):
            return get_transport_1d(self, tran_getMixDiffCoeffs)

    property mix_diff_coeffs_mass:
        """
        Mixture-averaged diffusion coefficients [m^2/s] relating the
        diffusive mass fluxes to gradients in the species mass fractions.
        """
        def __get__(self):
            return get_transport_1d(self, tran_getMixDiffCoeffsMass)

    property mix_diff_coeffs_mole:
        """
        Mixture-averaged diffusion coefficients [m^2/s] relating the
        molar diffusive fluxes to gradients in the species mole fractions.
        """
        def __get__(self):
            return get_transport_1d(self, tran_getMixDiffCoeffsMole)

    property thermal_diff_coeffs:
        """
        Return a one-dimensional array of the species thermal diffusion
        coefficients [kg/m/s].
        """
        def __get__(self):
            return get_transport_1d(self, tran_getThermalDiffCoeffs)

    property multi_diff_coeffs:
        """Multicomponent diffusion coefficients [m^2/s]."""
        def __get__(self):
            return get_transport_2d(self, tran_getMultiDiffCoeffs)

    property binary_diff_coeffs:
        """Binary diffusion coefficients [m^2/s]."""
        def __get__(self):
            return get_transport_2d(self, tran_getBinaryDiffCoeffs)


cdef class DustyGasTransport(Transport):
    """
    Implements the "dusty gas" model for transport in porous media.

    As implemented here, only species transport (`~Transport.multi_diff_coeffs`)
    is handled. The viscosity, thermal conductivity, and thermal diffusion
    coefficients are not implemented.
    """
    def __init__(self, *args, **kwargs):
        self.transport = newTransportMgr(stringify("DustyGas"), self.thermo)
        super().__init__(*args, **kwargs)

    property porosity:
        """Porosity of the porous medium [dimensionless]."""
        def __set__(self, value):
            (<CxxDustyGasTransport*>self.transport).setPorosity(value)

    property tortuosity:
        """Tortuosity of the porous medium [dimensionless]."""
        def __set__(self, value):
            (<CxxDustyGasTransport*>self.transport).setTortuosity(value)

    property mean_pore_radius:
        """Mean pore radius of the porous medium [m]."""
        def __set__(self, value):
            (<CxxDustyGasTransport*>self.transport).setMeanPoreRadius(value)

    property mean_particle_diameter:
        """Mean particle diameter of the porous medium [m]."""
        def __set__(self, value):
            (<CxxDustyGasTransport*>self.transport).setMeanParticleDiameter(value)

    property permeability:
        """Permeability of the porous medium [m^2]."""
        def __set__(self, value):
            (<CxxDustyGasTransport*>self.transport).setPermeability(value)

    def molar_fluxes(self, T1, T2, rho1, rho2, Y1, Y2, delta):
        """
        Get the molar fluxes [kmol/m^2/s], given the thermodynamic state at
        two nearby points.

        :param T1:
            Temperature [K] at the first point
        :param T2:
            Temperature [K] at the second point
        :param rho1:
            Density [kg/m^3] at the first point
        :param rho2:
            Density [kg/m^3] at the second point
        :param Y1:
            Array of mass fractions at the first point. Length `n_species`.
        :param Y2:
            Array of mass fractions at the second point. Length `n_species`.
        :param delta:
            Distance [m] between the two points.
        """

        cdef np.ndarray[np.double_t, ndim=1] state1 = np.empty(self.n_species + 2)
        cdef np.ndarray[np.double_t, ndim=1] state2 = np.empty(self.n_species + 2)
        cdef np.ndarray[np.double_t, ndim=1] fluxes = np.empty(self.n_species)

        state1[0] = T1
        state1[1] = rho1
        state1[2:] = Y1
        state2[0] = T2
        state2[1] = rho2
        state2[2:] = Y2

        (<CxxDustyGasTransport*>self.transport).getMolarFluxes(&state1[0],
            &state2[0], delta, &fluxes[0])
        return fluxes
