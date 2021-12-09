# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

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

cdef np.ndarray get_transport_polynomial(
        Transport tran, transportPolyMethod1i method, int index, int n_coeffs):
    cdef np.ndarray[np.double_t, ndim=1] data = np.empty(n_coeffs)
    method(tran.transport, index, &data[0])
    return data

cdef np.ndarray get_binary_transport_polynomial(
        Transport tran, transportPolyMethod2i method, int indexi, int indexj,
        int n_coeffs):
    cdef np.ndarray[np.double_t, ndim=1] data = np.empty(n_coeffs)
    method(tran.transport, indexi, indexj, &data[0])
    return data

cdef class GasTransportData:
    """
    Transport data for a single gas-phase species which can be used in
    mixture-averaged or multicomponent transport models.

    The arguments passed to the constructor are equivalent to the properties of
    the object, with values in MKS units. To set properties in non-MKS units,
    use the `set_customary_units` method.
    """
    def __cinit__(self, geometry='', diameter=-1, well_depth=-1,
                  dipole=0.0, polarizability=0.0, rotational_relaxation=0.0,
                  acentric_factor=0.0, dispersion_coefficient=0.0,
                  quadrupole_polarizability=0.0, *, init=True):
        if init:
            self._data.reset(new CxxGasTransportData(stringify(geometry),
                diameter, well_depth, dipole, polarizability,
                rotational_relaxation, acentric_factor,
                dispersion_coefficient, quadrupole_polarizability))
            self.data = <CxxGasTransportData*?>self._data.get()

    cdef _assign(self, shared_ptr[CxxTransportData] other):
        self._data = other
        self.data = <CxxGasTransportData*?>self._data.get()

    def set_customary_units(self, geometry, diameter, well_depth, dipole=0.0,
                            polarizability=0.0, rotational_relaxation=0.0,
                            acentric_factor=0.0, dispersion_coefficient=0.0,
                            quadrupole_polarizability=0.0):
        """
        Set the parameters using "customary" units: diameter in Angstroms, well
        depth in Kelvin, dipole in Debye, and polarizability in Angstroms^3.
        These are the units used in in CK-style input files.
        """
        self.data.setCustomaryUnits(stringify(geometry), diameter, well_depth,
            dipole, polarizability, rotational_relaxation, acentric_factor,
            dispersion_coefficient, quadrupole_polarizability)

    property input_data:
        def __get__(self):
            return anymap_to_dict(self.data.parameters(True))

    def update_user_data(self, data):
        """
        Add the contents of the provided `dict` as additional fields when generating
        YAML phase definition files with `Solution.write_yaml` or in the data returned
        by `input_data`. Existing keys with matching names are overwritten.
        """
        self.data.input.update(dict_to_anymap(data), False)

    def clear_user_data(self):
        """
        Clear all saved input data, so that the data given by `input_data` or
        `Solution.write_yaml` will only include values generated by Cantera based on
        the current object state.
        """
        self.data.input.clear()

    property geometry:
        """
        Get/Set the string specifying the molecular geometry. One of `atom`,
        `linear`, or `nonlinear`.
        """
        def __get__(self):
            return pystr(self.data.geometry)
        def __set__(self, geometry):
            self.data.geometry = stringify(geometry)

    property diameter:
        """ Get/Set the Lennard-Jones collision diameter [m] """
        def __get__(self):
            return self.data.diameter
        def __set__(self, diameter):
            self.data.diameter = diameter

    property well_depth:
        """ Get/Set the Lennard-Jones well depth [J] """
        def __get__(self):
            return self.data.well_depth
        def __set__(self, well_depth):
            self.data.well_depth = well_depth

    property dipole:
        """ Get/Set the permanent dipole moment of the molecule [Coulomb-m]. """
        def __get__(self):
            return self.data.dipole
        def __set__(self, dipole):
            self.data.dipole = dipole

    property polarizability:
        """ Get/Set the polarizability of the molecule [m^3]. """
        def __get__(self):
            return self.data.polarizability
        def __set__(self, polarizability):
            self.data.polarizability = polarizability

    property rotational_relaxation:
        """
        Get/Set the rotational relaxation number (the number of collisions it
        takes to equilibrate the rotational degrees of freedom with the
        temperature).
        """
        def __get__(self):
            return self.data.rotational_relaxation
        def __set__(self, rotational_relaxation):
            self.data.rotational_relaxation = rotational_relaxation

    property acentric_factor:
        """ Get/Set Pitzer's acentric factor. [dimensionless] """
        def __get__(self):
            return self.data.acentric_factor
        def __set__(self, acentric_factor):
            self.data.acentric_factor = acentric_factor

    property dispersion_coefficient:
        """ Get/Set dispersion coefficient. [m^5] """
        def __get__(self):
            return self.data.dispersion_coefficient
        def __set__(self, dispersion_coefficient):
            self.data.dispersion_coefficient = dispersion_coefficient

    property quadrupole_polarizability:
        """ Get/Set quadrupole polarizability. [m^5] """
        def __get__(self):
            return self.data.quadrupole_polarizability
        def __set__(self, quadrupole_polarizability):
            self.data.quadrupole_polarizability = quadrupole_polarizability


cdef class Transport(_SolutionBase):
    """
    This class is used to compute transport properties for a phase of matter.

    Not all transport properties are implemented in all transport models.
    """
    # The signature of this function causes warnings for Sphinx documentation
    def __init__(self, *args, **kwargs):
        if self.transport == NULL:
            # @todo ... after removal of CTI/XML, this should be handled by base.pyx
            if 'transport_model' not in kwargs:
                self.base.setTransport(newTransport(self.thermo, stringify("default")))
            else:
                model = kwargs['transport_model']
                if not model:
                    model = 'None'
                self.base.setTransport(newTransport(self.thermo, stringify(model)))
            self.transport = self.base.transport().get()

        super().__init__(*args, **kwargs)

    property transport_model:
        """
        Get/Set the transport model associated with this transport model.

        Setting a new transport model deletes the underlying C++ Transport
        object and replaces it with a new one implementing the specified model.
        """
        def __get__(self):
            return pystr(self.transport.transportType())

        def __set__(self, model):
            self.base.setTransport(newTransport(self.thermo, stringify(model)))
            self.transport = self.base.transport().get()

    property CK_mode:
        """Boolean to indicate if the chemkin interpretation is used."""
        def __get__(self):
            return self.transport.CKMode()

    property viscosity:
        """Viscosity [Pa-s]."""
        def __get__(self):
            return self.transport.viscosity()

    property species_viscosities:
        """Pure species viscosities [Pa-s]"""
        def __get__(self):
            return get_transport_1d(self, tran_getSpeciesViscosities)

    property electrical_conductivity:
        """Electrical conductivity. [S/m]."""
        def __get__(self):
            return self.transport.electricalConductivity()

    property thermal_conductivity:
        """
        Thermal conductivity. [W/m/K]
        Returns thermal conductivity of the ideal gas object using the multicomponent
        model. The value is not specific to the dusty gas model.
        """
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

    property mobilities:
        """
        Electrical mobilities of charged species [m^2/s-V]
        """
        def __get__(self):
            return get_transport_1d(self, tran_getMobilities)

    def get_viscosity_polynomial(self, i):
        """Get the polynomial fit to the logarithm of temperature for
        the viscosity of species i."""
        n_values = 4 if self.transport.CKMode() else 5
        return get_transport_polynomial(self, tran_getViscosityPolynomial, i, n_values)

    def get_thermal_conductivity_polynomial(self, i):
        """Get the polynomial fit to the logarithm of temperature for
        the thermal conductivity of species i."""
        n_values = 4 if self.transport.CKMode() else 5
        return get_transport_polynomial(self, tran_getConductivityPolynomial, i,
                                        n_values)

    def get_binary_diff_coeffs_polynomial(self, i, j):
        """Get the polynomial fit to the logarithm of temperature for
        the binary diffusion coefficient of species i and j."""
        n_values = 4 if self.transport.CKMode() else 5
        return get_binary_transport_polynomial(self, tran_getBinDiffusivityPolynomial,
                                               i, j, n_values)

    def get_collision_integral_polynomials(self, i, j):
        """Get the polynomial fit to the logarithm of temperature for
        the collision integral of species i and j."""
        n_values = 7 if self.transport.CKMode() else 9
        cdef np.ndarray[np.double_t, ndim=1] adata = np.empty(n_values)
        cdef np.ndarray[np.double_t, ndim=1] bdata = np.empty(n_values)
        cdef np.ndarray[np.double_t, ndim=1] cdata = np.empty(n_values)
        self.transport.getCollisionIntegralPolynomial(i, j, &adata[0], &bdata[0],
                                                      &cdata[0])
        return adata, bdata, cdata

    def set_viscosity_polynomial(self, i, values):
        """Set the polynomial fit to the logarithm of temperature for
        the viscosity of species i."""
        n_values = 4 if self.transport.CKMode() else 5
        if len(values) != n_values:
            raise ValueError(
                f"""Array has incorrect length: expected {n_values} but
                received {len(values)}.""")
        cdef np.ndarray[np.double_t, ndim=1] data = np.ascontiguousarray(values,
                                                                        dtype=np.double)
        tran_setViscosityPolynomial(self.transport, i, &data[0])

    def set_thermal_conductivity_polynomial(self, i, values):
        """Set the polynomial fit to the logarithm of temperature for
        the thermal conductivity of species i."""
        n_values = 4 if self.transport.CKMode() else 5
        if len(values) != n_values:
            raise ValueError(
                f"""Array has incorrect length: expected {n_values} but
                received {len(values)}.""")
        cdef np.ndarray[np.double_t, ndim=1] data = np.ascontiguousarray(values,
                                                                        dtype=np.double)
        tran_setConductivityPolynomial(self.transport, i, &data[0])

    def set_binary_diff_coeffs_polynomial(self, i, j, values):
        """Set the polynomial fit to the logarithm of temperature for
        the binary diffusion coefficient of species i and j."""
        n_values = 4 if self.transport.CKMode() else 5
        if len(values) != n_values:
            raise ValueError(
                f"""Array has incorrect length: expected {n_values} but
                received {len(values)}.""")
        cdef np.ndarray[np.double_t, ndim=1] data = np.ascontiguousarray(values,
                                                                        dtype=np.double)
        tran_setBinDiffusivityPolynomial(self.transport, i, j, &data[0])

    def set_collision_integral_polynomial(self, i, j, avalues, bvalues, cvalues,
                                          actualT=True):
        """Get the polynomial fit to the logarithm of temperature for
        the collision integral of species i and j."""
        n_values = 7 if self.transport.CKMode() else 9
        if len(avalues) != n_values:
            raise ValueError(
                f"""Array has incorrect length: expected {n_values} but
                received {len(avalues)}.""")
        if len(bvalues) != n_values:
            raise ValueError(
                f"""Array has incorrect length: expected {n_values} but
                received {len(bvalues)}.""")
        if len(cvalues) != n_values:
            raise ValueError(
                f"""Array has incorrect length: expected {n_values} but
                received {len(cvalues)}.""")
        cdef np.ndarray[np.double_t, ndim=1] adata = np.ascontiguousarray(avalues,
                                                                        dtype=np.double)
        cdef np.ndarray[np.double_t, ndim=1] bdata = np.ascontiguousarray(bvalues,
                                                                        dtype=np.double)
        cdef np.ndarray[np.double_t, ndim=1] cdata = np.ascontiguousarray(cvalues,
                                                                        dtype=np.double)
        self.transport.setCollisionIntegralPolynomial(i, j, &adata[0], &bdata[0],
                                                      &cdata[0], actualT)

cdef class DustyGasTransport(Transport):
    """
    Implements the "dusty gas" model for transport in porous media.

    As implemented here, only species transport (`~Transport.multi_diff_coeffs`)
    is handled. The viscosity, thermal conductivity, and thermal diffusion
    coefficients are not implemented.
    """
    # The signature of this function causes warnings for Sphinx documentation
    def __init__(self, *args, **kwargs):
        self.base.setTransport(newTransport(self.thermo, stringify("DustyGas")))
        self.transport = self.base.transport().get()
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

    property thermal_conductivity:
        def __get__(self):
            return (<CxxDustyGasTransport*>self.transport).gasTransport().thermalConductivity()


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
