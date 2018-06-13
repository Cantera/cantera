# This file is part of Cantera. See License.txt in the top-level directory or
# at http://www.cantera.org/license.txt for license and copyright information.

cdef extern from "cantera/thermo/speciesThermoTypes.h" namespace "Cantera":
    cdef int SPECIES_THERMO_CONSTANT_CP "CONSTANT_CP"
    cdef int SPECIES_THERMO_NASA2 "NASA2"
    cdef int SPECIES_THERMO_SHOMATE2 "SHOMATE2"


cdef class SpeciesThermo:
    """
    Base class for representing the reference-state thermodynamic properties of
    a pure species. These properties are a function of temperature. Derived
    classes implement a parameterization of this temperature dependence. This is
    a wrapper for the C++ class :ct:`SpeciesThermoInterpType`.

    :param T_low:
        The minimum temperature [K] at which the parameterization is valid
    :param T_high:
        The maximum temperature [K] at which the parameterization is valid
    :param P_ref:
        The reference pressure [Pa] for the parameterization
    :param coeffs:
        An array of coefficients for the parameterization. The length of this
        array and the meaning of each element depends on the specific
        parameterization.
    """
    def __cinit__(self, T_low=None, T_high=None, P_ref=None, coeffs=None, *args,
                  init=True, **kwargs):
        if not init:
            return

        if len(coeffs) != self.n_coeffs:
            raise ValueError("Coefficient array has incorrect length")
        cdef np.ndarray[np.double_t, ndim=1] data = np.ascontiguousarray(
            coeffs, dtype=np.double)
        self._spthermo.reset(CxxNewSpeciesThermo(self.derived_type, T_low,
                                                 T_high, P_ref, &data[0]))
        self.spthermo = self._spthermo.get()

    cdef _assign(self, shared_ptr[CxxSpeciesThermo] other):
        self._spthermo = other
        self.spthermo = self._spthermo.get()

    property min_temp:
        """ Minimum temperature [K] at which the parameterization is valid."""
        def __get__(self):
            return self.spthermo.minTemp()

    property max_temp:
        """ Maximum temperature [K] at which the parameterization is valid."""
        def __get__(self):
            return self.spthermo.maxTemp()

    property reference_pressure:
        """ Reference pressure [Pa] for the parameterization."""
        def __get__(self):
            return self.spthermo.refPressure()

    property coeffs:
        """
        Array of coefficients for the parameterization. The length of this
        array and the meaning of each element depends on the specific
        parameterization.
        """
        def __get__(self):
            cdef size_t index = 0
            cdef int thermo_type = 0
            cdef double T_low = 0, T_high = 0, P_ref = 0
            cdef np.ndarray[np.double_t, ndim=1] data = np.empty(self.n_coeffs)
            self.spthermo.reportParameters(index, thermo_type, T_low,
                                           T_high, P_ref, &data[0])
            return data

    def cp(self, T):
        """
        Molar heat capacity at constant pressure [J/kmol/K] at temperature *T*.
        """
        cdef double cp_r, h_rt, s_r
        self.spthermo.updatePropertiesTemp(T, &cp_r, &h_rt, &s_r)
        return cp_r * gas_constant

    def h(self, T):
        """ Molar enthalpy [J/kmol] at temperature *T* """
        cdef double cp_r, h_rt, s_r
        self.spthermo.updatePropertiesTemp(T, &cp_r, &h_rt, &s_r)
        return h_rt * gas_constant * T

    def s(self, T):
        """ Molar entropy [J/kmol/K] at temperature *T* """
        cdef double cp_r, h_rt, s_r
        self.spthermo.updatePropertiesTemp(T, &cp_r, &h_rt, &s_r)
        return s_r * gas_constant


cdef class ConstantCp(SpeciesThermo):
    r"""
    Thermodynamic properties for a species that has a constant specific heat
    capacity. This is a wrapper for the C++ class :ct:`ConstCpPoly`.

    :param coeffs:
        An array of 4 elements:

             - `coeffs[0]` = :math:`T_0` [K]
             - `coeffs[1]` = :math:`H^o(T_0, p_{ref})` [J/kmol]
             - `coeffs[2]` = :math:`S^o(T_0, p_{ref})` [J/kmol-K]
             - `coeffs[3]` = :math:`c_p^o(T_0, p_{ref})` [J/kmol-K]
    """
    derived_type = SPECIES_THERMO_CONSTANT_CP
    n_coeffs = 4


cdef class NasaPoly2(SpeciesThermo):
    """
    Thermodynamic properties for a species which is parameterized using the
    7-coefficient NASA polynomial form in two temperature ranges. This is a
    wrapper for the C++ class :ct:`NasaPoly2`.

    :param coeffs:
        An array of 15 elements, in the following order:

            - `coeffs[0]`: The mid-point temperature [K] between the two
              parameterizations
            - `coeffs[1:8]`: The 7 coefficients of the high-temperature
              parameterization
            - `coeffs[8:15]`: The 7 coefficients of the low-temperature
              parameterization

        This is the coefficient order used in the standard fixed-format NASA
        input files.
    """
    derived_type = SPECIES_THERMO_NASA2
    n_coeffs = 15


cdef class ShomatePoly2(SpeciesThermo):
    """
    Thermodynamic properties for a species which is parameterized using the
    Shomate equation in two temperature ranges. This is a wrapper for the C++
    class :ct:`ShomatePoly2`.

    :param coeffs:
        An array of 15 elements, in the following order:

            - `coeffs[0]`: The mid-point temperature [K] between the two
              parameterizations
            - `coeffs[1:8]`: The 7 coefficients of the low-temperature
              parameterization
            - `coeffs[8:15]`: The 7 coefficients of the high-temperature
              parameterization

        These coefficients should be provided in their customary units (i.e.
        such that :math:`c_p^o` is in J/gmol-K and :math:`H^o` is in kJ/gmol,
        as in the NIST Chemistry WebBook).
    """
    derived_type = SPECIES_THERMO_SHOMATE2
    n_coeffs = 15


cdef wrapSpeciesThermo(shared_ptr[CxxSpeciesThermo] spthermo):
    """
    Wrap a C++ SpeciesThermoInterpType object with a Python object of the
    correct derived type.
    """
    cdef int thermo_type = spthermo.get().reportType()

    if thermo_type == SPECIES_THERMO_NASA2:
        st = NasaPoly2(init=False)
    elif thermo_type == SPECIES_THERMO_CONSTANT_CP:
        st = ConstantCp(init=False)
    elif thermo_type == SPECIES_THERMO_SHOMATE2:
        st = ShomatePoly2(init=False)
    else:
        st = SpeciesThermo()

    st._assign(spthermo)
    return st
