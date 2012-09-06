cdef enum ThermoBasis:
    massBasis = 0
    molarBasis = 1

ctypedef void (*thermoMethod1d)(CxxThermoPhase*, double*) except +

cdef class ThermoPhase(_SolutionBase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.thermoBasis = massBasis

    def report(self, show_thermo=True):
        return pystr(self.thermo.report(bool(show_thermo)))

    def __call__(self):
        print(self.report())

    property name:
        def __get__(self):
            return pystr(self.thermo.name())
        def __set__(self, name):
            self.thermo.setName(stringify(name))

    property basis:
        def __get__(self):
            if self.thermoBasis == massBasis:
                return 'mass'
            else:
                return 'molar'

        def __set__(self, value):
            if value == 'mass':
                self.thermoBasis = massBasis
            elif value == 'molar':
                self.thermoBasis = molarBasis
            else:
                raise ValueError("Valid choices are 'mass' or 'molar'.")

    cdef double _massFactor(self):
        """ Conversion factor from current basis to kg """
        if self.thermoBasis == molarBasis:
            return self.thermo.meanMolecularWeight()
        else:
            return 1.0

    cdef double _moleFactor(self):
        """ Conversion factor from current basis to moles """
        if self.thermoBasis == massBasis:
            return 1.0/self.thermo.meanMolecularWeight()
        else:
            return 1.0

    ####### Composition, species, and elements ########

    property nElements:
        def __get__(self):
            return self.thermo.nElements()

    cpdef int elementIndex(self, element) except *:
        if isinstance(element, str):
            index = self.thermo.elementIndex(stringify(element))
        elif isinstance(element, (int, float)):
            index = <int>element
        else:
            raise TypeError("'element' must be a string or a number")

        if not 0 <= index < self.nElements:
            raise ValueError('No such element.')

        return index

    def elementName(self, m):
        return pystr(self.thermo.elementName(m))

    property elementNames:
        def __get__(self):
            return [self.elementName(m) for m in range(self.nElements)]

    property nSpecies:
        def __get__(self):
            return self.thermo.nSpecies()

    def speciesName(self, k):
        return pystr(self.thermo.speciesName(k))

    property speciesNames:
        def __get__(self):
            return [self.speciesName(k) for k in range(self.nSpecies)]

    cpdef int speciesIndex(self, species) except *:
        if isinstance(species, str):
            index = self.thermo.speciesIndex(stringify(species))
        elif isinstance(species, (int, float)):
            index = <int>species
        else:
            raise TypeError("'species' must be a string or a number")

        if not 0 <= index < self.nSpecies:
            raise ValueError('No such species.')

        return index

    def nAtoms(self, species, element):
        return self.thermo.nAtoms(self.speciesIndex(species),
                                  self.elementIndex(element))

    cdef np.ndarray _getArray1(self, thermoMethod1d method):
        cdef np.ndarray[np.double_t, ndim=1] data = np.empty(self.nSpecies)
        method(self.thermo, &data[0])
        return data

    cdef void _setArray1(self, thermoMethod1d method, values) except *:
        if len(values) != self.nSpecies:
            raise ValueError("Array has incorrect length")

        cdef np.ndarray[np.double_t, ndim=1] data = \
            np.ascontiguousarray(values, dtype=np.double)
        method(self.thermo, &data[0])

    property molecularWeights:
        def __get__(self):
            return self._getArray1(thermo_getMolecularWeights)

    def molecularWeight(self, species):
        return self.thermo.molecularWeight(self.speciesIndex(species))

    property meanMolecularWeight:
        def __get__(self):
            return self.thermo.meanMolecularWeight()

    property Y:
        def __get__(self):
            return self._getArray1(thermo_getMassFractions)
        def __set__(self, Y):
            if isinstance(Y, str):
                self.thermo.setMassFractionsByName(stringify(Y))
            else:
                self._setArray1(thermo_setMassFractions, Y)

    def massFraction(self, species):
        return self.thermo.massFraction(self.speciesIndex(species))

    property X:
        def __get__(self):
            return self._getArray1(thermo_getMoleFractions)
        def __set__(self, X):
            if isinstance(X, str):
                self.thermo.setMoleFractionsByName(stringify(X))
            else:
                self._setArray1(thermo_setMoleFractions, X)

    def moleFraction(self, species):
        return self.thermo.moleFraction(self.speciesIndex(species))

    property concentrations:
        def __get__(self):
            return self._getArray1(thermo_getConcentrations)
        def __set__(self, C):
            self._setArray1(thermo_setConcentrations, C)

    def concentration(self, species):
        return self.thermo.concentration(self.speciesIndex(species))

    ######## Read-only thermodynamic properties ########

    property P:
        def __get__(self):
            return self.thermo.pressure()

    property T:
        def __get__(self):
            return self.thermo.temperature()

    property density:
        def __get__(self):
            return self.thermo.density() / self._massFactor()
    property density_mass:
        def __get__(self):
            return self.thermo.density()
    property density_mole:
        def __get__(self):
            return self.thermo.molarDensity()

    property v:
        def __get__(self):
            return self._massFactor() / self.thermo.density()
    property volume_mass:
        def __get__(self):
            return 1.0 / self.thermo.density()
    property volume_mole:
        def __get__(self):
            return self.thermo.molarVolume()

    property u:
        def __get__(self):
            return self.thermo.intEnergy_mole() * self._moleFactor()
    property intEnergy_mole:
        def __get__(self):
            return self.thermo.intEnergy_mole()
    property intEnergy_mass:
        def __get__(self):
            return self.thermo.intEnergy_mass()

    property h:
        def __get__(self):
            return self.thermo.enthalpy_mole() * self._moleFactor()
    property enthalpy_mole:
        def __get__(self):
            return self.thermo.enthalpy_mole()
    property enthalpy_mass:
        def __get__(self):
            return self.thermo.enthalpy_mass()

    property s:
        def __get__(self):
            return self.thermo.entropy_mole() * self._moleFactor()
    property entropy_mole:
        def __get__(self):
            return self.thermo.entropy_mole()
    property entropy_mass:
        def __get__(self):
            return self.thermo.entropy_mass()

    property g:
        def __get__(self):
            return self.thermo.gibbs_mole() * self._moleFactor()
    property gibbs_mole:
        def __get__(self):
            return self.thermo.gibbs_mole()

    property gibbs_mass:
        def __get__(self):
            return self.thermo.gibbs_mass()

    property cv:
        def __get__(self):
            return self.thermo.cv_mole() * self._moleFactor()
    property cv_mole:
        def __get__(self):
            return self.thermo.cv_mole()
    property cv_mass:
        def __get__(self):
            return self.thermo.cv_mass()

    property cp:
        def __get__(self):
            return self.thermo.cp_mole() * self._moleFactor()
    property cp_mole:
        def __get__(self):
            return self.thermo.cp_mole()
    property cp_mass:
        def __get__(self):
            return self.thermo.cp_mass()

    ######## Methods to get/set the complete thermodynamic state ########

    property TD:
        def __get__(self):
            return self.T, self.density
        def __set__(self, values):
            self.thermo.setState_TR(values[0], values[1] * self._massFactor())

    property TDX:
        def __get__(self):
            return self.T, self.density, self.X
        def __set__(self, values):
            self.X = values[2]
            self.TD = values[:2]

    property TDY:
        def __get__(self):
            return self.T, self.density, self.Y
        def __set__(self, values):
            self.Y = values[2]
            self.TD = values[:2]

    property TP:
        def __get__(self):
            return self.T, self.P
        def __set__(self, values):
            self.thermo.setState_TP(values[0], values[1])

    property TPX:
        def __get__(self):
            return self.T, self.P, self.X
        def __set__(self, values):
            self.X = values[2]
            self.TP = values[:2]

    property TPY:
        def __get__(self):
            return self.T, self.P, self.Y
        def __set__(self, values):
            self.Y = values[2]
            self.TP = values[:2]

    property UV:
        def __get__(self):
            return self.u, self.v
        def __set__(self, values):
            self.thermo.setState_UV(values[0] / self._massFactor(),
                                    values[1] / self._massFactor())

    property UVX:
        def __get__(self):
            return self.u, self.v, self.X
        def __set__(self, values):
            self.X = values[2]
            self.UV = values[:2]

    property UVY:
        def __get__(self):
            return self.u, self.v, self.Y
        def __set__(self, values):
            self.Y = values[2]
            self.UV = values[:2]

    property HP:
        def __get__(self):
            return self.h, self.P
        def __set__(self, values):
            self.thermo.setState_HP(values[0] / self._massFactor(), values[1])

    property HPX:
        def __get__(self):
            return self.h, self.P, self.X
        def __set__(self, values):
            self.X = values[2]
            self.HP = values[:2]

    property HPY:
        def __get__(self):
            return self.h, self.P, self.Y
        def __set__(self, values):
            self.Y = values[2]
            self.HP = values[:2]

    property SP:
        def __get__(self):
            return self.s, self.P
        def __set__(self, values):
            self.thermo.setState_SP(values[0] / self._massFactor(), values[1])

    property SPX:
        def __get__(self):
            return self.s, self.P, self.X
        def __set__(self, values):
            self.X = values[2]
            self.SP = values[:2]

    property SPY:
        def __get__(self):
            return self.s, self.P, self.Y
        def __set__(self, values):
            self.Y = values[2]
            self.SP = values[:2]

    # partial molar / non-dimensional properties
    property partial_molar_enthalpies:
        def __get__(self):
            return self._getArray1(thermo_getPartialMolarEnthalpies)

    property partial_molar_entropies:
        def __get__(self):
            return self._getArray1(thermo_getPartialMolarEntropies)

    property partial_molar_int_energies:
        def __get__(self):
            return self._getArray1(thermo_getPartialMolarIntEnergies)

    property chem_potentials:
        def __get__(self):
            return self._getArray1(thermo_getChemPotentials)

    property electrochem_potentials:
        def __get__(self):
            return self._getArray1(thermo_getElectrochemPotentials)

    property partial_molar_cp:
        def __get__(self):
            return self._getArray1(thermo_getPartialMolarCp)

    property partial_molar_volumes:
        def __get__(self):
            return self._getArray1(thermo_getPartialMolarVolumes)

    property standard_enthalpies_RT:
        def __get__(self):
            return self._getArray1(thermo_getEnthalpy_RT)

    property standard_entropies_R:
        def __get__(self):
            return self._getArray1(thermo_getEntropy_R)

    property standard_intEnergies_RT:
        def __get__(self):
            return self._getArray1(thermo_getIntEnergy_RT)

    property standard_gibbs_RT:
        def __get__(self):
            return self._getArray1(thermo_getGibbs_RT)

    property standard_cp_R:
        def __get__(self):
            return self._getArray1(thermo_getCp_R)

    ######## Miscellaneous properties ########
    property isothermalCompressibility:
        def __get__(self):
            return self.thermo.isothermalCompressibility()

    property thermalExpansionCoeff:
        def __get__(self):
            return self.thermo.thermalExpansionCoeff()

    property minTemp:
        def __get__(self):
            return self.thermo.minTemp()

    property maxTemp:
        def __get__(self):
            return self.thermo.maxTemp()

    property refPressure:
        def __get__(self):
            return self.thermo.refPressure()

    property electricPotential:
        def __get__(self):
            return self.thermo.electricPotential()
        def __set__(self, double value):
            self.thermo.setElectricPotential(value)

    property elementPotentials:
        def __get__(self):
            cdef np.ndarray[np.double_t, ndim=1] data = np.zeros(self.nElements)
            self.thermo.getElementPotentials(&data[0])
            return data


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
            return self.thermo.critDensity() / self._massFactor()

    property Psat:
        def __get__(self):
            return self.thermo.satPressure(self.T)

    property Tsat:
        def __get__(self):
            return self.thermo.satTemperature(self.P)

    property X:
        def __get__(self):
            return self.thermo.vaporFraction()

    property TX:
        def __get__(self):
            return self.T, self.X
        def __set__(self, values):
            self.thermo.setState_Tsat(values[0], values[1])

    property PX:
        def __get__(self):
            return self.P, self.X
        def __set__(self, values):
            self.thermo.setState_Psat(values[0], values[1])
