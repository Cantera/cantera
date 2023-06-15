# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *
from .solutionbase cimport *

cdef extern from "cantera/thermo/Species.h" namespace "Cantera":
    cdef cppclass CxxSpeciesThermo "Cantera::SpeciesThermoInterpType"
    cdef cppclass CxxTransportData "Cantera::TransportData"

    cdef cppclass CxxSpecies "Cantera::Species":
        CxxSpecies()
        CxxSpecies(string, Composition)
        shared_ptr[CxxSpeciesThermo] thermo
        shared_ptr[CxxTransportData] transport

        string name
        Composition composition
        double charge
        double size
        double molecularWeight() except +translate_exception
        CxxAnyMap parameters(CxxThermoPhase*) except +translate_exception
        CxxAnyMap input

    cdef shared_ptr[CxxSpecies] CxxNewSpecies "newSpecies" (CxxAnyMap&) except +translate_exception
    cdef vector[shared_ptr[CxxSpecies]] CxxGetSpecies "getSpecies" (CxxAnyValue&) except +translate_exception


cdef extern from "cantera/thermo/ThermoPhase.h" namespace "Cantera":
    ctypedef enum ThermoBasis:
        mass "Cantera::ThermoBasis::mass",
        molar "Cantera::ThermoBasis::molar"

    cdef cppclass CxxThermoPhase "Cantera::ThermoPhase":
        CxxThermoPhase()

        # miscellaneous
        string type()
        string phaseOfMatter() except +translate_exception
        void getSpeciesParameters(string, CxxAnyMap&) except +translate_exception
        CxxAnyMap& input()
        string report(cbool, double) except +translate_exception
        cbool hasPhaseTransition()
        cbool isPure()
        cbool isCompressible()
        string nativeMode()
        stdmap[string, size_t] nativeState() except +translate_exception
        vector[string] fullStates()
        vector[string] partialStates()
        double minTemp() except +translate_exception
        double maxTemp() except +translate_exception
        double refPressure() except +translate_exception
        cbool getElementPotentials(double*) except +translate_exception
        void equilibrate(string, string, double, int, int, int, int) except +translate_exception
        size_t stateSize()
        void saveState(size_t, double*) except +translate_exception
        void restoreState(size_t, double*) except +translate_exception
        CxxUnits standardConcentrationUnits() except +translate_exception

        # initialization
        void addUndefinedElements() except +translate_exception
        cbool addSpecies(shared_ptr[CxxSpecies]) except +translate_exception
        void modifySpecies(size_t, shared_ptr[CxxSpecies]) except +translate_exception
        void initThermo() except +translate_exception
        void invalidateCache() except +translate_exception

        # basic thermodynamic properties
        double temperature() except +translate_exception
        double electronTemperature() except +translate_exception
        double pressure() except +translate_exception
        double density() except +translate_exception
        double molarDensity() except +translate_exception
        double molarVolume() except +translate_exception
        double isothermalCompressibility() except +translate_exception
        double thermalExpansionCoeff() except +translate_exception
        double soundSpeed() except +translate_exception
        double electricPotential() except +translate_exception
        void setElectricPotential(double) except +translate_exception

        # element properties
        size_t nElements()
        size_t elementIndex(string) except +translate_exception
        string elementName(size_t) except +translate_exception
        double atomicWeight(size_t) except +translate_exception

        # species properties
        size_t nSpecies()
        shared_ptr[CxxSpecies] species(string) except +translate_exception
        shared_ptr[CxxSpecies] species(size_t) except +translate_exception
        size_t speciesIndex(string) except +translate_exception
        string speciesName(size_t) except +translate_exception
        double nAtoms(size_t, size_t) except +translate_exception
        void getAtoms(size_t, double*) except +translate_exception
        cbool caseSensitiveSpecies()
        void setCaseSensitiveSpecies(cbool)
        void addSpeciesAlias(string, string) except +translate_exception
        vector[string] findIsomers(Composition&) except +translate_exception
        vector[string] findIsomers(string) except +translate_exception

        double molecularWeight(size_t) except +translate_exception
        double meanMolecularWeight()

        # composition
        void setMassFractionsByName(string) except +translate_exception
        void setMassFractionsByName(Composition&) except +translate_exception
        void setMassFractions_NoNorm(double*) except +translate_exception
        Composition getMassFractionsByName(double)
        double massFraction(size_t) except +translate_exception
        double massFraction(string) except +translate_exception

        void setMoleFractionsByName(string) except +translate_exception
        void setMoleFractionsByName(Composition&) except +translate_exception
        void setMoleFractions_NoNorm(double*) except +translate_exception
        void getMoleFractions(double*) except +translate_exception
        Composition getMoleFractionsByName(double)
        double moleFraction(size_t) except +translate_exception
        double moleFraction(string) except +translate_exception

        double concentration(size_t) except +translate_exception
        double elementalMassFraction(size_t) except +translate_exception
        double elementalMoleFraction(size_t) except +translate_exception

        # state setters
        void setState_TD(double, double) except +translate_exception
        void setState_TP(double, double) except +translate_exception
        void setState_HP(double, double) except +translate_exception
        void setState_UV(double, double) except +translate_exception
        void setState_SP(double, double) except +translate_exception
        void setState_SV(double, double) except +translate_exception
        void setState_DP(double, double) except +translate_exception
        void setState_ST(double, double) except +translate_exception
        void setState_TV(double, double) except +translate_exception
        void setState_PV(double, double) except +translate_exception
        void setState_UP(double, double) except +translate_exception
        void setState_VH(double, double) except +translate_exception
        void setState_TH(double, double) except +translate_exception
        void setState_SH(double, double) except +translate_exception

        # molar thermodynamic properties:
        double enthalpy_mole() except +translate_exception
        double intEnergy_mole() except +translate_exception
        double entropy_mole() except +translate_exception
        double gibbs_mole() except +translate_exception
        double cp_mole() except +translate_exception
        double cv_mole() except +translate_exception

        # specific (mass) properties:
        double enthalpy_mass() except +translate_exception
        double intEnergy_mass() except +translate_exception
        double entropy_mass() except +translate_exception
        double gibbs_mass() except +translate_exception
        double cp_mass() except +translate_exception
        double cv_mass() except +translate_exception

        # PureFluid properties
        double critTemperature() except +translate_exception
        double critPressure() except +translate_exception
        double critDensity() except +translate_exception

        double satTemperature(double P) except +translate_exception
        double satPressure(double T) except +translate_exception
        double vaporFraction() except +translate_exception

        void setState_Tsat(double T, double x) except +translate_exception
        void setState_Psat(double P, double x) except +translate_exception
        void setState_TPQ(double T, double P, double Q) except +translate_exception

        void setMixtureFraction(double mixFrac, const double* fuelComp, const double* oxComp, ThermoBasis basis) except +translate_exception
        double mixtureFraction(const double* fuelComp, const double* oxComp, ThermoBasis basis, string element) except +translate_exception
        void setEquivalenceRatio(double phi, const double* fuelComp, const double* oxComp, ThermoBasis basis) except +translate_exception
        double equivalenceRatio(const double* fuelComp, const double* oxComp, ThermoBasis basis) except +translate_exception
        double equivalenceRatio() except +translate_exception
        double stoichAirFuelRatio(const double* fuelComp, const double* oxComp, ThermoBasis basis) except +translate_exception


cdef extern from "cantera/thermo/SurfPhase.h":
    cdef cppclass CxxSurfPhase "Cantera::SurfPhase":
        CxxSurfPhase()
        double siteDensity()
        void setSiteDensity(double) except +translate_exception
        void setCoverages(double*) except +translate_exception
        void setCoveragesByName(Composition&) except +translate_exception
        void setCoveragesNoNorm(double*) except +translate_exception
        void getCoverages(double*) except +translate_exception


cdef extern from "cantera/thermo/PlasmaPhase.h":
    cdef cppclass CxxPlasmaPhase "Cantera::PlasmaPhase" (CxxThermoPhase):
        CxxPlasmaPhase()
        void setElectronTemperature(double) except +translate_exception
        void setElectronEnergyLevels(double*, size_t) except +translate_exception
        void getElectronEnergyLevels(double*)
        void setDiscretizedElectronEnergyDist(double*, double*, size_t) except +translate_exception
        void getElectronEnergyDistribution(double*)
        void setIsotropicShapeFactor(double) except +translate_exception
        void setElectronEnergyDistributionType(const string&) except +translate_exception
        string electronEnergyDistributionType()
        void setQuadratureMethod(const string&) except +translate_exception
        string quadratureMethod()
        void enableNormalizeElectronEnergyDist(cbool)
        cbool normalizeElectronEnergyDistEnabled()
        void setMeanElectronEnergy(double) except +translate_exception
        double isotropicShapeFactor()
        double meanElectronEnergy()
        size_t nElectronEnergyLevels()
        double electronPressure()


cdef extern from "cantera/cython/thermo_utils.h":
    # ThermoPhase composition
    cdef void thermo_getMassFractions(CxxThermoPhase*, double*) except +translate_exception
    cdef void thermo_setMassFractions(CxxThermoPhase*, double*) except +translate_exception
    cdef void thermo_getMoleFractions(CxxThermoPhase*, double*) except +translate_exception
    cdef void thermo_setMoleFractions(CxxThermoPhase*, double*) except +translate_exception
    cdef void thermo_getConcentrations(CxxThermoPhase*, double*) except +translate_exception
    cdef void thermo_setConcentrations(CxxThermoPhase*, double*) except +translate_exception

    # ThermoPhase partial molar properties
    cdef void thermo_getChemPotentials(CxxThermoPhase*, double*) except +translate_exception
    cdef void thermo_getElectrochemPotentials(CxxThermoPhase*, double*) except +translate_exception
    cdef void thermo_getPartialMolarEnthalpies(CxxThermoPhase*, double*) except +translate_exception
    cdef void thermo_getPartialMolarEntropies(CxxThermoPhase*, double*) except +translate_exception
    cdef void thermo_getPartialMolarIntEnergies(CxxThermoPhase*, double*) except +translate_exception
    cdef void thermo_getPartialMolarCp(CxxThermoPhase*, double*) except +translate_exception
    cdef void thermo_getPartialMolarVolumes(CxxThermoPhase*, double*) except +translate_exception

    # ThermoPhase partial non-dimensional properties
    void thermo_getEnthalpy_RT(CxxThermoPhase*, double*) except +translate_exception
    void thermo_getEntropy_R(CxxThermoPhase*, double*) except +translate_exception
    void thermo_getIntEnergy_RT(CxxThermoPhase*, double*) except +translate_exception
    void thermo_getGibbs_RT(CxxThermoPhase*, double*) except +translate_exception
    void thermo_getCp_R(CxxThermoPhase*, double*) except +translate_exception
    void thermo_getActivities(CxxThermoPhase*, double*) except +translate_exception
    void thermo_getActivityCoefficients(CxxThermoPhase*, double*) except +translate_exception

    # other ThermoPhase methods
    cdef void thermo_getMolecularWeights(CxxThermoPhase*, double*) except +translate_exception
    cdef void thermo_getCharges(CxxThermoPhase*, double*) except +translate_exception


ctypedef void (*thermoMethod1d)(CxxThermoPhase*, double*) except +translate_exception

cdef extern from "cantera/thermo/Elements.h" namespace "Cantera":
    vector[string] elementSymbols()
    vector[string] elementNames()
    double getElementWeight(string ename) except +translate_exception
    double getElementWeight(int atomicNumber) except +translate_exception
    int numElementsDefined()
    int getAtomicNumber(string ename) except +translate_exception
    string getElementSymbol(string ename) except +translate_exception
    string getElementSymbol(int atomicNumber) except +translate_exception
    string getElementName(string ename) except +translate_exception
    string getElementName(int atomicNumber) except +translate_exception


cdef class Species:
    cdef shared_ptr[CxxSpecies] _species
    cdef CxxSpecies* species
    cdef _SolutionBase _phase

    cdef _assign(self, shared_ptr[CxxSpecies] other)

cdef class ThermoPhase(_SolutionBase):
    cdef double _mass_factor(self)
    cdef double _mole_factor(self)
    cpdef int element_index(self, element) except *
    cpdef int species_index(self, species) except *
    cdef np.ndarray _getArray1(self, thermoMethod1d method)
    cdef void _setArray1(self, thermoMethod1d method, values) except *
    cdef CxxPlasmaPhase* plasma
    cdef public object _enable_plasma

cdef class InterfacePhase(ThermoPhase):
    cdef CxxSurfPhase* surf
