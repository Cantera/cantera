# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.
# cython: language_level=3

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.map cimport map as stdmap
from libcpp.cast cimport dynamic_cast
from libcpp.unordered_map cimport unordered_map
from libcpp.pair cimport pair
from libcpp cimport bool as cbool
from libcpp.functional cimport function
from cpython cimport bool as pybool
from cpython.ref cimport PyObject

import numpy as np
cimport numpy as np

ctypedef stdmap[string,double] Composition

cdef extern from "<map>" namespace "std":
    cdef cppclass multimap[T, U]:
        cppclass iterator:
            pair[T, U]& operator*() nogil
            iterator operator++() nogil
            iterator operator--() nogil
            bint operator==(iterator) nogil
            bint operator!=(iterator) nogil
        multimap() nogil except +
        U& operator[](T&) nogil
        iterator begin() nogil
        iterator end() nogil
        pair[iterator, bint] insert(pair[T, U]) nogil
        iterator find(T&) nogil

cdef extern from "<array>" namespace "std" nogil:
    cdef cppclass size_array1 "std::array<size_t, 1>":
        size_array1() except+
        size_t& operator[](size_t)

    cdef cppclass size_array2 "std::array<size_t, 2>":
        size_array2() except+
        size_t& operator[](size_t)

    cdef cppclass size_array3 "std::array<size_t, 3>":
        size_array3() except+
        size_t& operator[](size_t)

cdef extern from "cantera/cython/funcWrapper.h":
    ctypedef double (*callback_wrapper)(double, void*, void**)
    cdef int translate_exception()

    cdef cppclass CxxFunc1 "Func1Py":
        CxxFunc1(callback_wrapper, void*)
        double eval(double) except +translate_exception

    cdef cppclass PyFuncInfo:
        PyFuncInfo()
        PyObject* func()
        void setFunc(PyObject*)
        PyObject* exceptionType()
        void setExceptionType(PyObject*)
        PyObject* exceptionValue()
        void setExceptionValue(PyObject*)

    # pyOverride is actually a templated function, but we have to specify the individual
    # instantiations because Cython doesn't understand variadic templates
    cdef function[void(double)] pyOverride(PyObject*, void(PyFuncInfo&, double))
    cdef function[void(cbool)] pyOverride(PyObject*, void(PyFuncInfo&, cbool))
    cdef function[void()] pyOverride(PyObject*, void(PyFuncInfo&))
    cdef function[void(size_array1, double*)] pyOverride(
        PyObject*, void(PyFuncInfo&, size_array1, double*))
    cdef function[void(size_array1, double, double*)] pyOverride(
        PyObject*, void(PyFuncInfo&, size_array1, double, double*))
    cdef function[void(size_array2, double, double*, double*)] pyOverride(
        PyObject*, void(PyFuncInfo&, size_array2, double, double*, double*))
    cdef function[void(size_array3, double*, double*, double*)] pyOverride(
        PyObject*, void(PyFuncInfo&, size_array3, double*, double*, double*))
    cdef function[int(string&, size_t)] pyOverride(PyObject*, int(PyFuncInfo&, string&, size_t))
    cdef function[int(size_t&, const string&)] pyOverride(
        PyObject*, int(PyFuncInfo&, size_t&, const string&))

cdef extern from "cantera/numerics/Func1.h":
    cdef cppclass CxxTabulated1 "Cantera::Tabulated1":
        CxxTabulated1(int, double*, double*, string) except +translate_exception
        double eval(double) except +translate_exception

cdef extern from "cantera/numerics/eigen_sparse.h" namespace "Eigen":
    cdef cppclass CxxSparseMatrix "Eigen::SparseMatrix<double>":
        CxxSparseMatrix()
        size_t nonZeros()
        size_t rows()
        size_t cols()
        size_t outerSize()

cdef extern from "cantera/base/xml.h" namespace "Cantera":
    cdef cppclass XML_Node:
        XML_Node* findByName(string)
        XML_Node* findID(string)
        int nChildren()

cdef extern from "cantera/base/Units.h" namespace "Cantera":
    cdef cppclass CxxUnits "Cantera::Units":
        CxxUnits()
        CxxUnits(CxxUnits)
        CxxUnits(string, cbool) except +translate_exception
        string str()
        double factor()

    cdef cppclass CxxUnitSystem "Cantera::UnitSystem":
        CxxUnitSystem()
        stdmap[string, string] defaults()
        void setDefaults(stdmap[string, string]&) except +translate_exception

cdef extern from "cantera/base/AnyMap.h" namespace "Cantera":
    cdef cppclass CxxAnyValue "Cantera::AnyValue"

    cdef cppclass CxxAnyMap "Cantera::AnyMap":
        cppclass Iterator:
            pair[string, CxxAnyValue]& operator*()
            Iterator& operator++()
            cbool operator!=(Iterator&)

        cppclass OrderedIterator:
            pair[string, CxxAnyValue]& operator*()
            OrderedIterator& operator++()
            cbool operator!=(OrderedIterator&)

        cppclass OrderedProxy:
            OrderedIterator begin()
            OrderedIterator end()

        CxxAnyMap()
        Iterator begin()
        Iterator end()
        OrderedProxy ordered() except +translate_exception
        CxxAnyValue& operator[](string) except +translate_exception
        cbool empty()
        cbool hasKey(string)
        void clear()
        void update(CxxAnyMap& other, cbool)
        string keys_str()
        void applyUnits()

    cdef cppclass CxxAnyValue "Cantera::AnyValue":
        CxxAnyValue()
        CxxAnyValue& operator=(string) except +translate_exception
        CxxAnyValue& operator=(double) except +translate_exception
        CxxAnyValue& operator=(cbool) except +translate_exception
        CxxAnyValue& operator=(int) except +translate_exception
        CxxAnyValue& operator=(long) except +translate_exception
        CxxAnyValue& operator=(CxxAnyMap) except +translate_exception
        CxxAnyValue& operator=[T](vector[T]) except +translate_exception
        unordered_map[string, CxxAnyMap*] asMap(string) except +translate_exception
        CxxAnyMap& getMapWhere(string, string) except +translate_exception
        T& asType "as" [T]() except +translate_exception
        string type_str()
        cbool empty()
        cbool isType "is" [T]()
        cbool isScalar()

    CxxAnyMap AnyMapFromYamlFile "Cantera::AnyMap::fromYamlFile" (string) except +translate_exception
    CxxAnyMap AnyMapFromYamlString "Cantera::AnyMap::fromYamlString" (string) except +translate_exception

cdef extern from "cantera/base/stringUtils.h" namespace "Cantera":
    cdef Composition parseCompString(string) except +translate_exception

cdef extern from "cantera/base/global.h" namespace "Cantera":
    cdef void CxxAddDirectory "Cantera::addDirectory" (string)
    cdef string CxxGetDataDirectories "Cantera::getDataDirectories" (string)
    cdef size_t CxxNpos "Cantera::npos"
    cdef void CxxAppdelete "Cantera::appdelete" ()
    cdef XML_Node* CxxGetXmlFile "Cantera::get_XML_File" (string) except +translate_exception
    cdef XML_Node* CxxGetXmlFromString "Cantera::get_XML_from_string" (string) except +translate_exception
    cdef void Cxx_make_deprecation_warnings_fatal "Cantera::make_deprecation_warnings_fatal" ()
    cdef void Cxx_suppress_deprecation_warnings "Cantera::suppress_deprecation_warnings" ()
    cdef void Cxx_suppress_thermo_warnings "Cantera::suppress_thermo_warnings" (cbool)
    cdef void Cxx_use_legacy_rate_constants "Cantera::use_legacy_rate_constants" (cbool)
    cdef string CxxGitCommit "Cantera::gitCommit" ()
    cdef cbool CxxDebugModeEnabled "Cantera::debugModeEnabled" ()

cdef extern from "<memory>":
    cppclass shared_ptr "std::shared_ptr" [T]:
        T* get()
        void reset(T*)

cdef extern from "cantera/base/Array.h" namespace "Cantera":
    cdef cppclass CxxArray2D "Cantera::Array2D":
        CxxArray2D()
        CxxArray2D(size_t, size_t)
        void resize(size_t, size_t)
        double operator()(size_t, size_t)
        vector[double]& data()
        size_t nRows()
        size_t nColumns()

cdef extern from "cantera/base/Delegator.h" namespace "Cantera":
    cdef cppclass CxxDelegator "Cantera::Delegator":
        Delegator()

        void setDelegate(string&, function[void()], string&) except +translate_exception
        void setDelegate(string&, function[void(cbool)], string&) except +translate_exception
        void setDelegate(string&, function[void(double)], string&) except +translate_exception
        void setDelegate(string&, function[void(size_array1, double*)], string&) except +translate_exception
        void setDelegate(string&, function[void(size_array1, double, double*)], string&) except +translate_exception
        void setDelegate(string&, function[void(size_array2, double, double*, double*)], string&) except +translate_exception
        void setDelegate(string&, function[void(size_array3, double*, double*, double*)], string&) except +translate_exception
        void setDelegate(string&, function[int(string&, size_t)], string&) except +translate_exception
        void setDelegate(string&, function[int(size_t&, string&)], string&) except +translate_exception

ctypedef CxxDelegator* CxxDelegatorPtr

cdef extern from "cantera/thermo/SpeciesThermoInterpType.h":
    cdef cppclass CxxSpeciesThermo "Cantera::SpeciesThermoInterpType":
        CxxSpeciesThermo()
        int reportType()
        void updatePropertiesTemp(double, double*, double*, double*) except +translate_exception
        double minTemp()
        double maxTemp()
        double refPressure()
        void reportParameters(size_t&, int&, double&, double&, double&, double* const) except +translate_exception
        int nCoeffs() except +translate_exception
        CxxAnyMap parameters(cbool) except +translate_exception
        CxxAnyMap& input()

cdef extern from "cantera/thermo/SpeciesThermoFactory.h":
    cdef CxxSpeciesThermo* CxxNewSpeciesThermo "Cantera::newSpeciesThermoInterpType"\
        (int, double, double, double, double*) except +translate_exception

cdef extern from "cantera/thermo/Species.h" namespace "Cantera":
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
        CxxAnyMap parameters(CxxThermoPhase*) except +translate_exception
        CxxAnyMap input

    cdef shared_ptr[CxxSpecies] CxxNewSpecies "newSpecies" (XML_Node&)
    cdef vector[shared_ptr[CxxSpecies]] CxxGetSpecies "getSpecies" (XML_Node&)
    cdef shared_ptr[CxxSpecies] CxxNewSpecies "newSpecies" (CxxAnyMap&) except +translate_exception
    cdef vector[shared_ptr[CxxSpecies]] CxxGetSpecies "getSpecies" (CxxAnyValue&) except +translate_exception


cdef extern from "cantera/base/Solution.h" namespace "Cantera":
    cdef cppclass CxxSolution "Cantera::Solution":
        CxxSolution()
        string name()
        void setName(string)
        string source()
        void setSource(string)
        CxxAnyMap& header() except +translate_exception
        shared_ptr[CxxThermoPhase] thermo()
        void setThermo(shared_ptr[CxxThermoPhase])
        shared_ptr[CxxKinetics] kinetics()
        void setKinetics(shared_ptr[CxxKinetics])
        shared_ptr[CxxTransport] transport()
        void setTransport(shared_ptr[CxxTransport])
        CxxAnyMap parameters(cbool) except +translate_exception
        size_t nAdjacent()
        shared_ptr[CxxSolution] adjacent(size_t)

    cdef shared_ptr[CxxSolution] CxxNewSolution "Cantera::Solution::create" ()
    cdef shared_ptr[CxxSolution] newSolution (
        string, string, string, vector[shared_ptr[CxxSolution]]) except +translate_exception
    cdef shared_ptr[CxxSolution] newSolution (
        CxxAnyMap&, CxxAnyMap&, string, vector[shared_ptr[CxxSolution]]) except +translate_exception


cdef extern from "cantera/base/Interface.h" namespace "Cantera":
    cdef shared_ptr[CxxSolution] newInterface(string, string) except +translate_exception
    cdef shared_ptr[CxxSolution] newInterface(
        string, string, vector[string]) except +translate_exception


cdef extern from "cantera/thermo/ThermoPhase.h" namespace "Cantera":
    ctypedef enum ThermoBasis:
        mass "Cantera::ThermoBasis::mass",
        molar "Cantera::ThermoBasis::molar"

cdef extern from "cantera/thermo/ThermoPhase.h" namespace "Cantera":
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
        void setState_TR(double, double) except +translate_exception
        void setState_TP(double, double) except +translate_exception
        void setState_HP(double, double) except +translate_exception
        void setState_UV(double, double) except +translate_exception
        void setState_SP(double, double) except +translate_exception
        void setState_SV(double, double) except +translate_exception
        void setState_RP(double, double) except +translate_exception
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

cdef extern from "cantera/thermo/IdealGasPhase.h":
    cdef cppclass CxxIdealGasPhase "Cantera::IdealGasPhase"


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
    cdef cppclass CxxPlasmaPhase "Cantera::PlasmaPhase":
        CxxPlasmaPhase()
        double electronTemperature() except +translate_exception
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
        void enableNormalizeElectronEnergyDist(cbool) except +translate_exception
        cbool normalizeElectronEnergyDistEnabled()
        void setMeanElectronEnergy(double) except +translate_exception
        double isotropicShapeFactor()
        double meanElectronEnergy()
        size_t nElectronEnergyLevels() except +translate_exception

cdef extern from "cantera/kinetics/ReactionRateFactory.h" namespace "Cantera":
    cdef shared_ptr[CxxReactionRate] CxxNewReactionRate "newReactionRate" (string) except +translate_exception
    cdef shared_ptr[CxxReactionRate] CxxNewReactionRate "newReactionRate" (CxxAnyMap&) except +translate_exception

cdef extern from "cantera/kinetics/ReactionRate.h" namespace "Cantera":
    cdef cppclass CxxReactionRate "Cantera::ReactionRate":
        CxxReactionRate()
        string type()
        double eval(double) except +translate_exception
        double eval(double, double) except +translate_exception
        double eval(double, vector[double]&) except +translate_exception
        CxxAnyMap parameters() except +translate_exception

cdef extern from "cantera/kinetics/Arrhenius.h" namespace "Cantera":
    cdef cppclass CxxArrheniusBase "Cantera::ArrheniusBase" (CxxReactionRate):
        CxxArrheniusBase()
        double preExponentialFactor()
        double temperatureExponent()
        double activationEnergy()
        cbool allowNegativePreExponentialFactor()
        void setAllowNegativePreExponentialFactor(bool)

    cdef cppclass CxxArrheniusRate "Cantera::ArrheniusRate" (CxxArrheniusBase):
        CxxArrheniusRate(double, double, double)
        CxxArrheniusRate(CxxAnyMap) except +translate_exception
        double evalRate(double, double)

    cdef cppclass CxxArrhenius2 "Cantera::Arrhenius2" (CxxArrheniusRate):
        CxxArrhenius2(double, double, double)

    cdef cppclass CxxBlowersMasel "Cantera::BlowersMasel" (CxxArrheniusBase):
        CxxBlowersMasel(double, double, double, double)
        double evalRate(double, double)
        double bondEnergy()
        double deltaH()
        void setDeltaH(double)

    cdef cppclass CxxBlowersMaselRate "Cantera::BlowersMaselRate" (CxxBlowersMasel):
        CxxBlowersMaselRate(CxxAnyMap) except +translate_exception
        CxxBlowersMaselRate(double, double, double, double)

cdef extern from "cantera/kinetics/TwoTempPlasmaRate.h" namespace "Cantera":
    cdef cppclass CxxTwoTempPlasmaRate "Cantera::TwoTempPlasmaRate" (CxxArrheniusBase):
        CxxTwoTempPlasmaRate(CxxAnyMap) except +translate_exception
        CxxTwoTempPlasmaRate(double, double, double, double)
        double activationElectronEnergy()

cdef extern from "cantera/kinetics/Reaction.h" namespace "Cantera":
    cdef shared_ptr[CxxReaction] CxxNewReaction "Cantera::newReaction" (string) except +translate_exception
    cdef shared_ptr[CxxReaction] CxxNewReaction "newReaction" (XML_Node&) except +translate_exception
    cdef shared_ptr[CxxReaction] CxxNewReaction "newReaction" (CxxAnyMap&, CxxKinetics&) except +translate_exception
    cdef vector[shared_ptr[CxxReaction]] CxxGetReactions "getReactions" (XML_Node&) except +translate_exception
    cdef vector[shared_ptr[CxxReaction]] CxxGetReactions "getReactions" (CxxAnyValue&, CxxKinetics&) except +translate_exception

    cdef cppclass CxxFalloffRate "Cantera::FalloffRate" (CxxReactionRate):
        CxxFalloffRate()
        CxxFalloffRate(CxxAnyMap) except +translate_exception
        cbool allowNegativePreExponentialFactor()
        void setAllowNegativePreExponentialFactor(bool)
        cbool chemicallyActivated()
        void setChemicallyActivated(bool)
        CxxArrheniusRate& lowRate()
        void setLowRate(CxxArrheniusRate&) except +translate_exception
        CxxArrheniusRate& highRate()
        void setHighRate(CxxArrheniusRate&) except +translate_exception
        void getFalloffCoeffs(vector[double]&)
        void setFalloffCoeffs(vector[double]&) except +translate_exception
        double evalF(double, double) except +translate_exception

    cdef cppclass CxxLindemannRate "Cantera::LindemannRate" (CxxFalloffRate):
        CxxLindemannRate(CxxAnyMap) except +translate_exception

    cdef cppclass CxxTroeRate "Cantera::TroeRate" (CxxFalloffRate):
        CxxTroeRate(CxxAnyMap) except +translate_exception

    cdef cppclass CxxSriRate "Cantera::SriRate" (CxxFalloffRate):
        CxxSriRate(CxxAnyMap) except +translate_exception

    cdef cppclass CxxTsangRate "Cantera::TsangRate" (CxxFalloffRate):
        CxxTsangRate(CxxAnyMap) except +translate_exception

    cdef cppclass CxxPlogRate "Cantera::PlogRate" (CxxReactionRate):
        CxxPlogRate()
        CxxPlogRate(CxxAnyMap) except +translate_exception
        CxxPlogRate(multimap[double, CxxArrheniusRate])
        multimap[double, CxxArrheniusRate] getRates()

    cdef cppclass CxxChebyshevRate "Cantera::ChebyshevRate" (CxxReactionRate):
        CxxChebyshevRate()
        CxxChebyshevRate(CxxAnyMap) except +translate_exception
        CxxChebyshevRate(double, double, double, double, CxxArray2D)
        double Tmin()
        double Tmax()
        double Pmin()
        double Pmax()
        size_t nTemperature()
        size_t nPressure()
        CxxArray2D& data()

    cdef cppclass CxxCustomFunc1Rate "Cantera::CustomFunc1Rate" (CxxReactionRate):
        CxxCustomFunc1Rate()
        void setRateFunction(shared_ptr[CxxFunc1]) except +translate_exception

    cdef cppclass CxxInterfaceRateBase "Cantera::InterfaceRateBase":
        void getCoverageDependencies(CxxAnyMap)
        void setCoverageDependencies(CxxAnyMap) except +translate_exception
        void setSpecies(vector[string]&) except +translate_exception
        double siteDensity()
        void setSiteDensity(double)
        cbool usesElectrochemistry()
        double beta()

    cdef cppclass CxxStickingCoverage "Cantera::StickingCoverage" (CxxInterfaceRateBase):
        cbool motzWiseCorrection()
        void setMotzWiseCorrection(cbool)
        string stickingSpecies()
        void setStickingSpecies(string)
        double stickingOrder()
        void setStickingOrder(double)
        double stickingWeight()
        void setStickingWeight(double)

    cdef cppclass CxxInterfaceArrheniusRate "Cantera::InterfaceArrheniusRate" (CxxReactionRate, CxxArrheniusRate, CxxInterfaceRateBase):
        CxxInterfaceArrheniusRate()
        CxxInterfaceArrheniusRate(CxxAnyMap) except +translate_exception
        CxxInterfaceArrheniusRate(double, double, double)

    cdef cppclass CxxStickingArrheniusRate "Cantera::StickingArrheniusRate" (CxxReactionRate, CxxArrheniusRate, CxxStickingCoverage):
        CxxStickingArrheniusRate()
        CxxStickingArrheniusRate(CxxAnyMap) except +translate_exception
        CxxStickingArrheniusRate(double, double, double)

    cdef cppclass CxxInterfaceBlowersMaselRate "Cantera::InterfaceBlowersMaselRate" (CxxReactionRate, CxxBlowersMasel, CxxInterfaceRateBase):
        CxxInterfaceBlowersMaselRate()
        CxxInterfaceBlowersMaselRate(CxxAnyMap) except +translate_exception
        CxxInterfaceBlowersMaselRate(double, double, double, double)

    cdef cppclass CxxStickingBlowersMaselRate "Cantera::StickingBlowersMaselRate" (CxxReactionRate, CxxBlowersMasel, CxxStickingCoverage):
        CxxStickingBlowersMaselRate()
        CxxStickingBlowersMaselRate(CxxAnyMap) except +translate_exception
        CxxStickingBlowersMaselRate(double, double, double, double)

    cdef cppclass CxxThirdBody "Cantera::ThirdBody":
        CxxThirdBody()
        CxxThirdBody(double)
        double efficiency(string)
        Composition efficiencies
        double default_efficiency

    cdef cppclass CxxReaction "Cantera::Reaction":
        CxxReaction()

        string reactantString()
        string productString()
        string equation()
        void setEquation(const string&) except +translate_exception
        string type()
        void validate() except +translate_exception
        CxxAnyMap parameters(cbool) except +translate_exception
        CxxAnyMap input
        int reaction_type
        Composition reactants
        Composition products
        Composition orders
        string id
        cbool reversible
        cbool duplicate
        cbool allow_nonreactant_orders
        cbool allow_negative_orders
        shared_ptr[CxxThirdBody] thirdBody()
        cbool usesLegacy()
        CxxUnits rate_units

        shared_ptr[CxxReactionRate] rate()
        void setRate(shared_ptr[CxxReactionRate])

    cdef cppclass CxxElementaryReaction2 "Cantera::ElementaryReaction2" (CxxReaction):
        CxxElementaryReaction2()
        CxxArrhenius2 rate
        cbool allow_negative_pre_exponential_factor

    cdef cppclass CxxThreeBodyReaction2 "Cantera::ThreeBodyReaction2" (CxxElementaryReaction2):
        CxxThreeBodyReaction2()
        CxxThirdBody third_body

    cdef cppclass CxxFalloff "Cantera::FalloffRate":
        CxxFalloff()
        void updateTemp(double, double*)
        double F(double, double*)
        size_t workSize()

        size_t nParameters()
        string type()
        void getParameters(double*)

    cdef cppclass CxxFalloffReaction2 "Cantera::FalloffReaction2" (CxxReaction):
        CxxFalloffReaction2()

        CxxArrhenius2 low_rate
        CxxArrhenius2 high_rate
        CxxThirdBody third_body
        shared_ptr[CxxFalloff] falloff
        cbool allow_negative_pre_exponential_factor

    cdef cppclass CxxChemicallyActivatedReaction "Cantera::ChemicallyActivatedReaction" (CxxFalloffReaction2):
        CxxChemicallyActivatedReaction()

    cdef cppclass CxxPlog "Cantera::Plog":
        CxxPlog(multimap[double,CxxArrhenius2])
        vector[pair[double, CxxArrhenius2]] rates()
        void update_C(double*)
        double updateRC(double, double)

    cdef cppclass CxxPlogReaction2 "Cantera::PlogReaction2" (CxxReaction):
        CxxPlog rate

    cdef cppclass CxxChebyshev "Cantera::ChebyshevRate":
        CxxChebyshev(double, double, double, double, CxxArray2D)
        double Tmin()
        double Tmax()
        double Pmin()
        double Pmax()
        size_t nTemperature()
        size_t nPressure()
        CxxArray2D& data()
        void update_C(double*)
        double updateRC(double, double)

    cdef cppclass CxxChebyshevReaction2 "Cantera::ChebyshevReaction2" (CxxReaction):
        CxxChebyshev rate

    cdef cppclass CxxCoverageDependency "Cantera::CoverageDependency":
        CxxCoverageDependency(double, double, double)
        double a
        double E
        double m

    cdef cppclass CxxInterfaceReaction2 "Cantera::InterfaceReaction2" (CxxElementaryReaction2):
        stdmap[string, CxxCoverageDependency] coverage_deps
        cbool is_sticking_coefficient
        cbool use_motz_wise_correction
        string sticking_species

    cdef cppclass CxxThreeBodyReaction3 "Cantera::ThreeBodyReaction3" (CxxReaction):
        CxxThreeBodyReaction3()

    cdef cppclass CxxFalloffReaction3 "Cantera::FalloffReaction3" (CxxReaction):
        CxxFalloffReaction3()

    cdef cppclass CxxCustomFunc1Reaction "Cantera::CustomFunc1Reaction" (CxxReaction):
        CxxCustomFunc1Reaction()

cdef extern from "cantera/kinetics/FalloffFactory.h" namespace "Cantera":
    cdef shared_ptr[CxxFalloff] CxxNewFalloff "Cantera::newFalloff" (string, vector[double]) except +translate_exception

cdef extern from "cantera/kinetics/Kinetics.h" namespace "Cantera":
    cdef cppclass CxxKinetics "Cantera::Kinetics":
        CxxKinetics()
        string kineticsType()
        int nTotalSpecies()
        int nReactions()
        int nPhases()
        int reactionPhaseIndex()
        int phaseIndex(string)
        int kineticsSpeciesIndex(int, int)
        int kineticsSpeciesIndex(string)
        string kineticsSpeciesName(int)

        CxxThermoPhase& thermo(int)

        void addPhase(CxxThermoPhase&) except +translate_exception
        void init() except +translate_exception
        void skipUndeclaredThirdBodies(cbool)
        void addReaction(shared_ptr[CxxReaction]) except +translate_exception
        void addReaction(shared_ptr[CxxReaction], cbool) except +translate_exception
        void modifyReaction(int, shared_ptr[CxxReaction]) except +translate_exception
        void invalidateCache() except +translate_exception
        void resizeReactions()

        shared_ptr[CxxReaction] reaction(size_t) except +translate_exception
        double reactantStoichCoeff(int, int) except +translate_exception
        double productStoichCoeff(int, int) except +translate_exception

        double multiplier(int)
        void setMultiplier(int, double)

        void getDerivativeSettings(CxxAnyMap&) except +translate_exception
        void setDerivativeSettings(CxxAnyMap&) except +translate_exception


cdef extern from "cantera/kinetics/InterfaceKinetics.h":
    cdef cppclass CxxInterfaceKinetics "Cantera::InterfaceKinetics":
        void advanceCoverages(double, double, double, double, size_t, size_t) except +translate_exception
        void solvePseudoSteadyStateProblem() except +translate_exception


cdef extern from "cantera/transport/TransportBase.h" namespace "Cantera":
    cdef cppclass CxxTransport "Cantera::Transport":
        CxxTransport(CxxThermoPhase*)
        string transportType()
        cbool CKMode() except +translate_exception
        double viscosity() except +translate_exception
        double thermalConductivity() except +translate_exception
        double electricalConductivity() except +translate_exception
        void getSpeciesViscosities(double*) except +translate_exception
        void getCollisionIntegralPolynomial(size_t i, size_t j, double* dataA, double* dataB, double* dataC) except +translate_exception
        void setCollisionIntegralPolynomial(size_t i, size_t j, double* dataA, double* dataB, double* dataC, cbool flag) except +translate_exception


cdef extern from "cantera/transport/DustyGasTransport.h" namespace "Cantera":
    cdef cppclass CxxDustyGasTransport "Cantera::DustyGasTransport":
        void setPorosity(double) except +translate_exception
        void setTortuosity(double) except +translate_exception
        void setMeanPoreRadius(double) except +translate_exception
        void setMeanParticleDiameter(double) except +translate_exception
        void setPermeability(double) except +translate_exception
        void getMolarFluxes(double*, double*, double, double*) except +translate_exception
        CxxTransport& gasTransport() except +translate_exception


cdef extern from "cantera/transport/TransportData.h" namespace "Cantera":
    cdef cppclass CxxTransportData "Cantera::TransportData":
        CxxTransportData()
        CxxAnyMap parameters(cbool) except +translate_exception
        CxxAnyMap input

    cdef cppclass CxxGasTransportData "Cantera::GasTransportData" (CxxTransportData):
        CxxGasTransportData()
        CxxGasTransportData(string, double, double, double, double, double, double, double, double)
        void setCustomaryUnits(string, double, double, double, double, double, double, double, double)

        string geometry
        double diameter
        double well_depth
        double dipole
        double polarizability
        double rotational_relaxation
        double acentric_factor
        double dispersion_coefficient
        double quadrupole_polarizability

cdef extern from "cantera/base/YamlWriter.h" namespace "Cantera":
    cdef cppclass CxxYamlWriter "Cantera::YamlWriter":
        CxxYamlWriter()
        void setHeader(CxxAnyMap) except +translate_exception
        void addPhase(shared_ptr[CxxSolution]) except +translate_exception
        string toYamlString() except +translate_exception
        void toYamlFile(string&) except +translate_exception
        void setPrecision(int)
        void skipUserDefined(cbool)
        void setUnitSystem(CxxUnitSystem&) except +translate_exception

cdef extern from "cantera/equil/MultiPhase.h" namespace "Cantera":
    cdef cppclass CxxMultiPhase "Cantera::MultiPhase":
        CxxMultiPhase()
        void addPhase(CxxThermoPhase*, double) except +translate_exception
        void init() except +translate_exception
        void updatePhases() except +translate_exception

        void equilibrate(string, string, double, int, int, int, int) except +translate_exception

        size_t nSpecies()
        size_t nElements()
        size_t nPhases()
        size_t elementIndex(string) except +translate_exception
        size_t speciesIndex(size_t, size_t) except +translate_exception
        string speciesName(size_t) except +translate_exception
        double nAtoms(size_t, size_t) except +translate_exception

        double phaseMoles(size_t) except +translate_exception
        void setPhaseMoles(size_t, double) except +translate_exception
        void setMoles(double*) except +translate_exception
        void setMolesByName(string) except +translate_exception

        double speciesMoles(size_t) except +translate_exception
        double elementMoles(size_t) except +translate_exception

        void setTemperature(double) except +translate_exception
        double temperature()
        void setPressure(double) except +translate_exception
        double pressure()

        double minTemp() except +translate_exception
        double maxTemp() except +translate_exception
        double charge() except +translate_exception
        double phaseCharge(size_t) except +translate_exception
        void getChemPotentials(double*) except +translate_exception
        double enthalpy() except +translate_exception
        double entropy() except +translate_exception
        double gibbs() except +translate_exception
        double cp() except +translate_exception
        double volume() except +translate_exception


cdef extern from "cantera/zerodim.h" namespace "Cantera":

    cdef cppclass CxxWall "Cantera::Wall"
    cdef cppclass CxxReactorSurface "Cantera::ReactorSurface"
    cdef cppclass CxxFlowDevice "Cantera::FlowDevice"

    # factories

    cdef CxxReactorBase* newReactor(string) except +translate_exception
    cdef CxxFlowDevice* newFlowDevice(string) except +translate_exception
    cdef CxxWallBase* newWall(string) except +translate_exception

    # reactors

    cdef cppclass CxxReactorBase "Cantera::ReactorBase":
        CxxReactorBase()
        string type()
        void setThermoMgr(CxxThermoPhase&) except +translate_exception
        void restoreState() except +translate_exception
        void syncState() except +translate_exception
        double volume()
        string name()
        void setName(string)
        void setInitialVolume(double)

    cdef cppclass CxxReactor "Cantera::Reactor" (CxxReactorBase):
        CxxReactor()
        void setKineticsMgr(CxxKinetics&)
        void setChemistry(cbool)
        cbool chemistryEnabled()
        void setEnergy(int)
        cbool energyEnabled()
        size_t componentIndex(string&) except +translate_exception
        string componentName(size_t) except +translate_exception
        size_t neq()
        void getState(double*) except +translate_exception
        void addSurface(CxxReactorSurface*)
        void setAdvanceLimit(string&, double) except +translate_exception
        void addSensitivityReaction(size_t) except +translate_exception
        void addSensitivitySpeciesEnthalpy(size_t) except +translate_exception
        size_t nSensParams()

    cdef cppclass CxxFlowReactor "Cantera::FlowReactor" (CxxReactor):
        CxxFlowReactor()
        void setMassFlowRate(double) except +translate_exception
        double speed()
        double distance()


    # walls

    cdef cppclass CxxWallBase "Cantera::WallBase":
        CxxWallBase()
        string type()
        cbool install(CxxReactorBase&, CxxReactorBase&)
        double area()
        void setArea(double)
        void setKinetics(CxxKinetics*, CxxKinetics*)
        void setCoverages(int, double*)
        void setCoverages(int, Composition&) except +translate_exception
        void syncCoverages(int)
        double vdot(double)
        double Q(double)

        void addSensitivityReaction(int, size_t) except +translate_exception
        size_t nSensParams(int)

    cdef cppclass CxxWall "Cantera::Wall" (CxxWallBase):
        CxxWall()
        void setExpansionRateCoeff(double)
        double getExpansionRateCoeff()
        void setHeatTransferCoeff(double)
        double getHeatTransferCoeff()
        void setEmissivity(double) except +translate_exception
        double getEmissivity()
        void setVelocity(CxxFunc1*)
        void setHeatFlux(CxxFunc1*)

    # reactor surface

    cdef cppclass CxxReactorSurface "Cantera::ReactorSurface":
        CxxReactorSurface()
        double area()
        void setArea(double)
        void setKinetics(CxxKinetics*)
        void setCoverages(double*)
        void setCoverages(Composition&) except +translate_exception
        void syncState()
        void addSensitivityReaction(size_t) except +translate_exception
        size_t nSensParams()

    # flow devices

    cdef cppclass CxxFlowDevice "Cantera::FlowDevice":
        CxxFlowDevice()
        string type()
        double massFlowRate() except +translate_exception
        double massFlowRate(double) except +translate_exception
        cbool install(CxxReactorBase&, CxxReactorBase&) except +translate_exception
        void setPressureFunction(CxxFunc1*) except +translate_exception
        void setTimeFunction(CxxFunc1*) except +translate_exception

    cdef cppclass CxxMassFlowController "Cantera::MassFlowController" (CxxFlowDevice):
        CxxMassFlowController()
        void setMassFlowRate(double)
        void setMassFlowCoeff(double)
        double getMassFlowCoeff()

    cdef cppclass CxxValve "Cantera::Valve" (CxxFlowDevice):
        CxxValve()
        double getValveCoeff()
        void setValveCoeff(double)

    cdef cppclass CxxPressureController "Cantera::PressureController" (CxxFlowDevice):
        CxxPressureController()
        double getPressureCoeff()
        void setPressureCoeff(double)
        void setMaster(CxxFlowDevice*)

    # reactor net

    cdef cppclass CxxReactorNet "Cantera::ReactorNet":
        CxxReactorNet()
        void addReactor(CxxReactor&)
        double advance(double, cbool) except +translate_exception
        double step() except +translate_exception
        void initialize() except +translate_exception
        void reinitialize() except +translate_exception
        double time()
        void setInitialTime(double)
        void setTolerances(double, double)
        double rtol()
        double atol()
        double maxTimeStep()
        void setMaxTimeStep(double)
        void setMaxErrTestFails(int)
        void setMaxSteps(int)
        int maxSteps()
        cbool verbose()
        void setVerbose(cbool)
        size_t neq()
        void getState(double*)
        void getDerivative(int, double *) except +translate_exception
        void setAdvanceLimits(double*)
        cbool getAdvanceLimits(double*)
        string componentName(size_t) except +translate_exception
        size_t globalComponentIndex(string&, int) except +translate_exception
        void setSensitivityTolerances(double, double)
        double rtolSensitivity()
        double atolSensitivity()
        double sensitivity(size_t, size_t) except +translate_exception
        double sensitivity(string&, size_t, int) except +translate_exception
        size_t nparams()
        string sensitivityParameterName(size_t) except +translate_exception

cdef extern from "cantera/zeroD/ReactorDelegator.h" namespace "Cantera":
    cdef cppclass CxxReactorAccessor "Cantera::ReactorAccessor":
        CxxReactorAccessor()
        void setNEq(size_t)
        double vdot()
        void setVdot(double)
        double qdot()
        void setQdot(double)
        void restoreThermoState() except +translate_exception
        void restoreSurfaceState(size_t) except +translate_exception

ctypedef CxxReactorAccessor* CxxReactorAccessorPtr

cdef extern from "cantera/thermo/ThermoFactory.h" namespace "Cantera":
    cdef CxxThermoPhase* newPhase(string, string) except +translate_exception
    cdef CxxThermoPhase* newPhase(XML_Node&) except +translate_exception
    cdef shared_ptr[CxxThermoPhase] newPhase(CxxAnyMap&, CxxAnyMap&) except +translate_exception
    cdef CxxThermoPhase* newThermoPhase(string) except +translate_exception
    cdef shared_ptr[CxxThermoPhase] newThermo(string) except +translate_exception

cdef extern from "cantera/kinetics/KineticsFactory.h" namespace "Cantera":
    cdef CxxKinetics* newKineticsMgr(XML_Node&, vector[CxxThermoPhase*]) except +translate_exception
    cdef shared_ptr[CxxKinetics] newKinetics(vector[CxxThermoPhase*], CxxAnyMap&, CxxAnyMap&) except +translate_exception
    cdef CxxKinetics* CxxNewKinetics "Cantera::newKineticsMgr" (string) except +translate_exception
    cdef shared_ptr[CxxKinetics] newKinetics (string) except +translate_exception

cdef extern from "cantera/transport/TransportFactory.h" namespace "Cantera":
    cdef CxxTransport* newDefaultTransportMgr(CxxThermoPhase*) except +translate_exception
    cdef CxxTransport* newTransportMgr(string, CxxThermoPhase*) except +translate_exception
    cdef CxxTransport* newTransportMgr(string) except +translate_exception
    cdef shared_ptr[CxxTransport] newTransport(CxxThermoPhase*, string) except +translate_exception

cdef extern from "cantera/oneD/Domain1D.h":
    cdef cppclass CxxDomain1D "Cantera::Domain1D":
        size_t domainIndex()
        size_t nComponents()
        size_t nPoints()
        string componentName(size_t) except +translate_exception
        size_t componentIndex(string) except +translate_exception
        void setBounds(size_t, double, double)
        double upperBound(size_t)
        double lowerBound(size_t)
        void setTransientTolerances(double, double)
        void setTransientTolerances(double, double, size_t)
        void setSteadyTolerances(double, double)
        void setSteadyTolerances(double, double, size_t)
        double rtol(size_t)
        double atol(size_t)
        double steady_rtol(size_t)
        double steady_atol(size_t)
        double transient_rtol(size_t)
        double transient_atol(size_t)
        double grid(size_t)
        void setupGrid(size_t, double*) except +translate_exception
        void setID(string)
        string& id()


cdef extern from "cantera/oneD/Boundary1D.h":
    cdef cppclass CxxBoundary1D "Cantera::Boundary1D":
        double temperature()
        void setTemperature(double)
        double mdot()
        void setMdot(double)
        size_t nSpecies()
        void setMoleFractions(double*) except +translate_exception
        void setMoleFractions(string) except +translate_exception
        double massFraction(size_t)

    cdef cppclass CxxInlet1D "Cantera::Inlet1D":
        CxxInlet1D()
        double spreadRate()
        void setSpreadRate(double)

    cdef cppclass CxxOutlet1D "Cantera::Outlet1D":
        CxxOutlet1D()

    cdef cppclass CxxOutletRes1D "Cantera::OutletRes1D":
        CxxOutletRes1D()

    cdef cppclass CxxSymm1D "Cantera::Symm1D":
        CxxSymm1D()

    cdef cppclass CxxSurf1D "Cantera::Surf1D":
        CxxSurf1D()

    cdef cppclass CxxReactingSurf1D "Cantera::ReactingSurf1D":
        CxxReactingSurf1D()
        void setKineticsMgr(CxxInterfaceKinetics*) except +translate_exception
        void enableCoverageEquations(cbool) except +translate_exception
        cbool coverageEnabled()


cdef extern from "cantera/oneD/StFlow.h":
    cdef cppclass CxxStFlow "Cantera::StFlow":
        CxxStFlow(CxxIdealGasPhase*, int, int)
        void setKinetics(CxxKinetics&) except +translate_exception
        void setTransport(CxxTransport&, cbool) except +translate_exception
        void setTransport(CxxTransport&) except +translate_exception
        void setPressure(double)
        void enableRadiation(cbool)
        cbool radiationEnabled()
        double radiativeHeatLoss(size_t)
        double pressure()
        void setFixedTempProfile(vector[double]&, vector[double]&)
        void setBoundaryEmissivities(double, double)
        double leftEmissivity()
        double rightEmissivity()
        void solveEnergyEqn()
        void fixTemperature()
        cbool doEnergy(size_t)
        void enableSoret(cbool) except +translate_exception
        cbool withSoret()
        void setFreeFlow()
        void setAxisymmetricFlow()
        string flowType()


cdef extern from "cantera/oneD/IonFlow.h":
    cdef cppclass CxxIonFlow "Cantera::IonFlow":
        CxxIonFlow(CxxIdealGasPhase*, int, int)
        void setSolvingStage(int)
        void solveElectricField()
        void fixElectricField()
        cbool doElectricField(size_t)


cdef extern from "cantera/oneD/Sim1D.h":
    cdef cppclass CxxSim1D "Cantera::Sim1D":
        CxxSim1D(vector[CxxDomain1D*]&) except +translate_exception
        void setValue(size_t, size_t, size_t, double) except +translate_exception
        void setProfile(size_t, size_t, vector[double]&, vector[double]&) except +translate_exception
        void setFlatProfile(size_t, size_t, double) except +translate_exception
        void showSolution() except +translate_exception
        void setTimeStep(double, size_t, int*) except +translate_exception
        void restoreTimeSteppingSolution() except +translate_exception
        void restoreSteadySolution() except +translate_exception
        void setMaxTimeStepCount(int)
        int maxTimeStepCount()
        void getInitialSoln() except +translate_exception
        void solve(int, cbool) except +translate_exception
        void refine(int) except +translate_exception
        void setRefineCriteria(size_t, double, double, double, double) except +translate_exception
        vector[double] getRefineCriteria(int) except +translate_exception
        void save(string, string, string, int) except +translate_exception
        void restore(string, string, int) except +translate_exception
        void writeStats(int) except +translate_exception
        void clearStats()
        void resize() except +translate_exception
        vector[size_t]& gridSizeStats()
        vector[double]& jacobianTimeStats()
        vector[double]& evalTimeStats()
        vector[int]& jacobianCountStats()
        vector[int]& evalCountStats()
        vector[int]& timeStepStats()

        int domainIndex(string) except +translate_exception
        double value(size_t, size_t, size_t) except +translate_exception
        double workValue(size_t, size_t, size_t) except +translate_exception
        void eval(double ) except +translate_exception
        size_t size()
        void solveAdjoint(const double*, double*) except +translate_exception
        void getResidual(double, double*) except +translate_exception
        void setJacAge(int, int)
        void setTimeStepFactor(double)
        void setMinTimeStep(double)
        void setMaxTimeStep(double)
        void setMaxGridPoints(int, size_t) except +translate_exception
        size_t maxGridPoints(size_t) except +translate_exception
        void setGridMin(int, double) except +translate_exception
        void setFixedTemperature(double) except +translate_exception
        double fixedTemperature()
        double fixedTemperatureLocation()
        void setInterrupt(CxxFunc1*) except +translate_exception
        void setTimeStepCallback(CxxFunc1*)
        void setSteadyCallback(CxxFunc1*)

cdef extern from "<sstream>":
    cdef cppclass CxxStringStream "std::stringstream":
        string str()

cdef extern from "cantera/kinetics/ReactionPath.h":
    cdef enum CxxFlow_t "flow_t":
        CxxNetFlow "Cantera::NetFlow"
        CxxOneWayFlow "Cantera::OneWayFlow"

    cdef cppclass CxxReactionPathDiagram "Cantera::ReactionPathDiagram":
        cbool show_details
        double threshold
        string bold_color
        string normal_color
        string dashed_color
        string dot_options
        double bold_min
        double dashed_max
        double label_min
        double scale
        double arrow_width
        CxxFlow_t flow_type
        string title
        void setFont(string)
        string m_font
        void add(CxxReactionPathDiagram&) except +translate_exception
        void exportToDot(CxxStringStream&)
        void writeData(CxxStringStream&)
        void displayOnly(size_t)

    cdef cppclass CxxReactionPathBuilder "Cantera::ReactionPathBuilder":
        void init(CxxStringStream&, CxxKinetics&) except +translate_exception
        void build(CxxKinetics&, string&, CxxStringStream&, CxxReactionPathDiagram&, cbool)


cdef extern from "cantera/cython/wrappers.h":
    # config definitions
    cdef string get_cantera_version()
    cdef int get_sundials_version()

    cdef cppclass CxxPythonLogger "PythonLogger":
        pass

    cdef void CxxSetLogger "setLogger" (CxxPythonLogger*)

    cdef size_t CxxSparseTriplets "sparseTriplets" (CxxSparseMatrix, int*, int*, double*, size_t) except +translate_exception
    cdef void CxxSparseCscData "sparseCscData" (CxxSparseMatrix, double*, int*, int*) except +translate_exception

    # workaround for Cython assignment limitations
    cdef void CxxArray2D_set(CxxArray2D, size_t, size_t, double)

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

    # Kinetics per-reaction properties
    cdef void kin_getFwdRatesOfProgress(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getRevRatesOfProgress(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getNetRatesOfProgress(CxxKinetics*, double*) except +translate_exception

    cdef void kin_getFwdRateConstants_ddT(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getFwdRateConstants_ddP(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getFwdRateConstants_ddC(CxxKinetics*, double*) except +translate_exception

    cdef void kin_getFwdRatesOfProgress_ddT(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getRevRatesOfProgress_ddT(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getNetRatesOfProgress_ddT(CxxKinetics*, double*) except +translate_exception

    cdef void kin_getFwdRatesOfProgress_ddP(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getRevRatesOfProgress_ddP(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getNetRatesOfProgress_ddP(CxxKinetics*, double*) except +translate_exception

    cdef void kin_getFwdRatesOfProgress_ddC(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getRevRatesOfProgress_ddC(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getNetRatesOfProgress_ddC(CxxKinetics*, double*) except +translate_exception

    cdef void kin_getEquilibriumConstants(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getFwdRateConstants(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getRevRateConstants(CxxKinetics*, double*) except +translate_exception

    cdef void kin_getDeltaEnthalpy(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getDeltaGibbs(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getDeltaEntropy(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getDeltaSSEnthalpy(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getDeltaSSGibbs(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getDeltaSSEntropy(CxxKinetics*, double*) except +translate_exception

    cdef void kin_getThirdBodyConcentrations(CxxKinetics*, double*) except +translate_exception

    # Kinetics per-species properties
    cdef void kin_getCreationRates(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getDestructionRates(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getNetProductionRates(CxxKinetics*, double*) except +translate_exception

    cdef void kin_getCreationRates_ddT(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getDestructionRates_ddT(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getNetProductionRates_ddT(CxxKinetics*, double*) except +translate_exception

    cdef void kin_getCreationRates_ddP(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getDestructionRates_ddP(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getNetProductionRates_ddP(CxxKinetics*, double*) except +translate_exception

    cdef void kin_getCreationRates_ddC(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getDestructionRates_ddC(CxxKinetics*, double*) except +translate_exception
    cdef void kin_getNetProductionRates_ddC(CxxKinetics*, double*) except +translate_exception

    # Kinetics sparse matrices
    cdef CxxSparseMatrix kin_reactantStoichCoeffs(CxxKinetics*) except +translate_exception
    cdef CxxSparseMatrix kin_productStoichCoeffs(CxxKinetics*) except +translate_exception
    cdef CxxSparseMatrix kin_revProductStoichCoeffs(CxxKinetics*) except +translate_exception

    cdef CxxSparseMatrix kin_fwdRatesOfProgress_ddX(CxxKinetics*) except +translate_exception
    cdef CxxSparseMatrix kin_revRatesOfProgress_ddX(CxxKinetics*) except +translate_exception
    cdef CxxSparseMatrix kin_netRatesOfProgress_ddX(CxxKinetics*) except +translate_exception

    cdef CxxSparseMatrix kin_creationRates_ddX(CxxKinetics*) except +translate_exception
    cdef CxxSparseMatrix kin_destructionRates_ddX(CxxKinetics*) except +translate_exception
    cdef CxxSparseMatrix kin_netProductionRates_ddX(CxxKinetics*) except +translate_exception

    # Transport properties
    cdef void tran_getMixDiffCoeffs(CxxTransport*, double*) except +translate_exception
    cdef void tran_getMixDiffCoeffsMass(CxxTransport*, double*) except +translate_exception
    cdef void tran_getMixDiffCoeffsMole(CxxTransport*, double*) except +translate_exception
    cdef void tran_getThermalDiffCoeffs(CxxTransport*, double*) except +translate_exception
    cdef void tran_getSpeciesViscosities(CxxTransport*, double*) except +translate_exception
    cdef void tran_getMobilities(CxxTransport*, double*) except +translate_exception

    cdef void tran_getMultiDiffCoeffs(CxxTransport*, size_t, double*) except +translate_exception
    cdef void tran_getBinaryDiffCoeffs(CxxTransport*, size_t, double*) except +translate_exception

    cdef void tran_getViscosityPolynomial(CxxTransport*, size_t, double*) except +translate_exception
    cdef void tran_getConductivityPolynomial(CxxTransport*, size_t, double*) except +translate_exception
    cdef void tran_getBinDiffusivityPolynomial(CxxTransport*, size_t, size_t, double*) except +translate_exception

    cdef void tran_setViscosityPolynomial(CxxTransport*, size_t, double*) except +translate_exception
    cdef void tran_setConductivityPolynomial(CxxTransport*, size_t, double*) except +translate_exception
    cdef void tran_setBinDiffusivityPolynomial(CxxTransport*, size_t, size_t, double*) except +translate_exception

# typedefs
ctypedef void (*thermoMethod1d)(CxxThermoPhase*, double*) except +translate_exception
ctypedef void (*transportMethod1d)(CxxTransport*, double*) except +translate_exception
ctypedef void (*transportMethod2d)(CxxTransport*, size_t, double*) except +translate_exception
ctypedef void (*transportPolyMethod1i)(CxxTransport*, size_t, double*) except +translate_exception
ctypedef void (*transportPolyMethod2i)(CxxTransport*, size_t, size_t, double*) except +translate_exception
ctypedef void (*kineticsMethod1d)(CxxKinetics*, double*) except +translate_exception
ctypedef CxxSparseMatrix (*kineticsMethodSparse)(CxxKinetics*) except +translate_exception

# classes
cdef class Units:
    cdef CxxUnits units
    @staticmethod
    cdef copy(CxxUnits)

cdef class UnitSystem:
    cdef CxxUnitSystem unitsystem

cdef class SpeciesThermo:
    cdef shared_ptr[CxxSpeciesThermo] _spthermo
    cdef CxxSpeciesThermo* spthermo
    cdef _assign(self, shared_ptr[CxxSpeciesThermo] other)

cdef class GasTransportData:
    cdef shared_ptr[CxxTransportData] _data
    cdef CxxGasTransportData* data
    cdef _assign(self, shared_ptr[CxxTransportData] other)

cdef class _SolutionBase:
    cdef shared_ptr[CxxSolution] _base
    cdef CxxSolution* base
    cdef CxxThermoPhase* thermo
    cdef CxxKinetics* kinetics
    cdef CxxTransport* transport
    cdef int thermo_basis
    cdef np.ndarray _selected_species
    cdef object parent
    cdef object _adjacent
    cdef public object _references

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

cdef class ReactionRate:
    cdef shared_ptr[CxxReactionRate] _rate
    cdef CxxReactionRate* rate
    @staticmethod
    cdef wrap(shared_ptr[CxxReactionRate])
    cdef set_cxx_object(self)

cdef class ArrheniusRateBase(ReactionRate):
    cdef CxxArrheniusBase* base

cdef class FalloffRate(ReactionRate):
    cdef CxxFalloffRate* falloff
    cdef set_cxx_object(self)

cdef class CustomRate(ReactionRate):
    cdef CxxCustomFunc1Rate* cxx_object(self)
    cdef Func1 _rate_func

cdef class InterfaceRateBase(ArrheniusRateBase):
    cdef CxxInterfaceRateBase* interface

cdef class StickRateBase(InterfaceRateBase):
    cdef CxxStickingCoverage* stick

cdef class Reaction:
    cdef shared_ptr[CxxReaction] _reaction
    cdef CxxReaction* reaction
    @staticmethod
    cdef wrap(shared_ptr[CxxReaction])

cdef class CustomReaction(Reaction):
    cdef CustomRate _rate

cdef class Arrhenius:
    cdef CxxArrhenius2* legacy # used by legacy objects only
    cdef CxxArrheniusRate* base
    cdef cbool own_rate
    cdef Reaction reaction # parent reaction, to prevent garbage collection
    @staticmethod
    cdef wrap(CxxArrheniusRate*)

cdef class Falloff:
    cdef shared_ptr[CxxFalloff] _falloff
    cdef CxxFalloff* falloff

cdef class Kinetics(_SolutionBase):
    pass

cdef class InterfaceKinetics(Kinetics):
    pass

cdef class Transport(_SolutionBase):
     pass

cdef class DustyGasTransport(Transport):
     pass

cdef class YamlWriter:
    cdef shared_ptr[CxxYamlWriter] _writer
    cdef CxxYamlWriter* writer
    @staticmethod
    cdef CxxUnitSystem _get_unitsystem(UnitSystem units)

cdef class Mixture:
    cdef CxxMultiPhase* mix
    cdef list _phases
    cdef object _weakref_proxy
    cpdef int element_index(self, element) except *

cdef class Func1:
    cdef shared_ptr[CxxFunc1] _func
    cdef CxxFunc1* func
    cdef object callable
    cdef object exception
    cpdef void _set_callback(self, object) except *

cdef class TabulatedFunction(Func1):
    cpdef void _set_tables(self, object, object, string) except *

cdef class ReactorBase:
    cdef CxxReactorBase* rbase
    cdef object _thermo
    cdef list _inlets
    cdef list _outlets
    cdef list _walls
    cdef list _surfaces
    cdef object _weakref_proxy

cdef class Reactor(ReactorBase):
    cdef CxxReactor* reactor
    cdef object _kinetics

cdef class Reservoir(ReactorBase):
    pass

cdef class ConstPressureReactor(Reactor):
    pass

cdef class IdealGasReactor(Reactor):
    pass

cdef class IdealGasConstPressureReactor(Reactor):
    pass

cdef class FlowReactor(Reactor):
    pass

cdef class ExtensibleReactor(Reactor):
    cdef CxxReactorAccessor* accessor

cdef class ReactorSurface:
    cdef CxxReactorSurface* surface
    cdef Kinetics _kinetics

cdef class WallBase:
    cdef CxxWallBase* wall
    cdef object _velocity_func
    cdef object _heat_flux_func
    cdef ReactorBase _left_reactor
    cdef ReactorBase _right_reactor
    cdef str name

cdef class Wall(WallBase):
    pass

cdef class FlowDevice:
    cdef CxxFlowDevice* dev
    cdef Func1 _rate_func
    cdef Func1 _time_func
    cdef str name
    cdef ReactorBase _upstream
    cdef ReactorBase _downstream

cdef class MassFlowController(FlowDevice):
    pass

cdef class Valve(FlowDevice):
    pass

cdef class PressureController(FlowDevice):
    pass

cdef class ReactorNet:
    cdef CxxReactorNet net
    cdef list _reactors

cdef class Domain1D:
    cdef CxxDomain1D* domain
    cdef _SolutionBase gas
    cdef object _weakref_proxy
    cdef public pybool have_user_tolerances

cdef class Boundary1D(Domain1D):
    cdef CxxBoundary1D* boundary

cdef class Inlet1D(Boundary1D):
    cdef CxxInlet1D* inlet

cdef class Outlet1D(Boundary1D):
    cdef CxxOutlet1D* outlet

cdef class OutletReservoir1D(Boundary1D):
    cdef CxxOutletRes1D* outlet

cdef class SymmetryPlane1D(Boundary1D):
    cdef CxxSymm1D* symm

cdef class Surface1D(Boundary1D):
    cdef CxxSurf1D* surf

cdef class ReactingSurface1D(Boundary1D):
    cdef CxxReactingSurf1D* surf
    cdef public Kinetics surface

cdef class _FlowBase(Domain1D):
    cdef CxxStFlow* flow

cdef class IdealGasFlow(_FlowBase):
    pass

cdef class FreeFlow(IdealGasFlow):
    pass

cdef class IonFlow(_FlowBase):
    pass

cdef class AxisymmetricStagnationFlow(IdealGasFlow):
    pass

cdef class Sim1D:
    cdef CxxSim1D* sim
    cdef readonly object domains
    cdef object _initialized
    cdef object _initial_guess_args
    cdef object _initial_guess_kwargs
    cdef public Func1 _interrupt
    cdef public Func1 _time_step_callback
    cdef public Func1 _steady_callback

cdef class ReactionPathDiagram:
    cdef CxxReactionPathDiagram diagram
    cdef CxxReactionPathBuilder builder
    cdef Kinetics kinetics
    cdef str element
    cdef pybool built
    cdef CxxStringStream* _log

# free functions
cdef string stringify(x) except *
cdef pystr(string x)
cdef np.ndarray get_species_array(Kinetics kin, kineticsMethod1d method)
cdef np.ndarray get_reaction_array(Kinetics kin, kineticsMethod1d method)
cdef np.ndarray get_dense(Kinetics kin, kineticsMethodSparse method)
cdef tuple get_sparse(Kinetics kin, kineticsMethodSparse method)
cdef np.ndarray get_transport_1d(Transport tran, transportMethod1d method)
cdef np.ndarray get_transport_2d(Transport tran, transportMethod2d method)
cdef np.ndarray get_transport_polynomial(Transport tran, transportPolyMethod1i method, int index, int n_coeffs)
cdef np.ndarray get_binary_transport_polynomial(Transport tran, transportPolyMethod2i method, int indexi, int indexj, int n_coeffs)
cdef CxxIdealGasPhase* getIdealGasPhase(ThermoPhase phase) except *
cdef wrapSpeciesThermo(shared_ptr[CxxSpeciesThermo] spthermo)
cdef int assign_delegates(object, CxxDelegator*) except -1

cdef extern from "cantera/thermo/Elements.h" namespace "Cantera":
    double getElementWeight(string ename) except +translate_exception
    double getElementWeight(int atomicNumber) except +translate_exception
    int numElementsDefined()
    int getAtomicNumber(string ename) except +translate_exception
    string getElementSymbol(string ename) except +translate_exception
    string getElementSymbol(int atomicNumber) except +translate_exception
    string getElementName(string ename) except +translate_exception
    string getElementName(int atomicNumber) except +translate_exception
