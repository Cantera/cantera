#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *
from .func1 cimport *

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


cdef extern from "cantera/base/Array.h" namespace "Cantera":
    cdef cppclass CxxArray2D "Cantera::Array2D":
        CxxArray2D()
        CxxArray2D(size_t, size_t)
        void resize(size_t, size_t)
        double operator()(size_t, size_t)
        vector[double]& data()
        size_t nRows()
        size_t nColumns()


cdef extern from "cantera/cython/wrappers.h":
    # workaround for Cython assignment limitations
    cdef void CxxArray2D_set(CxxArray2D, size_t, size_t, double)


cdef extern from "cantera/kinetics/Reaction.h" namespace "Cantera":
    cdef cppclass CxxKinetics "Cantera::Kinetics"
    cdef shared_ptr[CxxReaction] CxxNewReaction "Cantera::newReaction" (string) except +translate_exception
    cdef shared_ptr[CxxReaction] CxxNewReaction "newReaction" (CxxAnyMap&, CxxKinetics&) except +translate_exception
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
        CxxAnyMap parameters(cbool) except +translate_exception
        CxxAnyMap input
        Composition reactants
        Composition products
        Composition orders
        string id
        cbool reversible
        cbool duplicate
        cbool allow_nonreactant_orders
        cbool allow_negative_orders
        shared_ptr[CxxThirdBody] thirdBody()
        CxxUnits rate_units

        shared_ptr[CxxReactionRate] rate()
        void setRate(shared_ptr[CxxReactionRate])

    cdef cppclass CxxFalloff "Cantera::FalloffRate":
        CxxFalloff()
        void updateTemp(double, double*)
        double F(double, double*)
        size_t workSize()

        size_t nParameters()
        string type()
        void getParameters(double*)

    cdef cppclass CxxChebyshev "Cantera::ChebyshevRate":
        CxxChebyshev(double, double, double, double, CxxArray2D)
        double Tmin()
        double Tmax()
        double Pmin()
        double Pmax()
        size_t nTemperature()
        size_t nPressure()
        CxxArray2D& data()

    cdef cppclass CxxThreeBodyReaction "Cantera::ThreeBodyReaction" (CxxReaction):
        CxxThreeBodyReaction()

    cdef cppclass CxxFalloffReaction "Cantera::FalloffReaction" (CxxReaction):
        CxxFalloffReaction()

    cdef cppclass CxxCustomFunc1Reaction "Cantera::CustomFunc1Reaction" (CxxReaction):
        CxxCustomFunc1Reaction()


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
    cdef CxxArrheniusRate* base
    cdef cbool own_rate
    cdef Reaction reaction # parent reaction, to prevent garbage collection
    @staticmethod
    cdef wrap(CxxArrheniusRate*)

cdef class Falloff:
    cdef shared_ptr[CxxFalloff] _falloff
    cdef CxxFalloff* falloff
