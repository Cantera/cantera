# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

#cython: language_level=3
#distutils: language = c++

from .ctcxx cimport *
from .func1 cimport *
from .delegator cimport CxxDelegator

# @todo Replace with import from libcpp.map once Cython 0.29.x is no longer supported
cdef extern from "<map>" namespace "std" nogil:
    cdef cppclass multimap[T, U]:
        cppclass iterator:
            pair[T, U]& operator*()
            iterator operator++()
            iterator operator--()
            bint operator==(iterator)
            bint operator!=(iterator)
        multimap() except +
        U& operator[](T&)
        iterator begin()
        iterator end()
        pair[iterator, bint] insert(pair[T, U])
        iterator find(T&)


cdef extern from "cantera/kinetics/ReactionRateFactory.h" namespace "Cantera":
    cdef shared_ptr[CxxReactionRate] CxxNewReactionRate "newReactionRate" (CxxAnyMap&) except +translate_exception


cdef extern from "cantera/kinetics/ReactionRate.h" namespace "Cantera":
    cdef cppclass CxxReactionRate "Cantera::ReactionRate":
        CxxReactionRate()
        string type()
        string subType()
        double eval(double) except +translate_exception
        double eval(double, double) except +translate_exception
        double eval(double, vector[double]&) except +translate_exception
        CxxAnyMap parameters() except +translate_exception
        CxxUnits conversionUnits()


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

cdef extern from "cantera/kinetics/ElectronCollisionPlasmaRate.h" namespace "Cantera":
    cdef cppclass CxxElectronCollisionPlasmaRate "Cantera::ElectronCollisionPlasmaRate" (CxxReactionRate):
        CxxElectronCollisionPlasmaRate(CxxAnyMap) except +translate_exception
        vector[double]& energyLevels()
        vector[double]& crossSections()

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
    cdef shared_ptr[CxxReaction] CxxNewReaction "newReaction" (CxxAnyMap&, CxxKinetics&) except +translate_exception
    cdef vector[shared_ptr[CxxReaction]] CxxGetReactions "getReactions" (CxxAnyValue&, CxxKinetics&) except +translate_exception

    cdef cppclass CxxThirdBody "Cantera::ThirdBody":
        CxxThirdBody()
        CxxThirdBody(string)
        double efficiency(string)
        string name()
        Composition efficiencies
        double default_efficiency
        cbool mass_action

    cdef cppclass CxxReaction "Cantera::Reaction":
        CxxReaction()
        CxxReaction(Composition&, Composition&, shared_ptr[CxxReactionRate]) except +translate_exception
        CxxReaction(Composition&, Composition&, shared_ptr[CxxReactionRate], shared_ptr[CxxThirdBody]) except +translate_exception
        CxxReaction(string&, shared_ptr[CxxReactionRate]) except +translate_exception
        CxxReaction(string&, shared_ptr[CxxReactionRate], shared_ptr[CxxThirdBody]) except +translate_exception

        string reactantString()
        string productString()
        string equation()
        void setEquation(const string&) except +translate_exception
        string type()
        CxxAnyMap parameters(cbool) except +translate_exception
        cbool usesThirdBody()
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


cdef extern from "cantera/kinetics/Falloff.h" namespace "Cantera":
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


cdef extern from "cantera/kinetics/PlogRate.h" namespace "Cantera":
    cdef cppclass CxxPlogRate "Cantera::PlogRate" (CxxReactionRate):
        CxxPlogRate()
        CxxPlogRate(CxxAnyMap) except +translate_exception
        CxxPlogRate(multimap[double, CxxArrheniusRate])
        multimap[double, CxxArrheniusRate] getRates()

cdef extern from "cantera/kinetics/LmrRate.h" namespace "Cantera":
    cdef cppclass CxxLmrRate "Cantera::LmrRate" (CxxReactionRate):
        CxxLmrRate()
        CxxLmrRate(CxxAnyMap) except +translate_exception

cdef extern from "cantera/kinetics/ChebyshevRate.h" namespace "Cantera":
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


cdef extern from "cantera/kinetics/Custom.h" namespace "Cantera":
    cdef cppclass CxxCustomFunc1Rate "Cantera::CustomFunc1Rate" (CxxReactionRate):
        CxxCustomFunc1Rate()
        void setRateFunction(shared_ptr[CxxFunc1]) except +translate_exception


cdef extern from "cantera/kinetics/ReactionRateDelegator.h" namespace "Cantera":
    cdef cppclass CxxReactionDataDelegator "Cantera::ReactionDataDelegator" (CxxDelegator):
        CxxReactionDataDelegator()

    cdef cppclass CxxReactionRateDelegator "Cantera::ReactionRateDelegator" (CxxDelegator, CxxReactionRate):
        CxxReactionRateDelegator()
        void setType(string&)


cdef extern from "cantera/kinetics/InterfaceRate.h" namespace "Cantera":
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


cdef class ReactionRate:
    cdef shared_ptr[CxxReactionRate] _rate
    cdef CxxReactionRate* rate
    @staticmethod
    cdef wrap(shared_ptr[CxxReactionRate])
    cdef set_cxx_object(self)

cdef class ArrheniusRateBase(ReactionRate):
    cdef CxxArrheniusBase* base

cdef class ElectronCollisionPlasmaRate(ReactionRate):
    cdef CxxElectronCollisionPlasmaRate* base
    cdef set_cxx_object(self)

cdef class FalloffRate(ReactionRate):
    cdef CxxFalloffRate* falloff
    cdef set_cxx_object(self)

cdef class CustomRate(ReactionRate):
    cdef CxxCustomFunc1Rate* cxx_object(self)
    cdef Func1 _rate_func  # prevent premature garbage collection

cdef class ExtensibleRate(ReactionRate):
    cdef public list _delegates
    cdef set_cxx_object(self, CxxReactionRate* rate=*)

cdef class ExtensibleRateData:
    cdef public list _delegates
    cdef set_cxx_object(self, CxxReactionDataDelegator* rate)

cdef class InterfaceRateBase(ArrheniusRateBase):
    cdef CxxInterfaceRateBase* interface

cdef class StickRateBase(InterfaceRateBase):
    cdef CxxStickingCoverage* stick

cdef class ThirdBody:
    cdef shared_ptr[CxxThirdBody] _third_body
    cdef CxxThirdBody* third_body
    @staticmethod
    cdef wrap(shared_ptr[CxxThirdBody])

cdef class Reaction:
    cdef shared_ptr[CxxReaction] _reaction
    cdef CxxReaction* reaction
    @staticmethod
    cdef wrap(shared_ptr[CxxReaction])
    cdef ReactionRate _rate

cdef class Arrhenius:
    cdef CxxArrheniusRate* base
    cdef cbool own_rate
    cdef Reaction reaction # parent reaction, to prevent garbage collection
    @staticmethod
    cdef wrap(CxxArrheniusRate*)
