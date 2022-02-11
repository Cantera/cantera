/**
 * @file Coverage.h
 * Header for reaction rates that occur at interfaces.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_COVERAGE_H
#define CT_COVERAGE_H

#include "cantera/kinetics/Arrhenius.h"
#include "MultiRate.h"

namespace Cantera
{

class AnyMap;

/**
 *  @defgroup surfaceGroup  Coverage-dependent rate parameterizations
 *
 *  This section describes the parameterizations used to describe rate
 *  parameterization that involve interfaces.
 *
 *  @ingroup chemkinetics
 */

//! Base class for reaction rate parameterizations that involve interfaces
class Coverage
{
public:
    Coverage();

    //! Perform object setup based on AnyMap node information
    /*!
     *  @param dependencies  Coverage dependencies
     *  @param units  Unit system
     */
    void setCoverageDependencies(const AnyMap& dependencies,
                                 const UnitSystem& units=UnitSystem());

    void getCoverageDependencies(AnyMap& dependencies, bool asVector=false) const;

    //! Add a coverage dependency for species *sp*, with exponential dependence
    //! *a*, power-law exponent *m*, and activation energy dependence *e*,
    //! where *e* is in Kelvin, i.e. energy divided by the molar gas constant.
    void addCoverageDependence(std::string sp, double a, double m, double e);

    //! Set species indices within coverage array
    void setSpecies(const Kinetics& kin);
    void setSpecies(const std::vector<std::string>& species);

    //! Update reaction rate parameters
    //! @param shared_data  data shared by all reactions of a given type
    void updateFromStruct(const CoverageData& shared_data) {
        if (shared_data.ready) {
            m_siteDensity = shared_data.siteDensity;
        }
        if (m_indices.size() != m_cov.size()) {
            // object is not set up correctly (setSpecies needs to be run)
            m_acov = NAN;
            m_ecov = NAN;
            m_mcov = NAN;
            return;
        }
        m_acov = 0.0;
        m_ecov = 0.0;
        m_mcov = 0.0;
        for (auto& item : m_indices) {
            m_acov += m_ac[item.first] * shared_data.coverages[item.second];
            m_ecov += m_ec[item.first] * shared_data.coverages[item.second];
            m_mcov += m_mc[item.first] * shared_data.logCoverages[item.second];
        }
    }

    //! Return site density [kmol/m^2]
    double siteDensity() const {
        return m_siteDensity;
    }

    //! Set site density [kmol/m^2]
    //! @internal  used for testing purposes only
    //! Note that this quantity is not an independent variable and will be
    //! overwritten during an update of the thermodynamic state.
    void setSiteDensity(double siteDensity) {
        m_siteDensity = siteDensity;
    }

protected:
    double m_siteDensity; //!< Site density [kmol/m^2]
    double m_acov; //!< Coverage contribution to pre-exponential factor
    double m_ecov; //!< Coverage contribution to activation energy
    double m_mcov; //!< Coverage term in reaction rate
    std::map<size_t, size_t> m_indices; //!< Map holding indices of coverage species
    std::vector<std::string> m_cov; //!< Vector holding names of coverage species
    vector_fp m_ac; //!< Vector holding coverage-specific exponential dependence
    vector_fp m_ec; //!< Vector holding coverage-specific activation energy dependence
    vector_fp m_mc; //!< Vector holding coverage-specific power-law exponents
};


//! A class template for interface reaction rate specifications
template <class RateType, class DataType>
class InterfaceRate : public RateType, public Coverage
{
public:
    InterfaceRate() = default;
    using RateType::RateType; // inherit constructors

    //! Constructor based on AnyMap content
    InterfaceRate(const AnyMap& node, const UnitStack& rate_units={}) {
        setParameters(node, rate_units);
    }

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return unique_ptr<MultiRateBase>(
            new MultiRate<InterfaceRate<RateType, DataType>, DataType>);
    }

    //! Identifier of reaction rate type
    virtual const std::string type() const override {
        return RateType::type() + "-interface";
    }

    virtual void setParameters(
        const AnyMap& node, const UnitStack& rate_units) override
    {
        if (node.hasKey("coverage-dependencies")) {
            setCoverageDependencies(
                node["coverage-dependencies"].as<AnyMap>(), node.units());
        }
        RateType::m_negativeA_ok = node.getBool("negative-A", false);
        if (!node.hasKey("rate-constant")) {
            RateType::setRateParameters(AnyValue(), node.units(), rate_units);
            return;
        }
        RateType::setRateParameters(node["rate-constant"], node.units(), rate_units);
    }

    virtual void getParameters(AnyMap& node) const override {
        node["type"] = type();
        if (RateType::m_negativeA_ok) {
            node["negative-A"] = true;
        }
        AnyMap rateNode;
        RateType::getRateParameters(rateNode);
        if (!rateNode.empty()) {
            // RateType object is configured
            node["rate-constant"] = std::move(rateNode);
        }
        if (!m_cov.empty()) {
            AnyMap deps;
            getCoverageDependencies(deps);
            node["coverage-dependencies"] = std::move(deps);
        }
    }

    virtual void setContext(const Reaction& rxn, const Kinetics& kin) override {
        RateType::setContext(rxn, kin);
        setSpecies(kin);
    }

    //! Evaluate reaction rate
    //! @param shared_data  data shared by all reactions of a given type
    double evalFromStruct(const DataType& shared_data) const {
        return RateType::evalRate(shared_data.logT, shared_data.recipT) *
            std::exp(std::log(10.0) * m_acov - m_ecov * shared_data.recipT + m_mcov);
    }

    //! Evaluate derivative of reaction rate with respect to temperature
    //! divided by reaction rate
    //! @param shared_data  data shared by all reactions of a given type
    double ddTScaledFromStruct(const DataType& shared_data) const {
        throw NotImplementedError("InterfaceRate<>::ddTScaledFromStruct");
    }

    //! Return the effective pre-exponential factor *A* (in m, kmol, s to powers
    //! depending on the reaction order) accounting for coverage dependent terms.
    double effectivePreExponentialFactor() const {
        return RateType::preExponentialFactor() *
            std::exp(std::log(10.0) * m_acov + m_mcov);
    }

    //! Return the effective activation energy *Ea* [J/kmol]
    double effectiveActivationEnergy() const {
        return RateType::activationEnergy() + m_ecov * GasConstant;
    }

    //! Return the effective activation temperature (energy divided by the gas
    //! constant) [K]
    double effectiveActivationEnergy_R() const {
        return RateType::activationEnergy_R() + m_ecov;
    }
};

using ArrheniusInterfaceRate = InterfaceRate<Arrhenius3, CoverageData>;

}
#endif
