/**
 * @file BulkRate.h
 * Header for reaction rates that occur in bulk phases.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_BULKRATE_H
#define CT_BULKRATE_H

#include "Arrhenius.h"

namespace Cantera
{

//! A class template for bulk phase reaction rate specifications
template <class RateType, class DataType>
class BulkRate : public RateType
{
public:
    BulkRate() = default;
    using RateType::RateType; // inherit constructors

    //! Constructor based on AnyMap content
    BulkRate(const AnyMap& node, const UnitStack& rate_units={}) {
        setParameters(node, rate_units);
    }

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return unique_ptr<MultiRateBase>(
            new MultiRate<BulkRate<RateType, DataType>, DataType>);
    }

    using RateType::setParameters;
    using RateType::getParameters;
};


using ArrheniusRate = BulkRate<Arrhenius3, ReactionData>;
using TwoTempPlasmaRate = BulkRate<TwoTempPlasma, TwoTempPlasmaData>;
using BlowersMaselRate = BulkRate<BlowersMasel, BulkData>;


class ThreeBodyBase
{
public:
    ThreeBodyBase();

    //! Perform object setup based on AnyMap node information
    //! @param node  AnyMap object containing reaction rate specification
    void setParameters(const AnyMap& node);

    //! Store parameters needed to reconstruct an identical object
    //! @param node  AnyMap object receiving reaction rate specification
    void getParameters(AnyMap& node) const;

    //! Get third-body collision efficiency parameters
    //! @param efficiencies  AnyMap receiving coverage information
    void getEfficiencies(AnyMap& efficiencies) const;

    //! Get the default efficiency
    double defaultEfficiency() const {
        return m_defaultEfficiency;
    }

    //! Get the third-body efficiency for species *k*
    double efficiency(const std::string& k) const;

    //! Build rate-specific parameters based on Reaction and Kinetics context
    //! @param rxn  Associated reaction object
    //! @param kin  Kinetics object holding the rate evaluator
    void setContext(const Reaction& rxn, const Kinetics& kin);

    //! Update reaction rate parameters
    //! @param shared_data  data shared by all reactions of a given type
    void updateFromStruct(const BulkData& shared_data) {
        if (shared_data.ready) {
            m_thirdBodyConc = m_defaultEfficiency * shared_data.molarDensity;
            for (const auto& eff : m_efficiencies) {
                m_thirdBodyConc += shared_data.concentrations[eff.first] * eff.second;
            }
        }
    }

    //! Third-body concentration
    double thirdBodyConcentration() const {
        return m_thirdBodyConc;
    }

    //! Apply correction
    double applyCorrection(double value) const {
        if (m_massAction) {
            return value * m_thirdBodyConc;
        }
        return value;
    }

protected:
    double m_thirdBodyConc; //!< Effective third-body concentration

    //! The default third body efficiency for species not listed in
    double m_defaultEfficiency;

    //! Input explicitly specifies collision partner
    bool m_specifiedCollisionPartner;

    //! Third body is used by law of mass action
    //! (`true` for three-body reactions, `false` for falloff reactions)
    bool m_massAction;

    Composition m_efficiencyMap; //!< Map of species to third body efficiency

private:
    //! Vector of pairs containing indices and efficiencies
    std::vector<std::pair<size_t, double>> m_efficiencies;
};


template <class RateType, class DataType>
class ThreeBodyRate : public RateType, public ThreeBodyBase
{
    CT_DEFINE_HAS_MEMBER(has_update, updateFromStruct)

public:
    ThreeBodyRate() = default;
    using RateType::RateType; // inherit constructors

    //! Constructor based on AnyMap content
    ThreeBodyRate(const AnyMap& node, const UnitStack& rate_units={}) {
        setParameters(node, rate_units);
    }

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return unique_ptr<MultiRateBase>(
            new MultiRate<ThreeBodyRate<RateType, DataType>, DataType>);
    }

    virtual const std::string type() const override {
        return "three-body-" + RateType::type();
    }

    virtual void setParameters(
        const AnyMap& node, const UnitStack& rate_units) override
    {
        RateType::setParameters(node, rate_units);
        ThreeBodyBase::setParameters(node);
    }

    virtual void getParameters(AnyMap& node) const override {
        RateType::getParameters(node);
        node["type"] = type();
        ThreeBodyBase::getParameters(node);
    }

    virtual void setContext(const Reaction& rxn, const Kinetics& kin) override {
        RateType::setContext(rxn, kin);
        ThreeBodyBase::setContext(rxn, kin);
    }

    //! Update reaction rate parameters
    //! @param shared_data  data shared by all reactions of a given type
    void updateFromStruct(const DataType& shared_data) {
        _update(shared_data);
        ThreeBodyBase::updateFromStruct(shared_data);
    }

    //! Evaluate reaction rate
    //! @param shared_data  data shared by all reactions of a given type
    double evalFromStruct(const DataType& shared_data) const {
        return applyCorrection(RateType::evalFromStruct(shared_data));
    }

protected:
    //! Helper function to process updates for rate types that implement the
    //! `updateFromStruct` method.
    template <typename T=RateType,
        typename std::enable_if<has_update<T>::value, bool>::type = true>
    void _update(const DataType& shared_data) {
        T::updateFromStruct(shared_data);
    }

    //! Helper function for rate types that do not implement `updateFromStruct`.
    //! Does nothing, but exists to allow generic implementations of update().
    template <typename T=RateType,
        typename std::enable_if<!has_update<T>::value, bool>::type = true>
    void _update(const DataType& shared_data) {
    }
};

using ThreeBodyArrheniusRate = ThreeBodyRate<Arrhenius3, BulkData>;
using ThreeBodyBlowersMaselRate = ThreeBodyRate<BlowersMasel, BulkData>;

}

#endif
