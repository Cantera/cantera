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

    virtual void setParameters(
        const AnyMap& node, const UnitStack& rate_units) override
    {
        RateType::m_negativeA_ok = node.getBool("negative-A", false);
        if (!node.hasKey("rate-constant")) {
            RateType::setRateParameters(AnyValue(), node.units(), rate_units);
            return;
        }
        RateType::setRateParameters(node["rate-constant"], node.units(), rate_units);
    }

    virtual void getParameters(AnyMap& node) const override {
        if (RateType::m_negativeA_ok) {
            node["negative-A"] = true;
        }
        AnyMap rateNode;
        RateType::getRateParameters(rateNode);
        if (!rateNode.empty()) {
            // RateType object is configured
            node["rate-constant"] = std::move(rateNode);
        }
        if (RateType::type() != "Arrhenius") {
            node["type"] = RateType::type();
        }
    }
};

typedef BulkRate<Arrhenius3, ArrheniusData> ArrheniusRate;
typedef BulkRate<TwoTempPlasma, TwoTempPlasmaData> TwoTempPlasmaRate;
typedef BulkRate<BlowersMasel, BlowersMaselData> BlowersMaselRate;

}

#endif
