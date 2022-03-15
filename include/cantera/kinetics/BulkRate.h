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

typedef BulkRate<Arrhenius3, ArrheniusData> ArrheniusRate;
typedef BulkRate<TwoTempPlasma, TwoTempPlasmaData> TwoTempPlasmaRate;
typedef BulkRate<BlowersMasel, BlowersMaselData> BlowersMaselRate;

}

#endif
