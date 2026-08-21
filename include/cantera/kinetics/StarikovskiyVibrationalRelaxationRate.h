//! @file StarikovskiyVibrationalRelaxationRate.h
//! Starikovskiy vibrational relaxation reaction rate.
//! @since New in %Cantera 4.0

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_STARIKOVSKIYVIBRATIONALRELAXATIONRATE_H
#define CT_STARIKOVSKIYVIBRATIONALRELAXATIONRATE_H

#include "cantera/kinetics/MultiRate.h"
#include "cantera/kinetics/VibrationalRelaxationRate.h"

namespace Cantera
{

//! Starikovskiy vibrational relaxation rate.
class StarikovskiyVibrationalRelaxationRate final
    : public VibrationalRelaxationRate
{
public:
    StarikovskiyVibrationalRelaxationRate() = default;

    explicit StarikovskiyVibrationalRelaxationRate(
        const AnyMap& node, const UnitStack& rate_units = {});

    void setParameters(
        const AnyMap& node, const UnitStack& rate_units) override;

    void getParameters(AnyMap& node) const override;

    unique_ptr<MultiRateBase> newMultiRate() const override
    {
        return make_unique<
            MultiRate<StarikovskiyVibrationalRelaxationRate,
                      VibrationalRelaxationData>>();
    }

    const string type() const override
    {
        return "Starikovskiy-vibrational-relaxation";
    }
};

} // namespace Cantera

#endif