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

//! Starikovskiy vibrational relaxation rate implementation.
/**
 * Maps the user-facing parameters b, C0, C13, Cm, m, Cn, and n onto the generic
 * coefficients stored by VibrationalRelaxationRate.
 *
 * @ingroup otherRateGroup
 */
class StarikovskiyVibrationalRelaxationRate final
    : public VibrationalRelaxationRate
{
public:
    StarikovskiyVibrationalRelaxationRate() = default;

    explicit StarikovskiyVibrationalRelaxationRate(
        double A, double b = 0.0, double C0 = 0.0, double C13 = 0.0, double Cm = 0.0,
        double m = 1.0, double Cn = 0.0, double n = 1.0);

    explicit StarikovskiyVibrationalRelaxationRate(
        const AnyMap& node, const UnitStack& rate_units = {});

    void setParameters(
        const AnyMap& node, const UnitStack& rate_units) override;

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
    
protected:
    void getRateParameters(AnyMap& rateNode) const override;
};

} // namespace Cantera

#endif
