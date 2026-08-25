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
 * Maps the user-facing parameters n, K, B, C, m, D, and z onto the generic
 * coefficients stored by VibrationalRelaxationRate.
 *
 * @ingroup otherRateGroup
 */
class StarikovskiyVibrationalRelaxationRate final
    : public VibrationalRelaxationRate
{
public:
    StarikovskiyVibrationalRelaxationRate() = default;

    explicit StarikovskiyVibrationalRelaxationRate(double A, double n = 0.0, double K = 0.0,
        double B = 0.0, double C = 0.0, double m = 1.0, double D = 0.0,  double z = 1.0);

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