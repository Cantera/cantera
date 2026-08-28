//! @file MultiStateResolvedVibrationalRelaxationRate.h
//! Multi-state-resolved vibrational relaxation reaction rate.
//! @since New in %Cantera 4.0

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MULTISTATERESOLVEDVIBRATIONALRELAXATIONRATE_H
#define CT_MULTISTATERESOLVEDVIBRATIONALRELAXATIONRATE_H

#include "cantera/kinetics/MultiRate.h"
#include "cantera/kinetics/VibrationalRelaxationRate.h"

namespace Cantera
{

//! Multi-state-resolved vibrational relaxation rate implementation.
/**
 * Maps the coefficients A, b, C0, C13, and C23 onto the common
 * VibrationalRelaxationRate representation, with m fixed to 2/3.
 *
 * @ingroup otherRateGroup
 */
class MultiStateResolvedVibrationalRelaxationRate final
    : public VibrationalRelaxationRate
{
public:
    MultiStateResolvedVibrationalRelaxationRate() = default;

    MultiStateResolvedVibrationalRelaxationRate(double A, double b,
        double C0 = 0.0, double C13 = 0.0, double C23 = 0.0);

    explicit MultiStateResolvedVibrationalRelaxationRate(
        const AnyMap& node, const UnitStack& rate_units = {});

    void setParameters(
        const AnyMap& node, const UnitStack& rate_units) override;

    void getParameters(AnyMap& node) const override;

    unique_ptr<MultiRateBase> newMultiRate() const override
    {
        return make_unique<
            MultiRate<MultiStateResolvedVibrationalRelaxationRate,
                      VibrationalRelaxationData>>();
    }

    const string type() const override
    {
        return "multi-state-resolved-vibrational-relaxation";
    }
};

} // namespace Cantera

#endif