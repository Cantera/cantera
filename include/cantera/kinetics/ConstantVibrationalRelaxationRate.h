//! @file ConstantVibrationalRelaxationRate.h
//! Constant vibrational relaxation reaction rate.
//! @since New in %Cantera 4.0

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CONSTANTVIBRATIONALRELAXATIONRATE_H
#define CT_CONSTANTVIBRATIONALRELAXATIONRATE_H

#include "cantera/kinetics/MultiRate.h"
#include "cantera/kinetics/VibrationalRelaxationRate.h"

namespace Cantera
{

//! Constant vibrational relaxation rate implementation.
/**
 * Maps the user-facing constant rate coefficient onto the common
 * VibrationalRelaxationRate representation with all temperature-dependent
 * terms set to zero.
 * 
 * Implements
 *
 * @f[
 * k_f = A.
 * @f]
 * 
 * @ingroup otherRateGroup
 */
class ConstantVibrationalRelaxationRate final
    : public VibrationalRelaxationRate
{
public:
    ConstantVibrationalRelaxationRate() = default;

    explicit ConstantVibrationalRelaxationRate(double A);

    explicit ConstantVibrationalRelaxationRate(
        const AnyMap& node, const UnitStack& rate_units = {});

    void setParameters(
        const AnyMap& node, const UnitStack& rate_units) override;

    void getParameters(AnyMap& node) const override;

    unique_ptr<MultiRateBase> newMultiRate() const override
    {
        return make_unique<
            MultiRate<ConstantVibrationalRelaxationRate,
                      VibrationalRelaxationData>>();
    }

    const string type() const override
    {
        return "constant-vibrational-relaxation";
    }
};

} // namespace Cantera

#endif