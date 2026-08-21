//! @file CastelaVibrationalRelaxationRate.h
//! Castela vibrational relaxation reaction rate.
//! @since New in %Cantera 4.0

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CASTELAVIBRATIONALRELAXATIONRATE_H
#define CT_CASTELAVIBRATIONALRELAXATIONRATE_H

#include "cantera/kinetics/MultiRate.h"
#include "cantera/kinetics/VibrationalRelaxationRate.h"

namespace Cantera
{

//! Castela vibrational relaxation rate.
class CastelaVibrationalRelaxationRate final
    : public VibrationalRelaxationRate
{
public:
    CastelaVibrationalRelaxationRate() = default;

    explicit CastelaVibrationalRelaxationRate(
        const AnyMap& node, const UnitStack& rate_units = {});

    void setParameters(
        const AnyMap& node, const UnitStack& rate_units) override;

    void getParameters(AnyMap& node) const override;

    unique_ptr<MultiRateBase> newMultiRate() const override
    {
        return make_unique<
            MultiRate<CastelaVibrationalRelaxationRate,
                      VibrationalRelaxationData>>();
    }

    const string type() const override
    {
        return "Castela-vibrational-relaxation";
    }

private:
    //! Castela coefficient a.
    double m_castela_a = 0.0;

    //! Castela coefficient b.
    double m_castela_b = 0.0;

    //! Reference pressure.
    double m_referencePressure = OneAtm;
};

} // namespace Cantera

#endif