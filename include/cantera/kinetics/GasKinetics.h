/**
 * @file GasKinetics.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_GASKINETICS_H
#define CT_GASKINETICS_H

#include "BulkKinetics.h"
#include "Reaction.h"

#ifndef CT_SKIP_DEPRECATION_WARNINGS
#pragma message("warning: GasKinetics.h and class GasKinetics are deprecated and will " \
                "be removed after Cantera 3.0. Replace with class BulkKinetics.")
#endif

namespace Cantera
{

/**
 * Kinetics manager for elementary gas-phase chemistry. This kinetics manager
 * implements standard mass-action reaction rate expressions for low-density
 * gases.
 * @ingroup kineticsmgr
 * @deprecated Replace with BulkKinetics. To be removed after %Cantera 3.0.
 */
class GasKinetics : public BulkKinetics
{
public:

    //! Constructor.
    GasKinetics() {}

    //! @deprecated To be removed after %Cantera 3.0; code base only uses default.
    GasKinetics(ThermoPhase* thermo) : GasKinetics() {
        warn_deprecated("GasKinetics::GasKinetics(ThermoPhase*)",
            "To be removed after Cantera 3.0. Use default constructor instead.");
        addPhase(*thermo);
    }

    string kineticsType() const override {
        return "gas";
    }

    //! Update temperature-dependent portions of reaction rates and falloff
    //! functions.
    virtual void update_rates_T() {
        warn_deprecated("GasKinetics::update_rates_T",
            "Class GasKinetics has been merged into class BulkKinetics, and the "
            "update_rates_T() method is now part of updateROP(). Class GasKinetics "
            "will be removed after Cantera 3.0.");
        updateROP();
    }

    //! Update properties that depend on concentrations.
    //! Currently the enhanced collision partner concentrations are updated
    //! here, as well as the pressure-dependent portion of P-log and Chebyshev
    //! reactions.
    virtual void update_rates_C() {
        warn_deprecated("GasKinetics::update_rates_C",
            "Class GasKinetics has been merged into class BulkKinetics, and the "
            "update_rates_T() method is now part of updateROP(). Class GasKinetics "
            "will be removed after Cantera 3.0.");
        updateROP();
    }

};

}

#endif
