/**
 *  @file WaterTransport.h Header file defining class WaterTransport
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_WATERTRAN_H
#define CT_WATERTRAN_H

#include "cantera/thermo/WaterPropsIAPWS.h"
#include "cantera/transport/Transport.h"

namespace Cantera
{

//! Transport Parameters for pure water
//! @ingroup tranprops
class WaterTransport : public Transport
{
public:
    //! default constructor
    /*!
     *  @param thermo   ThermoPhase object that represents the phase.
     *                  Defaults to zero
     *  @param ndim     Number of dimensions of the flux expressions.
     *                  Defaults to a value of one.
     *
     * @deprecated The `thermo` and `ndim` parameters will be removed after %Cantera 3.0.
     *     The ThermoPhase object should be specifed when calling the `init` method.
     */
    WaterTransport(ThermoPhase* thermo = 0, int ndim = -1);

    string transportModel() const override {
        return "Water";
    }

    //! Returns the viscosity of water at the current conditions (kg/m/s)
    /*!
     * This function calculates the value of the viscosity of pure water at the
     * current T and P.
     *
     * The formulas used are from Sengers and Watson @cite sengers1986.
     *
     * The formulation is accurate for all temperatures and pressures, for steam
     * and for water, even near the critical point. Pressures above 500 MPa and
     * temperature above 900 C are suspect.
     */
    double viscosity() override;

    double bulkViscosity() override {
        return 0.0;
    }

    //! Returns the thermal conductivity of water at the current conditions
    //! (W/m/K)
    /*!
     * This function calculates the value of the thermal conductivity of water
     * at the current T and P.
     *
     * The formulas used are from Sengers and Watson @cite sengers1986.
     *
     * The formulation is accurate for all temperatures and pressures, for steam
     * and for water, even near the critical point. Pressures above 500 MPa and
     * temperature above 900 C are suspect.
     */
    double thermalConductivity() override;

    void init(ThermoPhase* thermo, int mode=0, int log_level=0) override;
};
}
#endif
