/**
 *  @file WaterTransport.h Header file defining class WaterTransport
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_WATERTRAN_H
#define CT_WATERTRAN_H

#include "cantera/thermo/WaterPropsIAPWS.h"
#include "cantera/transport/TransportBase.h"

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
     */
    WaterTransport(thermo_t* thermo = 0, int ndim = 1);

    virtual std::string transportType() const {
        return "Water";
    }

    //! Returns the viscosity of water at the current conditions (kg/m/s)
    /*!
     * This function calculates the value of the viscosity of pure water at the
     * current T and P.
     *
     * The formulas used are from the paper: J. V. Sengers, J. T. R. Watson,
     * "Improved International Formulations for the Viscosity and Thermal
     * Conductivity of Water Substance", J. Phys. Chem. Ref. Data, 15, 1291
     * (1986).
     *
     * The formulation is accurate for all temperatures and pressures, for steam
     * and for water, even near the critical point. Pressures above 500 MPa and
     * temperature above 900 C are suspect.
     */
    virtual doublereal viscosity();

    virtual doublereal bulkViscosity() {
        return 0.0;
    }

    //! Returns the thermal conductivity of water at the current conditions
    //! (W/m/K)
    /*!
     * This function calculates the value of the thermal conductivity of water
     * at the current T and P.
     *
     * The formulas used are from the paper: J. V. Sengers, J. T. R. Watson,
     * "Improved International Formulations for the Viscosity and Thermal
     * Conductivity of Water Substance", J. Phys. Chem. Ref. Data, 15, 1291
     * (1986).
     *
     * The formulation is accurate for all temperatures and pressures, for steam
     * and for water, even near the critical point. Pressures above 500 MPa and
     * temperature above 900 C are suspect.
     */
    virtual doublereal thermalConductivity();

    virtual void init(thermo_t* thermo, int mode=0, int log_level=0);
};
}
#endif
