/**
 *  @file PDSS_IdealGas.h
 *   Declarations for the class PDSS_IdealGas (pressure dependent standard state)
 *    which handles calculations for a single ideal gas species in a phase
 *    (see @ref pdssthermo and class @link Cantera::PDSS_IdealGas PDSS_IdealGas@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PDSS_IDEALGAS_H
#define CT_PDSS_IDEALGAS_H

#include "PDSS.h"

namespace Cantera
{
//! Derived class for pressure dependent standard states of an ideal gas species
/*!
 * This class is for a single Ideal Gas species.
 *
 * @ingroup pdssthermo
 * @deprecated To be removed after %Cantera 3.0.
 */
class PDSS_IdealGas : public PDSS_Nondimensional
{
public:
    //! Default Constructor
    PDSS_IdealGas();

    //! @name Molar Thermodynamic Properties of the Species Standard State
    //! @{

    // See PDSS.h for documentation of functions overridden from Class PDSS

    double intEnergy_mole() const override;
    double cv_mole() const override;

    //! @}
    //! @name Mechanical Equation of State Properties
    //! @{

    double pressure() const override;
    void setPressure(double pres) override;
    void setTemperature(double temp) override;
    void setState_TP(double temp, double pres) override;
    void setState_TR(double temp, double rho) override;

    //! @}
    //! @name Initialization of the Object
    //! @{

    void initThermo() override;
    void getParameters(AnyMap& eosNode) const override;
    //! @}
};
}

#endif
