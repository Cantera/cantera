/**
 *  @file PDSS_ConstVol.h
 *    Declarations for the class PDSS_ConstVol (pressure dependent standard state)
 *    which handles calculations for a single species with a constant molar volume in a phase
 *    (see class @ref pdssthermo and @link Cantera::PDSS_ConstVol PDSS_ConstVol@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PDSS_CONSTVOL_H
#define CT_PDSS_CONSTVOL_H

#include "PDSS.h"

namespace Cantera
{
//! Class for pressure dependent standard states that use a constant volume model
/*!
 * @ingroup pdssthermo
 */
class PDSS_ConstVol : public PDSS_Nondimensional
{
public:
    //! Default Constructor
    PDSS_ConstVol() = default;

    //! @name Molar Thermodynamic Properties of the Species Standard State
    //! @{

    // See PDSS.h for documentation of functions overridden from Class PDSS
    double intEnergy_mole() const override;
    double cv_mole() const override;

    //! @}
    //! @name Mechanical Equation of State Properties
    //! @{

    void setPressure(double pres) override;
    void setTemperature(double temp) override;
    void setState_TP(double temp, double pres) override;
    double satPressure(double t) override;

    //! @}
    //! @name Initialization of the Object
    //! @{

    void initThermo() override;
    void getParameters(AnyMap& eosNode) const override;

    //! Set the (constant) molar volume [m3/kmol] of the species. Must be called before
    //! initThermo().
    void setMolarVolume(double v) {
        m_constMolarVolume = v;
    }

    //! @}

private:
    //! Value of the constant molar volume for the species
    /*!
     *    m3 / kmol
     */
    double m_constMolarVolume;
};

}

#endif
