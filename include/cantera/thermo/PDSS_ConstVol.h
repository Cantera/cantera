/**
 *  @file PDSS_ConstVol.h
 *    Declarations for the class PDSS_ConstVol (pressure dependent standard state)
 *    which handles calculations for a single species with a constant molar volume in a phase
 *    (see class \ref pdssthermo and \link Cantera::PDSS_ConstVol PDSS_ConstVol\endlink).
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
    PDSS_ConstVol();

    //! @name Molar Thermodynamic Properties of the Species Standard State in
    //!     the Solution
    //! @{

    // See PDSS.h for documentation of functions overridden from Class PDSS
    virtual doublereal intEnergy_mole() const;
    virtual doublereal cv_mole() const;

    //! @}
    //! @name Mechanical Equation of State Properties
    //! @{

    virtual void setPressure(doublereal pres);
    virtual void setTemperature(doublereal temp);
    virtual void setState_TP(doublereal temp, doublereal pres);
    virtual void setState_TR(doublereal temp, doublereal rho);

    //! @}
    //!  @name  Miscellaneous properties of the standard state
    //! @{

    virtual doublereal satPressure(doublereal t);

    //! @}
    //! @name Initialization of the Object
    //! @{

    virtual void initThermo();
    virtual void setParametersFromXML(const XML_Node& speciesNode);

    //! Set the (constant) molar volume [m3/kmol] of the species. Must be called before
    //! initThermo().
    void setMolarVolume(double v) {
        m_constMolarVolume = v;
    }

    //@}

private:
    //! Value of the constant molar volume for the species
    /*!
     *    m3 / kmol
     */
    doublereal m_constMolarVolume;
};

}

#endif
