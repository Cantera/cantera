/**
 *  @file PDSS_Water.h
 * Implementation of a pressure dependent standard state
 * virtual function for a Pure Water Phase
 * (see \ref pdssthermo and class \link Cantera::PDSS_Water PDSS_Water\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PDSS_WATER_H
#define CT_PDSS_WATER_H

#include "PDSS.h"
#include "WaterPropsIAPWS.h"
#include "WaterProps.h"

namespace Cantera
{
//!  Class for the liquid water pressure dependent
//!  standard state
/*!
 * Notes:
 *
 * Base state for thermodynamic properties:
 *
 * The thermodynamic base state for water is set to the NIST basis here by
 * specifying constants EW_Offset and SW_Offset. These offsets are specified so
 * that the following properties hold:
 *
 * Delta_Hfo_gas(298.15) = -241.826 kJ/gmol
 * So_gas(298.15, 1bar)  = 188.835 J/gmolK
 *
 * (http://webbook.nist.gov)
 *
 * The "o" here refers to a hypothetical ideal gas state. The way we achieve
 * this in practice is to evaluate at a very low pressure and then use the
 * theoretical ideal gas results to scale up to higher pressures:
 *
 * Ho(1bar) = H(P0)
 *
 * So(1bar) = S(P0) + RT ln(1bar/P0)
 *
 * The offsets used in the steam tables are different than NIST's. They assume
 * u_liq(TP) = 0.0, s_liq(TP) = 0.0, where TP is the triple point conditions.
 *
 * @ingroup pdssthermo
 */
class PDSS_Water : public PDSS_Molar
{
public:
    //! Default constructor
    PDSS_Water();

    //! @name  Molar Thermodynamic Properties of the Species Standard State in the Solution
    //! @{

    // See PDSS.h for documentation of functions overridden from Class PDSS

    virtual doublereal enthalpy_mole() const;
    virtual doublereal intEnergy_mole() const;
    virtual doublereal entropy_mole() const;
    virtual doublereal gibbs_mole() const;
    virtual doublereal cp_mole() const;
    virtual doublereal cv_mole() const;
    virtual doublereal molarVolume() const;
    virtual doublereal density() const;

    //! @}
    //! @name Properties of the Reference State of the Species in the Solution
    //! @{

    //! Returns a reference pressure value that can be safely calculated by the
    //! underlying real equation of state for water
    /*!
     * Note, this function is needed because trying to calculate a one atm value
     * around the critical point will cause a crash
     *
     * @param temp  Temperature (Kelvin)
     */
    doublereal pref_safe(doublereal temp) const;

    virtual doublereal gibbs_RT_ref() const;
    virtual doublereal enthalpy_RT_ref() const;
    virtual doublereal entropy_R_ref() const;
    virtual doublereal cp_R_ref() const;
    virtual doublereal molarVolume_ref() const;

    //! @}
    //! @name Mechanical Equation of State Properties
    //! @{

    virtual doublereal pressure() const;
    virtual void setPressure(doublereal pres);
    virtual void setTemperature(doublereal temp);
    virtual void setState_TP(doublereal temp, doublereal pres);
    virtual void setState_TR(doublereal temp, doublereal rho);

    //! Set the density of the water phase
    /*!
     *  This is a non-virtual function because it specific to this object.
     *
     * @param dens Density of the water (kg/m3)
     */
    void setDensity(doublereal dens);

    virtual doublereal thermalExpansionCoeff() const;

    //! Return the derivative of the volumetric thermal expansion coefficient.
    //! Units: 1/K2.
    /*!
     * The thermal expansion coefficient is defined as
     * \f[
     * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
     * \f]
     */
    virtual doublereal dthermalExpansionCoeffdT() const;

    //! Returns the isothermal compressibility. Units: 1/Pa.
    /*!
     * The isothermal compressibility is defined as
     * \f[
     * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
     * \f]
     *  or
     * \f[
     * \kappa_T = \frac{1}{\rho}\left(\frac{\partial \rho}{\partial P}\right)_T
     * \f]
     */
    virtual doublereal isothermalCompressibility() const;

    //! @}
    //! @name Miscellaneous properties of the standard state
    //! @{

    virtual doublereal critTemperature() const;
    virtual doublereal critPressure() const;
    virtual doublereal critDensity() const;
    virtual doublereal satPressure(doublereal t);

    //! Get a pointer to a changeable WaterPropsIAPWS object
    WaterPropsIAPWS* getWater() {
        return &m_sub;
    }

    //! Get a pointer to a changeable WaterPropsIAPWS object
    WaterProps* getWaterProps() {
        return &m_waterProps;
    }

    //! @}

private:
    //! Pointer to the WaterPropsIAPWS object, which does the actual calculations
    //! for the real equation of state
    /*!
     * This object owns m_sub
     */
    mutable WaterPropsIAPWS m_sub;

    //! Pointer to the WaterProps object
    /*!
     *   This class is used to house several approximation
     *   routines for properties of water.
     *
     * This object owns m_waterProps, and the WaterPropsIAPWS object used by
     * WaterProps is m_sub, which is defined above.
     */
    WaterProps m_waterProps;

    //! State of the system - density
    /*!
     * Density is the independent variable here, but it's hidden behind the
     * object's interface.
     */
    doublereal m_dens;

    //! state of the fluid
    /*!
     *  @code
     *    0  WATER_GAS
     *    1  WATER_LIQUID
     *    2  WATER_SUPERCRIT
     *    3  WATER_UNSTABLELIQUID
     *    4  WATER_UNSTABLEGAS
     *  @endcode
     */
    int m_iState;

    /**
     *  Offset constants used to obtain consistency with the NIST database.
     *  This is added to all internal energy and enthalpy results.
     *  units = J kmol-1.
     */
    doublereal EW_Offset;

    /**
     *  Offset constant used to obtain consistency with NIST convention.
     *  This is added to all internal entropy results.
     *  units = J kmol-1 K-1.
     */
    doublereal SW_Offset;

public:
    /**
     *  Since this phase represents a liquid phase, it's an error to
     *  return a gas-phase answer. However, if the below is true, then
     *  a gas-phase answer is allowed. This is used to check the thermodynamic
     *  consistency with ideal-gas thermo functions for example.
     */
    bool m_allowGasPhase;
};

}

#endif
