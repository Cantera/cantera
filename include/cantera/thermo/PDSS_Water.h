/**
 *  @file PDSS_Water.h
 * Implementation of a pressure dependent standard state
 * virtual function for a Pure Water Phase
 * (see @ref pdssthermo and class @link Cantera::PDSS_Water PDSS_Water@endlink).
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

    //! @name  Molar Thermodynamic Properties of the Species Standard State
    //! @{

    // See PDSS.h for documentation of functions overridden from Class PDSS

    double enthalpy_mole() const override;
    double intEnergy_mole() const override;
    double entropy_mole() const override;
    double gibbs_mole() const override;
    double cp_mole() const override;
    double cv_mole() const override;
    double molarVolume() const override;
    double density() const override;

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
    double pref_safe(double temp) const;

    double gibbs_RT_ref() const override;
    double enthalpy_RT_ref() const override;
    double entropy_R_ref() const override;
    double cp_R_ref() const override;
    double molarVolume_ref() const override;

    //! @}
    //! @name Mechanical Equation of State Properties
    //! @{

    double pressure() const override;
    void setPressure(double pres) override;
    void setTemperature(double temp) override;
    void setState_TP(double temp, double pres) override;

    //! Set the density of the water phase
    /*!
     *  This is a non-virtual function because it specific to this object.
     *
     * @param dens Density of the water (kg/m3)
     */
    void setDensity(double dens);

    double thermalExpansionCoeff() const override;

    //! Return the derivative of the volumetric thermal expansion coefficient.
    //! Units: 1/K2.
    /*!
     * The thermal expansion coefficient is defined as
     * @f[
     * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
     * @f]
     */
    double dthermalExpansionCoeffdT() const;

    //! Returns the isothermal compressibility. Units: 1/Pa.
    /*!
     * The isothermal compressibility is defined as
     * @f[
     * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
     * @f]
     *  or
     * @f[
     * \kappa_T = \frac{1}{\rho}\left(\frac{\partial \rho}{\partial P}\right)_T
     * @f]
     */
    double isothermalCompressibility() const;

    //! @}
    //! @name Miscellaneous properties of the standard state
    //! @{

    double critTemperature() const override;
    double critPressure() const override;
    double critDensity() const override;
    double satPressure(double t) override;

    //! Get a pointer to a changeable WaterPropsIAPWS object
    WaterPropsIAPWS* getWater() {
        return &m_sub;
    }

    //! Get a pointer to a changeable WaterPropsIAPWS object
    WaterProps* getWaterProps() {
        return &m_waterProps;
    }

    void getParameters(AnyMap& eosNode) const override;

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
    double m_dens;

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
    int m_iState = WATER_LIQUID;

    /**
     *  Offset constants used to obtain consistency with the NIST database.
     *  This is added to all internal energy and enthalpy results.
     *  units = J kmol-1.
     */
    double EW_Offset = 0.0;

    /**
     *  Offset constant used to obtain consistency with NIST convention.
     *  This is added to all internal entropy results.
     *  units = J kmol-1 K-1.
     */
    double SW_Offset = 0.0;

public:
    /**
     *  Since this phase represents a liquid phase, it's an error to
     *  return a gas-phase answer. However, if the below is true, then
     *  a gas-phase answer is allowed. This is used to check the thermodynamic
     *  consistency with ideal-gas thermo functions for example.
     */
    bool m_allowGasPhase = false;
};

}

#endif
