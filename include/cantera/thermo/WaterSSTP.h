/**
 *  @file WaterSSTP.h
 * Declares a ThermoPhase class consisting of pure water (see @ref thermoprops
 * and class @link Cantera::WaterSSTP WaterSSTP@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_WATERSSTP_H
#define CT_WATERSSTP_H

#include "SingleSpeciesTP.h"
#include "cantera/thermo/WaterPropsIAPWS.h"
#include "cantera/thermo/WaterProps.h"

namespace Cantera
{

class WaterPropsIAPWS;
class WaterProps;
//! Class for single-component water. This is designed to cover just the liquid
//! and supercritical phases of water.
/*!
 * The reference is W. Wagner, A. Pruss, "The IAPWS Formulation 1995 for the
 * Thermodynamic Properties of Ordinary Water Substance for General and
 * Scientific Use," J. Phys. Chem. Ref. Dat, 31, 387, 2002.
 *
 * ## Specification of Species Standard State Properties
 *
 * The offsets used in the steam tables are different than NIST's. They assume
 * u_liq(TP) = 0.0, s_liq(TP) = 0.0, where TP is the triple point conditions:
 *
 *      -  u(273.16, rho)    = 0.0
 *      -  s(273.16, rho)    = 0.0
 *      -  psat(273.16)      = 611.655 Pascal
 *      -  rho(273.16, psat) = 999.793 kg m-3
 *
 * These "steam table" assumptions are used by the WaterPropsIAPWS class.
 * Therefore, offsets must be calculated to make the thermodynamic properties
 * calculated within this class to be consistent with thermo properties within
 * Cantera.
 *
 * The thermodynamic base state for water is set to the NIST basis here by
 * specifying constants, #EW_Offset and #SW_Offset, one for energy quantities
 * and one for entropy quantities. The offsets are specified so that the
 * following properties hold:
 *
 * - Delta_Hfo_idealgas(298.15) = -241.826 kJ/gmol
 * - So_idealgas(298.15, 1bar)  =  188.835 J/gmolK
 *
 * (From http://webbook.nist.gov)
 *
 * The "o" here refers to a hypothetical ideal gas state. The way we achieve
 * this in practice is to evaluate at a very low pressure and then use the
 * theoretical ideal gas results to scale up to higher pressures:
 *
 * Ho(1bar) = H(P0)
 *
 * So(1bar) = S(P0) + RT ln(1bar/P0)
 *
 * ## Application within Kinetics Managers
 *
 * This is unimplemented.
 *
 * @ingroup thermoprops
 */
class WaterSSTP : public SingleSpeciesTP
{
public:
    //! Full constructor for a water phase
    /*!
     * @param inputFile String name of the input file
     * @param id        string id of the phase name
     */
    explicit WaterSSTP(const string& inputFile="", const string& id="");

    string type() const override {
        return "liquid-water-IAPWS95";
    }

    string phaseOfMatter() const override;

    //! @name  Molar Thermodynamic Properties of the Solution
    //! @{

    double cv_mole() const override;

    //! @}
    //! @name Mechanical Equation of State Properties
    //! @{

    double pressure() const override;
    void setPressure(double p) override;
    double isothermalCompressibility() const override;
    double thermalExpansionCoeff() const override;

    //! Return the derivative of the volumetric thermal expansion coefficient.
    //! Units: 1/K2.
    double dthermalExpansionCoeffdT() const;

    //! @}
    //! @name Properties of the Standard State of the Species in the Solution
    //! @{

    void getStandardChemPotentials(double* gss) const override;
    void getGibbs_RT(double* grt) const override;
    void getEnthalpy_RT(double* hrt) const override;
    void getEntropy_R(double* sr) const override;
    void getCp_R(double* cpr) const override;
    void getIntEnergy_RT(double* urt) const override;

    //! @}
    //! @name Thermodynamic Values for the Species Reference State
    //!
    //! All functions in this group need to be overridden, because the
    //! m_spthermo MultiSpeciesThermo function is not adequate for the real
    //! equation of state.
    //! @{

    void getEnthalpy_RT_ref(double* hrt) const override;
    void getGibbs_RT_ref(double* grt) const override;
    void getGibbs_ref(double* g) const override;
    void getEntropy_R_ref(double* er) const override;
    void getCp_R_ref(double* cprt) const override;
    void getStandardVolumes_ref(double* vol) const override;
    //! @}

    double critTemperature() const override;
    double critPressure() const override;
    double critDensity() const override;

    double satPressure(double t) override;

    bool compatibleWithMultiPhase() const override{
        return false;
    }

    //! Return the fraction of vapor at the current conditions
    /*!
     * Below Tcrit, this routine will always return 0, by definition of the
     * functionality of the routine. Above Tcrit, we query the density to toggle
     * between 0 and 1.
     */
    double vaporFraction() const override;

    //! Set the temperature of the phase
    /*!
     * The density and composition of the phase is constant during this
     * operator.
     *
     * @param temp Temperature (Kelvin)
     */
    void setTemperature(const double temp) override;

    //! Set the density of the phase
    /*!
     * The temperature and composition of the phase is constant during this
     * operator.
     *
     * @param dens value of the density in kg m-3
     */
    void setDensity(const double dens) override;

    void initThermo() override;

    //! Get a pointer to a changeable WaterPropsIAPWS object
    WaterPropsIAPWS* getWater() {
        return &m_sub;
    }

    //! Get a pointer to a changeable WaterPropsIAPWS object
    WaterProps* getWaterProps() {
        return m_waterProps.get();
    }

    //! Switch that enables calculations in the gas phase
    /**
     *  Since this phase represents a liquid (or supercritical) phase, it is an
     *  error to return a gas-phase answer. The sole intended use for this
     *  member function is to check the thermodynamic consistency of the
     *  underlying WaterProps class with ideal-gas thermo functions.
     */
    void _allowGasPhase(bool flag) { m_allowGasPhase = flag; }

protected:
    //! This routine must be overridden because it is not applicable.
    void _updateThermo() const;

private:
    //! WaterPropsIAPWS that calculates the real properties of water.
    mutable WaterPropsIAPWS m_sub;

    //! Pointer to the WaterProps object
    /*!
     * This class is used to house several approximation routines for properties
     * of water. This object owns m_waterProps, and the WaterPropsIAPWS object
     * used by WaterProps is m_sub, which is defined above.
     */
    unique_ptr<WaterProps> m_waterProps;

    //! Molecular weight of Water -> %Cantera assumption
    double m_mw = 0.0;

    //! Offset constants used to obtain consistency with the NIST database.
    /*!
     *  This is added to all internal energy and enthalpy results.
     *  units = J kmol-1.
     */
    double EW_Offset = 0.0;

    //! Offset constant used to obtain consistency with NIST convention.
    /*!
     *  This is added to all internal entropy results.
     *  units = J kmol-1 K-1.
     */
    double SW_Offset = 0.0;

    //! Boolean is true if object has been properly initialized for calculation
    bool m_ready = false;

    /**
     *  Since this phase represents a liquid (or supercritical) phase, it is an
     *  error to return a gas-phase answer. However, if the below is true, then
     *  a gas-phase answer is allowed. This is used to check the thermodynamic
     *  consistency with ideal-gas thermo functions for example.
     */
    bool m_allowGasPhase = false;
};

}

#endif
