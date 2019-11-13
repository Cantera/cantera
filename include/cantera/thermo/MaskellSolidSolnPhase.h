/**
 * @file MaskellSolidSolnPhase.h Header file for a solid solution model
 * following Maskell, Shaw, and Tye. Electrochimica Acta 1982
 *
 * This class inherits from the Cantera class ThermoPhase and implements a
 * non-ideal solid solution model with incompressible thermodynamics.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MASKELLSOLIDSOLNPHASE_H
#define CT_MASKELLSOLIDSOLNPHASE_H

#include "VPStandardStateTP.h"

namespace Cantera
{
/**
 * Class MaskellSolidSolnPhase represents a condensed phase non-ideal solution
 * with 2 species following the thermodynamic model described in Maskell, Shaw,
 * and Tye, Manganese Dioxide Electrode -- IX, Electrochimica Acta 28(2) pp
 * 231-235, 1983.
 *
 * @ingroup thermoprops
 */
class MaskellSolidSolnPhase : public VPStandardStateTP
{
public:
    MaskellSolidSolnPhase();

    virtual std::string type() const {
        return "MaskellSolidsoln";
    }

    virtual Units standardConcentrationUnits() const { return Units(1.0); }
    virtual void getActivityConcentrations(doublereal* c) const;
    virtual doublereal standardConcentration(size_t k=0) const { return 1.0; }
    virtual doublereal logStandardConc(size_t k=0) const { return 0.0; }

    //! @name Molar Thermodynamic Properties of the Solution
    //! @{

    virtual doublereal enthalpy_mole() const;
    virtual doublereal entropy_mole() const;

    //@}
    /** @name Mechanical Equation of State Properties
     *
     * In this equation of state implementation, the density is a function only
     * of the mole fractions. Therefore, it can't be an independent variable.
     * Instead, the pressure is used as the independent variable. Functions
     * which try to set the thermodynamic state by calling setDensity() will
     * cause an exception to be thrown.
     */
    //@{

    /**
     * Pressure. Units: Pa.
     * For this incompressible system, we return the internally stored
     * independent value of the pressure.
     */
    virtual doublereal pressure() const {
        return m_Pcurrent;
    }

    /**
     * Set the pressure at constant temperature. Units: Pa. This method sets a
     * constant within the object. The mass density is not a function of
     * pressure.
     *
     * @param p   Input Pressure (Pa)
     */
    virtual void setPressure(doublereal p);

    /**
     * Overridden setDensity() function is necessary because the density is not
     * an independent variable.
     *
     * This function will now throw an error condition
     *
     * @param rho  Input density
     * @deprecated Functionality merged with base function after Cantera 2.5.
     *             (superseded by isCompressible check in Phase::setDensity)
     */
    virtual void setDensity(const doublereal rho);

    virtual void calcDensity();

    /**
     * Overridden setMolarDensity() function is necessary because the density
     * is not an independent variable.
     *
     * This function will now throw an error condition.
     *
     * @param rho   Input Density
     * @deprecated Functionality merged with base function after Cantera 2.5.
     *             (superseded by isCompressible check in Phase::setDensity)
     */
    virtual void setMolarDensity(const doublereal rho);

    //@}

    /**
     * @name Chemical Potentials and Activities
     * @{
     */

    virtual void getActivityCoefficients(doublereal* ac) const;
    virtual void getChemPotentials(doublereal* mu) const;
    virtual void getChemPotentials_RT(doublereal* mu) const;

    //@}
    /// @name  Partial Molar Properties of the Solution
    //@{

    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;
    virtual void getPartialMolarEntropies(doublereal* sbar) const;
    virtual void getPartialMolarCp(doublereal* cpbar) const;
    virtual void getPartialMolarVolumes(doublereal* vbar) const;
    virtual void getPureGibbs(doublereal* gpure) const;
    virtual void getStandardChemPotentials(doublereal* mu) const;

    //@}
    /// @name Utility Functions
    //@{

    virtual void initThermo();
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

    void set_h_mix(const doublereal hmix) { h_mixing = hmix; }

    //! Set the product Species. Must be called after species have been added.
    void setProductSpecies(const std::string& name);
    //@}

private:
    /**
     * m_Pcurrent = The current pressure. Since the density isn't a function of
     * pressure, but only of the mole fractions, we need to independently
     * specify the pressure.
     */
    doublereal m_Pcurrent;

    //! Value of the enthalpy change on mixing due to protons changing from type
    //! B to type A configurations.
    doublereal h_mixing;

    //! Index of the species whose mole fraction defines the extent of reduction r
    int product_species_index;
    int reactant_species_index;

    // Functions to calculate some of the pieces of the mixing terms.
    doublereal s() const;
    doublereal fm(const doublereal r) const;
    doublereal p(const doublereal r) const;
};
}

#endif
