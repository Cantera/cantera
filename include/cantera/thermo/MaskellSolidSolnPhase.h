/**
 * @file MaskellSolidSolnPhase.h Header file for a solid solution model
 * following Maskell, Shaw, and Tye. Electrochimica Acta 1982
 *
 * This class inherits from the %Cantera class ThermoPhase and implements a
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
 *
 * @deprecated To be removed after %Cantera 3.0. This class has numerous thermodynamic
 *    inconsistencies. See https://github.com/Cantera/cantera/issues/1321.
 */
class MaskellSolidSolnPhase : public VPStandardStateTP
{
public:
    MaskellSolidSolnPhase();

    string type() const override {
        return "MaskellSolidsoln";
    }

    Units standardConcentrationUnits() const override { return Units(1.0); }
    void getActivityConcentrations(double* c) const override;
    double standardConcentration(size_t k=0) const override { return 1.0; }
    double logStandardConc(size_t k=0) const override { return 0.0; }

    //! @name Molar Thermodynamic Properties of the Solution
    //! @{

    double enthalpy_mole() const override;
    double entropy_mole() const override;

    //! @}
    //! @name Mechanical Equation of State Properties
    //!
    //! In this equation of state implementation, the density is a function only
    //! of the mole fractions. Therefore, it can't be an independent variable.
    //! Instead, the pressure is used as the independent variable. Functions
    //! which try to set the thermodynamic state by calling setDensity() will
    //! cause an exception to be thrown.
    //! @{

    /**
     * Pressure. Units: Pa.
     * For this incompressible system, we return the internally stored
     * independent value of the pressure.
     */
    double pressure() const override {
        return m_Pcurrent;
    }

    /**
     * Set the pressure at constant temperature. Units: Pa. This method sets a
     * constant within the object. The mass density is not a function of
     * pressure.
     *
     * @param p   Input Pressure (Pa)
     */
    void setPressure(double p) override;

    void calcDensity() override;

    //! @}
    //! @name Chemical Potentials and Activities
    //! @{

    void getActivityCoefficients(double* ac) const override;
    void getChemPotentials(double* mu) const override;
    //! @deprecated To be removed after %Cantera 3.0. Use getChemPotentials() instead.
    void getChemPotentials_RT(double* mu) const override;

    //! @}
    //! @name  Partial Molar Properties of the Solution
    //! @{

    void getPartialMolarEnthalpies(double* hbar) const override;
    void getPartialMolarEntropies(double* sbar) const override;
    void getPartialMolarCp(double* cpbar) const override;
    void getPartialMolarVolumes(double* vbar) const override;
    void getPureGibbs(double* gpure) const override;
    void getStandardChemPotentials(double* mu) const override;

    //! @}
    //! @name Utility Functions
    //! @{

    void initThermo() override;
    void getParameters(AnyMap& phaseNode) const override;

    void set_h_mix(const double hmix) { h_mixing = hmix; }

    //! Set the product Species. Must be called after species have been added.
    void setProductSpecies(const string& name);
    //! @}

private:
    /**
     * m_Pcurrent = The current pressure. Since the density isn't a function of
     * pressure, but only of the mole fractions, we need to independently
     * specify the pressure.
     */
    double m_Pcurrent = OneAtm;

    //! Value of the enthalpy change on mixing due to protons changing from type
    //! B to type A configurations.
    double h_mixing = 0.0;

    //! Index of the species whose mole fraction defines the extent of reduction r
    int product_species_index = -1;
    int reactant_species_index = -1;

    // Functions to calculate some of the pieces of the mixing terms.
    double s() const;
    double fm(const double r) const;
    double p(const double r) const;
};
}

#endif
