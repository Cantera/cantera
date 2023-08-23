/**
 *  @file PureFluidPhase.h
 *
 *   Header for a ThermoPhase class for a pure fluid phase consisting of
 *   gas, liquid, mixed-gas-liquid and supercritical fluid (see @ref thermoprops
 *   and class @link Cantera::PureFluidPhase PureFluidPhase@endlink).
 *
 * It inherits from ThermoPhase, but is built on top of the tpx package.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_EOS_TPX_H
#define CT_EOS_TPX_H

#include "ThermoPhase.h"
#include "cantera/tpx/Sub.h"

namespace Cantera
{
//! This phase object consists of a single component that can be a gas, a
//! liquid, a mixed gas-liquid fluid, or a fluid beyond its critical point
/*!
 * The object inherits from ThermoPhase. However, it's built on top of the tpx
 * package.
 *
 * @ingroup thermoprops
 */
class PureFluidPhase : public ThermoPhase
{
public:
    //! Empty Base Constructor
    PureFluidPhase() = default;

    string type() const override {
        return "pure-fluid";
    }

    //! String indicating the mechanical phase of the matter in this Phase.
    /*!
     * Options for the string are:
     *   * `supercritical`
     *   * `gas`
     *   * `liquid`
     *   * `liquid-gas-mix`
     *
     * If the temperature or pressure are greater than the critical temperature or
     * pressure, respectively, the mechanical phase is `supercritical`. If the
     * underlying tpx::TwoPhase() returns `True`, the mechanical phase is
     * `liquid-gas-mix`. If the temperature is greater than the saturation temperature
     * at the current pressure, the mechanical phase is `gas`. Otherwise, the mechanical
     * phase is `liquid`.
     */
    string phaseOfMatter() const override;

    //! Set the name of the TPX substance to use for the equation of state. This
    //! function should be called before initThermo().
    void setSubstance(const string& name) {
        m_tpx_name = name;
    }

    bool isPure() const override {
        return true;
    }

    bool hasPhaseTransition() const override {
        return true;
    }

    vector<string> fullStates() const override;
    vector<string> partialStates() const override;

    double minTemp(size_t k=npos) const override;
    double maxTemp(size_t k=npos) const override;

    double enthalpy_mole() const override;
    double intEnergy_mole() const override;
    double entropy_mole() const override;
    double gibbs_mole() const override;
    double cp_mole() const override;
    double cv_mole() const override;

    //! Return the thermodynamic pressure (Pa).
    /*!
     * This method calculates the current pressure consistent with the
     * independent variables, T, rho.
     */
    double pressure() const override;

    //! sets the thermodynamic pressure (Pa).
    /*!
     * This method calculates the density that is consistent with the
     * desired pressure, given the temperature.
     *
     * @param p  Pressure (Pa)
     */
    void setPressure(double p) override;
    void setTemperature(const double T) override;
    void setDensity(const double rho) override;

    void getChemPotentials(double* mu) const override{
        mu[0] = gibbs_mole();
    }

    void getPartialMolarEnthalpies(double* hbar) const override;
    void getPartialMolarEntropies(double* sbar) const override;
    void getPartialMolarIntEnergies(double* ubar) const override;
    void getPartialMolarCp(double* cpbar) const override;
    void getPartialMolarVolumes(double* vbar) const override;

    Units standardConcentrationUnits() const override;
    void getActivityConcentrations(double* c) const override;
    double standardConcentration(size_t k=0) const override;

    void getActivities(double* a) const override;

    double isothermalCompressibility() const override;
    double thermalExpansionCoeff() const override;

    //! Returns a reference to the substance object
    tpx::Substance& TPX_Substance();

    //! @name Properties of the Standard State of the Species in the Solution
    //!
    //! The standard state of the pure fluid is defined as the real properties
    //! of the pure fluid at the most stable state of the fluid at the current
    //! temperature and pressure of the solution. With this definition, the
    //! activity of the fluid is always then defined to be equal to one.
    //! @{

    void getStandardChemPotentials(double* mu) const override;
    void getEnthalpy_RT(double* hrt) const override;
    void getEntropy_R(double* sr) const override;
    void getGibbs_RT(double* grt) const override;

    //! @}
    //! @name Thermodynamic Values for the Species Reference States
    //!
    //! The species reference state for pure fluids is defined as an ideal gas at
    //! the reference pressure and current temperature of the fluid.
    //! @{

    void getEnthalpy_RT_ref(double* hrt) const override;
    void getGibbs_RT_ref(double* grt) const override;
    void getGibbs_ref(double* g) const override;
    void getEntropy_R_ref(double* er) const override;

    //! @}
    //! @name Setting the State
    //!
    //! These methods set all or part of the thermodynamic state.
    //! @{

    void setState_HP(double h, double p, double tol=1e-9) override;
    void setState_UV(double u, double v, double tol=1e-9) override;
    void setState_SV(double s, double v, double tol=1e-9) override;
    void setState_SP(double s, double p, double tol=1e-9) override;
    void setState_ST(double s, double t, double tol=1e-9) override;
    void setState_TV(double t, double v, double tol=1e-9) override;
    void setState_PV(double p, double v, double tol=1e-9) override;
    void setState_UP(double u, double p, double tol=1e-9) override;
    void setState_VH(double v, double h, double tol=1e-9) override;
    void setState_TH(double t, double h, double tol=1e-9) override;
    void setState_SH(double s, double h, double tol=1e-9) override;
    //! @}
    //! @name Critical State Properties
    //! @{

    double critTemperature() const override;
    double critPressure() const override;
    double critDensity() const override;

    //! @}
    //! @name Saturation properties.
    //! @{

    double satTemperature(double p) const override;
    double satPressure(double t) override;
    double vaporFraction() const override;

    void setState_Tsat(double t, double x) override;
    void setState_Psat(double p, double x) override;
    //! @}

    void initThermo() override;
    void getParameters(AnyMap& phaseNode) const override;

    string report(bool show_thermo=true, double threshold=1e-14) const override;

    bool compatibleWithMultiPhase() const override {
        return false;
    }

protected:
    //! Main call to the tpx level to set the state of the system
    /*!
     * @param n  Integer indicating which 2 thermo components are held constant
     * @param x  Value of the first component
     * @param y  Value of the second component
     */
    void Set(tpx::PropertyPair::type n, double x, double y) const;

private:
    //! Pointer to the underlying tpx object Substance that does the work
    mutable unique_ptr<tpx::Substance> m_sub;

    //! Name for this substance used by the TPX package
    string m_tpx_name;

    //! Molecular weight of the substance (kg kmol-1)
    double m_mw = -1.0;

    //! flag to turn on some printing.
    bool m_verbose = false;
};

}

#endif
