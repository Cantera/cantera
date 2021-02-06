/**
 *  @file PureFluidPhase.h
 *
 *   Header for a ThermoPhase class for a pure fluid phase consisting of
 *   gas, liquid, mixed-gas-liquid and supercrit fluid (see \ref thermoprops
 *   and class \link Cantera::PureFluidPhase PureFluidPhase\endlink).
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
    PureFluidPhase();

    virtual std::string type() const {
        return "PureFluid";
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
    virtual std::string phaseOfMatter() const;

    //! Set the name of the TPX substance to use for the equation of state. This
    //! function should be called before initThermo().
    void setSubstance(const std::string& name) {
        m_tpx_name = name;
    }

    virtual bool isPure() const {
        return true;
    }

    virtual bool hasPhaseTransition() const {
        return true;
    }

    virtual std::vector<std::string> fullStates() const;
    virtual std::vector<std::string> partialStates() const;

    virtual double minTemp(size_t k=npos) const;
    virtual double maxTemp(size_t k=npos) const;

    virtual doublereal enthalpy_mole() const;
    virtual doublereal intEnergy_mole() const;
    virtual doublereal entropy_mole() const;
    virtual doublereal gibbs_mole() const;
    virtual doublereal cp_mole() const;
    virtual doublereal cv_mole() const;

    //! Return the thermodynamic pressure (Pa).
    /*!
     * This method calculates the current pressure consistent with the
     * independent variables, T, rho.
     */
    virtual doublereal pressure() const;

    //! sets the thermodynamic pressure (Pa).
    /*!
     * This method calculates the density that is consistent with the
     * desired pressure, given the temperature.
     *
     * @param p  Pressure (Pa)
     */
    virtual void setPressure(doublereal p);
    virtual void setTemperature(const double T);
    virtual void setDensity(const double rho);

    virtual void getChemPotentials(doublereal* mu) const {
        mu[0] = gibbs_mole();
    }

    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;
    virtual void getPartialMolarEntropies(doublereal* sbar) const;
    virtual void getPartialMolarIntEnergies(doublereal* ubar) const;
    virtual void getPartialMolarCp(doublereal* cpbar) const;
    virtual void getPartialMolarVolumes(doublereal* vbar) const;

    virtual Units standardConcentrationUnits() const;
    virtual void getActivityConcentrations(doublereal* c) const;
    virtual doublereal standardConcentration(size_t k=0) const;

    virtual void getActivities(doublereal* a) const;

    virtual doublereal isothermalCompressibility() const;
    virtual doublereal thermalExpansionCoeff() const;

    //! Returns a reference to the substance object
    tpx::Substance& TPX_Substance();

    //@}
    /// @name Properties of the Standard State of the Species in the Solution
    /*!
     *  The standard state of the pure fluid is defined as the real properties
     *  of the pure fluid at the most stable state of the fluid at the current
     *  temperature and pressure of the solution. With this definition, the
     *  activity of the fluid is always then defined to be equal to one.
     */
    //@{

    virtual void getStandardChemPotentials(doublereal* mu) const;
    virtual void getEnthalpy_RT(doublereal* hrt) const;
    virtual void getEntropy_R(doublereal* sr) const;
    virtual void getGibbs_RT(doublereal* grt) const;

    //@}

    /// @name Thermodynamic Values for the Species Reference States
    /*!
     * The species reference state for pure fluids is defined as an ideal gas at
     * the reference pressure and current temperature of the fluid.
     */
    //@{

    virtual void getEnthalpy_RT_ref(doublereal* hrt) const;
    virtual void getGibbs_RT_ref(doublereal* grt) const;
    virtual void getGibbs_ref(doublereal* g) const;
    virtual void getEntropy_R_ref(doublereal* er) const;

    /**
     * @name Setting the State
     *
     * These methods set all or part of the thermodynamic state.
     * @{
     */

    virtual void setState_HP(double h, double p, double tol=1e-9);
    virtual void setState_UV(double u, double v, double tol=1e-9);
    virtual void setState_SV(double s, double v, double tol=1e-9);
    virtual void setState_SP(double s, double p, double tol=1e-9);
    virtual void setState_ST(double s, double t, double tol=1e-9);
    virtual void setState_TV(double t, double v, double tol=1e-9);
    virtual void setState_PV(double p, double v, double tol=1e-9);
    virtual void setState_UP(double u, double p, double tol=1e-9);
    virtual void setState_VH(double v, double h, double tol=1e-9);
    virtual void setState_TH(double t, double h, double tol=1e-9);
    virtual void setState_SH(double s, double h, double tol=1e-9);
    //@}

    //! @name Critical State Properties
    //@{

    virtual doublereal critTemperature() const;
    virtual doublereal critPressure() const;
    virtual doublereal critDensity() const;

    //@}

    //! @name Saturation properties.
    //@{

    virtual doublereal satTemperature(doublereal p) const;
    virtual doublereal satPressure(doublereal t);
    virtual doublereal vaporFraction() const;

    virtual void setState_Tsat(doublereal t, doublereal x);
    virtual void setState_Psat(doublereal p, doublereal x);
    //@}

    virtual void initThermo();
    virtual void setParametersFromXML(const XML_Node& eosdata);

    virtual std::string report(bool show_thermo=true,
                               doublereal threshold=1e-14) const;

    virtual bool compatibleWithMultiPhase() const {
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
    mutable std::unique_ptr<tpx::Substance> m_sub;

    //! Int indicating the type of the fluid
    /*!
     * The tpx package uses an int to indicate what fluid is being sought. Used
     * only if #m_tpx_name is not set.
     */
    int m_subflag;

    //! Name for this substance used by the TPX package. If this is not set,
    //! #m_subflag is used instead.
    std::string m_tpx_name;

    //! Molecular weight of the substance (kg kmol-1)
    doublereal m_mw;

    //! flag to turn on some printing.
    bool m_verbose;
};

}

#endif
