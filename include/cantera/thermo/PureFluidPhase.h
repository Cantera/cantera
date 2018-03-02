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
// at http://www.cantera.org/license.txt for license and copyright information.

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

    //! Set the name of the TPX substance to use for the equation of state. This
    //! function should be called before initThermo().
    void setSubstance(const std::string& name) {
        m_tpx_name = name;
    }

    virtual double minTemp(size_t k=npos) const;
    virtual double maxTemp(size_t k=npos) const;

    virtual double enthalpy_mole() const;
    virtual double intEnergy_mole() const;
    virtual double entropy_mole() const;
    virtual double gibbs_mole() const;
    virtual double cp_mole() const;
    virtual double cv_mole() const;

    //! Return the thermodynamic pressure (Pa).
    /*!
     * This method calculates the current pressure consistent with the
     * independent variables, T, rho.
     */
    virtual double pressure() const;

    //! sets the thermodynamic pressure (Pa).
    /*!
     * This method calculates the density that is consistent with the
     * desired pressure, given the temperature.
     *
     * @param p  Pressure (Pa)
     */
    virtual void setPressure(double p);

    virtual void getChemPotentials(double* mu) const {
        mu[0] = gibbs_mole();
    }

    virtual void getPartialMolarEnthalpies(double* hbar) const;
    virtual void getPartialMolarEntropies(double* sbar) const;
    virtual void getPartialMolarIntEnergies(double* ubar) const;
    virtual void getPartialMolarCp(double* cpbar) const;
    virtual void getPartialMolarVolumes(double* vbar) const;

    virtual void getActivityConcentrations(double* c) const;
    virtual double standardConcentration(size_t k=0) const;

    virtual void getActivities(double* a) const;

    virtual double isothermalCompressibility() const;
    virtual double thermalExpansionCoeff() const;

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

    virtual void getStandardChemPotentials(double* mu) const;
    virtual void getEnthalpy_RT(double* hrt) const;
    virtual void getEntropy_R(double* sr) const;
    virtual void getGibbs_RT(double* grt) const;

    //@}

    /// @name Thermodynamic Values for the Species Reference States
    /*!
     * The species reference state for pure fluids is defined as an ideal gas at
     * the reference pressure and current temperature of the fluid.
     */
    //@{

    virtual void getEnthalpy_RT_ref(double* hrt) const;
    virtual void getGibbs_RT_ref(double* grt) const;
    virtual void getGibbs_ref(double* g) const;
    virtual void getEntropy_R_ref(double* er) const;

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

    virtual double critTemperature() const;
    virtual double critPressure() const;
    virtual double critDensity() const;

    //@}

    //! @name Saturation properties.
    //@{

    virtual double satTemperature(double p) const;
    virtual double satPressure(double t);
    virtual double vaporFraction() const;

    virtual void setState_Tsat(double t, double x);
    virtual void setState_Psat(double p, double x);
    //@}

    virtual void initThermo();
    virtual void setParametersFromXML(const XML_Node& eosdata);

    virtual std::string report(bool show_thermo=true,
                               double threshold=1e-14) const;

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

    //! Sets the state using a TPX::TV call
    void setTPXState() const;

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
    double m_mw;

    //! flag to turn on some printing.
    bool m_verbose;
};

}

#endif
