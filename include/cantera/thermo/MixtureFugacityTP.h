/**
 *  @file MixtureFugacityTP.h
 *    Header file for a derived class of ThermoPhase that handles
 *    non-ideal mixtures based on the fugacity models (see @ref thermoprops and
 *    class @link Cantera::MixtureFugacityTP MixtureFugacityTP@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MIXTUREFUGACITYTP_H
#define CT_MIXTUREFUGACITYTP_H

#include "ThermoPhase.h"

namespace Cantera
{
//! Various states of the Fugacity object. In general there can be multiple liquid
//! objects for a single phase identified with each species.

#define FLUID_UNSTABLE -4
#define FLUID_UNDEFINED -3
#define FLUID_SUPERCRIT -2
#define FLUID_GAS -1
#define FLUID_LIQUID_0 0
#define FLUID_LIQUID_1 1
#define FLUID_LIQUID_2 2
#define FLUID_LIQUID_3 3
#define FLUID_LIQUID_4 4
#define FLUID_LIQUID_5 5
#define FLUID_LIQUID_6 6
#define FLUID_LIQUID_7 7
#define FLUID_LIQUID_8 8
#define FLUID_LIQUID_9 9

/**
 * @ingroup thermoprops
 *
 * This is a filter class for ThermoPhase that implements some preparatory steps
 * for efficiently handling mixture of gases that whose standard states are
 * defined as ideal gases, but which describe also non-ideal solutions. In
 * addition a multicomponent liquid phase below the critical temperature of the
 * mixture is also allowed. The main subclass is currently a mixture Redlich-
 * Kwong class.
 *
 * Several concepts are introduced. The first concept is there are temporary
 * variables for holding the species standard state values of Cp, H, S, G, and V
 * at the last temperature and pressure called. These functions are not
 * recalculated if a new call is made using the previous temperature and
 * pressure.
 *
 * The other concept is that the current state of the mixture is tracked. The
 * state variable is either GAS, LIQUID, or SUPERCRIT fluid.  Additionally, the
 * variable LiquidContent is used and may vary between 0 and 1.
 *
 * Typically, only one liquid phase is allowed to be formed within these
 * classes. Additionally, there is an inherent contradiction between three phase
 * models and the ThermoPhase class. The ThermoPhase class is really only meant
 * to represent a single instantiation of a phase. The three phase models may be
 * in equilibrium with multiple phases of the fluid in equilibrium with each
 * other. This has yet to be resolved.
 *
 *  This class is usually used for non-ideal gases.
 */
class MixtureFugacityTP : public ThermoPhase
{
public:
    //! @name Constructors and Duplicators for MixtureFugacityTP
    //! @{

    //! Constructor.
    MixtureFugacityTP() = default;

    //! @}
    //! @name  Utilities
    //! @{

    string type() const override {
        return "MixtureFugacity";
    }

    int standardStateConvention() const override;

    //! Set the solution branch to force the ThermoPhase to exist on one branch
    //! or another
    /*!
     *  @param solnBranch  Branch that the solution is restricted to. the value
     *       -1 means gas. The value -2 means unrestricted. Values of zero or
     *       greater refer to species dominated condensed phases.
     */
    void setForcedSolutionBranch(int solnBranch);

    //! Report the solution branch which the solution is restricted to
    /*!
     *  @return Branch that the solution is restricted to. the value -1 means
     *      gas. The value -2 means unrestricted. Values of zero or greater
     *      refer to species dominated condensed phases.
     */
    int forcedSolutionBranch() const;

    //! Report the solution branch which the solution is actually on
    /*!
     *  @return Branch that the solution is restricted to. the value -1 means
     *      gas. The value -2 means superfluid.. Values of zero or greater refer
     *      to species dominated condensed phases.
     */
    int reportSolnBranchActual() const;

    //! @}
    //! @name Molar Thermodynamic properties
    //! @{

    double enthalpy_mole() const override;
    double entropy_mole() const override;

    //! @}
    //! @name  Properties of the Standard State of the Species in the Solution
    //!
    //! Within MixtureFugacityTP, these properties are calculated via a common
    //! routine, _updateStandardStateThermo(), which must be overloaded in
    //! inherited objects. The values are cached within this object, and are not
    //! recalculated unless the temperature or pressure changes.
    //! @{

    //! Get the array of chemical potentials at unit activity.
    /*!
     * These are the standard state chemical potentials @f$ \mu^0_k(T,P)
     * @f$. The values are evaluated at the current temperature and pressure.
     *
     * For all objects with the Mixture Fugacity approximation, we define the
     * standard state as an ideal gas at the current temperature and pressure
     * of the solution.
     *
     * @param mu   Output vector of standard state chemical potentials.
     *             length = m_kk. units are J / kmol.
     */
    void getStandardChemPotentials(double* mu) const override;

    //! Get the nondimensional Enthalpy functions for the species at their
    //! standard states at the current *T* and *P* of the solution.
    /*!
    * For all objects with the Mixture Fugacity approximation, we define the
    * standard state as an ideal gas at the current temperature and pressure
    * of the solution.
    *
    * @param hrt     Output vector of standard state enthalpies.
    *                length = m_kk. units are unitless.
    */
    void getEnthalpy_RT(double* hrt) const override;

    //! Get the array of nondimensional Enthalpy functions for the standard
    //! state species at the current *T* and *P* of the solution.
    /*!
     * For all objects with the Mixture Fugacity approximation, we define the
     * standard state as an ideal gas at the current temperature and pressure of
     * the solution.
     *
     * @param sr     Output vector of nondimensional standard state entropies.
     *               length = m_kk.
     */
    void getEntropy_R(double* sr) const override;

    //! Get the nondimensional Gibbs functions for the species at their standard
    //! states of solution at the current T and P of the solution.
    /*!
     * For all objects with the Mixture Fugacity approximation, we define the
     * standard state as an ideal gas at the current temperature and pressure
     * of the solution.
     *
     * @param grt    Output vector of nondimensional standard state Gibbs free
     *               energies. length = m_kk.
     */
    void getGibbs_RT(double* grt) const override;

    //! Get the pure Gibbs free energies of each species. Species are assumed to
    //! be in their standard states.
    /*!
     * This is the same as getStandardChemPotentials().
     *
     * @param[out] gpure   Array of standard state Gibbs free energies. length =
     *     m_kk. units are J/kmol.
     */
    void getPureGibbs(double* gpure) const override;

    //! Returns the vector of nondimensional internal Energies of the standard
    //! state at the current temperature and pressure of the solution for each
    //! species.
    /*!
     * For all objects with the Mixture Fugacity approximation, we define the
     * standard state as an ideal gas at the current temperature and pressure
     * of the solution.
     *
     * @f[
     *  u^{ss}_k(T,P) = h^{ss}_k(T)  - P * V^{ss}_k
     * @f]
     *
     * @param urt    Output vector of nondimensional standard state internal
     *               energies. length = m_kk.
     */
    void getIntEnergy_RT(double* urt) const override;

    //! Get the nondimensional Heat Capacities at constant pressure for the
    //! standard state of the species at the current T and P.
    /*!
     * For all objects with the Mixture Fugacity approximation, we define the
     * standard state as an ideal gas at the current temperature and pressure of
     * the solution.
     *
     * @param cpr    Output vector containing the the nondimensional Heat
     *               Capacities at constant pressure for the standard state of
     *               the species. Length: m_kk.
     */
    void getCp_R(double* cpr) const override;

    //! Get the molar volumes of each species in their standard states at the
    //! current *T* and *P* of the solution.
    /*!
     * For all objects with the Mixture Fugacity approximation, we define the
     * standard state as an ideal gas at the current temperature and pressure of
     * the solution.
     *
     * units = m^3 / kmol
     *
     * @param vol Output vector of species volumes. length = m_kk.
     *            units =  m^3 / kmol
     */
    void getStandardVolumes(double* vol) const override;
    //! @}

    //! Set the temperature of the phase
    /*!
     * Currently this passes down to setState_TP(). It does not make sense to
     * calculate the standard state without first setting T and P.
     *
     * @param temp  Temperature (kelvin)
     */
    void setTemperature(const double temp) override;

    //! Set the internally stored pressure (Pa) at constant temperature and
    //! composition
    /*!
     * Currently this passes down to setState_TP().  It does not make sense to
     * calculate the standard state without first setting T and P.
     *
     *  @param p input Pressure (Pa)
     */
    void setPressure(double p) override;

protected:
    void compositionChanged() override;

    //! Updates the reference state thermodynamic functions at the current T of
    //! the solution.
    /*!
     * This function must be called for every call to functions in this class.
     * It checks to see whether the temperature has changed and thus the ss
     * thermodynamics functions for all of the species must be recalculated.
     *
     * This function is responsible for updating the following internal members:
     *
     *  -  m_h0_RT;
     *  -  m_cp0_R;
     *  -  m_g0_RT;
     *  -  m_s0_R;
     */
    virtual void _updateReferenceStateThermo() const;

    //! Temporary storage - length = m_kk.
    mutable vector<double> m_tmpV;
public:

    //! @name Thermodynamic Values for the Species Reference States
    //!
    //! There are also temporary variables for holding the species reference-
    //! state values of Cp, H, S, and V at the last temperature and reference
    //! pressure called. These functions are not recalculated if a new call is
    //! made using the previous temperature. All calculations are done within the
    //! routine _updateRefStateThermo().
    //! @{

    void getEnthalpy_RT_ref(double* hrt) const override;
    void getGibbs_RT_ref(double* grt) const override;

protected:
    //! Returns the vector of nondimensional Gibbs free energies of the
    //! reference state at the current temperature of the solution and the
    //! reference pressure for the species.
    /*!
     * @return  Output vector contains the nondimensional Gibbs free energies
     *          of the reference state of the species
     *          length = m_kk, units = dimensionless.
     */
    const vector<double>& gibbs_RT_ref() const;

public:
    void getGibbs_ref(double* g) const override;
    void getEntropy_R_ref(double* er) const override;
    void getCp_R_ref(double* cprt) const override;
    void getStandardVolumes_ref(double* vol) const override;

    //! @}
    //! @name Initialization Methods - For Internal use
    //!
    //! The following methods are used in the process of constructing
    //! the phase and setting its parameters from a specification in an
    //! input file. They are not normally used in application programs.
    //! To see how they are used, see importPhase().
    //! @{
    bool addSpecies(shared_ptr<Species> spec) override;
    //! @}

protected:
    //! @name Special Functions for fugacity classes
    //! @{

    //!  Calculate the value of z
    /*!
     *  @f[
     *        z = \frac{P v}{R T}
     *  @f]
     *
     *  returns the value of z
     */
    double z() const;

    //! Calculate the deviation terms for the total entropy of the mixture from
    //! the ideal gas mixture
    /**
     * Here we use the current state conditions
     *
     * @returns the change in entropy in units of J kmol-1 K-1.
     */
    virtual double sresid() const;

    //! Calculate the deviation terms for the total enthalpy of the mixture from
    //! the ideal gas mixture
    /**
     * Here we use the current state conditions
     *
     * @returns the change in entropy in units of J kmol-1.
     */
    virtual double hresid() const;

    //! Estimate for the saturation pressure
    /*!
     * Note: this is only used as a starting guess for later routines that
     * actually calculate an accurate value for the saturation pressure.
     *
     * @param TKelvin  temperature in kelvin
     * @return the estimated saturation pressure at the given temperature
     */
    virtual double psatEst(double TKelvin) const;

public:
    //! Estimate for the molar volume of the liquid
    /*!
     * Note: this is only used as a starting guess for later routines that
     * actually calculate an accurate value for the liquid molar volume. This
     * routine doesn't change the state of the system.
     *
     * @param TKelvin  temperature in kelvin
     * @param pres     Pressure in Pa. This is used as an initial guess. If the
     *                 routine needs to change the pressure to find a stable
     *                 liquid state, the new pressure is returned in this
     *                 variable.
     * @returns the estimate of the liquid volume. If the liquid can't be
     *          found, this routine returns -1.
     */
    virtual double liquidVolEst(double TKelvin, double& pres) const;

    //! Calculates the density given the temperature and the pressure and a
    //! guess at the density.
    /*!
     * Note, below T_c, this is a multivalued function. We do not cross the
     * vapor dome in this. This is protected because it is called during
     * setState_TP() routines. Infinite loops would result if it were not
     * protected.
     *
     * @param TKelvin   Temperature in Kelvin
     * @param pressure  Pressure in Pascals (Newton/m**2)
     * @param phaseRequested     int representing the phase whose density we are
     *     requesting. If we put a gas or liquid phase here, we will attempt to
     *     find a volume in that part of the volume space, only, in this
     *     routine. A value of FLUID_UNDEFINED means that we will accept
     *     anything.
     * @param rhoguess   Guessed density of the fluid. A value of -1.0 indicates
     *     that there is no guessed density
     * @return    We return the density of the fluid at the requested phase. If
     *            we have not found any acceptable density we return a -1. If we
     *            have found an acceptable density at a different phase, we
     *            return a -2.
     */
    virtual double densityCalc(double TKelvin, double pressure, int phaseRequested,
                               double rhoguess);

protected:
    //! Utility routine in the calculation of the saturation pressure
    /*!
     * @param TKelvin        temperature (kelvin)
     * @param pres           pressure (Pascal)
     * @param[out] densLiq   density of liquid
     * @param[out] densGas   density of gas
     * @param[out] liqGRT    deltaG/RT of liquid
     * @param[out] gasGRT    deltaG/RT of gas
     */
    int corr0(double TKelvin, double pres, double& densLiq,
              double& densGas, double& liqGRT, double& gasGRT);

public:
    //! Returns the Phase State flag for the current state of the object
    /*!
     * @param checkState If true, this function does a complete check to see
     *        where in parameters space we are
     *
     *  There are three values:
     *  - WATER_GAS   below the critical temperature but below the critical density
     *  - WATER_LIQUID  below the critical temperature but above the critical density
     *  - WATER_SUPERCRIT   above the critical temperature
     */
    int phaseState(bool checkState = false) const;

    //! Return the value of the density at the liquid spinodal point (on the
    //! liquid side) for the current temperature.
    /*!
     * @returns the density with units of kg m-3
     */
    virtual double densSpinodalLiquid() const;

    //! Return the value of the density at the gas spinodal point (on the gas
    //! side) for the current temperature.
    /*!
     * @returns the density with units of kg m-3
     */
    virtual double densSpinodalGas() const;

public:
    //! Calculate the saturation pressure at the current mixture content for the
    //! given temperature
    /*!
     * @param TKelvin         (input) Temperature (Kelvin)
     * @param molarVolGas     (return) Molar volume of the gas
     * @param molarVolLiquid  (return) Molar volume of the liquid
     * @returns the saturation pressure at the given temperature
     */
    double calculatePsat(double TKelvin, double& molarVolGas, double& molarVolLiquid);

public:
    //! Calculate the saturation pressure at the current mixture content for the
    //! given temperature
    /*!
     * @param TKelvin   Temperature (Kelvin)
     * @return          The saturation pressure at the given temperature
     */
    double satPressure(double TKelvin) override;
    void getActivityConcentrations(double* c) const override;

protected:
    //! Calculate the pressure and the pressure derivative given the temperature
    //! and the molar volume
    /*!
     *  Temperature and mole number are held constant
     *
     * @param   TKelvin   temperature in kelvin
     * @param   molarVol  molar volume ( m3/kmol)
     * @param   presCalc  Returns the pressure.
     * @returns the derivative of the pressure wrt the molar volume
     */
    virtual double dpdVCalc(double TKelvin, double molarVol, double& presCalc) const;

    virtual void updateMixingExpressions();

    //! @}
    //! @name Critical State Properties.
    //! @{

    double critTemperature() const override;
    double critPressure() const override;
    double critVolume() const override;
    double critCompressibility() const override;
    double critDensity() const override;
    virtual void calcCriticalConditions(double& pc, double& tc, double& vc) const;

    //! Solve the cubic equation of state
    /*!
     *
     * Returns the number of solutions found. For the gas phase solution, it returns
     * a positive number (1 or 2). If it only finds the liquid branch solution,
     * it will return -1 or -2 instead of 1 or 2.
     * If it returns 0, then there is an error.
     * The cubic equation is solved using Nickalls' method @cite nickalls1993.
     *
     * @param   T         temperature (kelvin)
     * @param   pres      pressure (Pa)
     * @param   a         "a" parameter in the non-ideal EoS [Pa-m^6/kmol^2]
     * @param   b         "b" parameter in the non-ideal EoS [m^3/kmol]
     * @param   aAlpha    a*alpha (temperature dependent function for P-R EoS, 1 for R-K EoS)
     * @param   Vroot     Roots of the cubic equation for molar volume (m3/kmol)
     * @param   an        constant used in cubic equation
     * @param   bn        constant used in cubic equation
     * @param   cn        constant used in cubic equation
     * @param   dn        constant used in cubic equation
     * @param   tc        Critical temperature (kelvin)
     * @param   vc        Critical volume
     * @returns the number of solutions found
     */
    int solveCubic(double T, double pres, double a, double b,
                   double aAlpha, double Vroot[3], double an,
                   double bn, double cn, double dn, double tc, double vc) const;

    //! @}

    //! Storage for the current values of the mole fractions of the species
    /*!
     * This vector is kept up-to-date when some the setState functions are called.
     */
    vector<double> moleFractions_;

    //! Current state of the fluid
    /*!
     *  There are three possible states of the fluid:
     *  - FLUID_GAS
     *  - FLUID_LIQUID
     *  - FLUID_SUPERCRIT
     */
    int iState_ = FLUID_GAS;

    //! Force the system to be on a particular side of the spinodal curve
    int forcedState_ = FLUID_UNDEFINED;

    //! Temporary storage for dimensionless reference state enthalpies
    mutable vector<double> m_h0_RT;

    //! Temporary storage for dimensionless reference state heat capacities
    mutable vector<double> m_cp0_R;

    //! Temporary storage for dimensionless reference state Gibbs energies
    mutable vector<double> m_g0_RT;

    //! Temporary storage for dimensionless reference state entropies
    mutable vector<double> m_s0_R;
};
}

#endif
