/**
 *  @file MixtureFugacityTP.h
 *    Header file for a derived class of ThermoPhase that handles
 *    non-ideal mixtures based on the fugacity models (see \ref thermoprops and
 *    class \link Cantera::MixtureFugacityTP MixtureFugacityTP\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MIXTUREFUGACITYTP_H
#define CT_MIXTUREFUGACITYTP_H

#include "ThermoPhase.h"
#include "cantera/numerics/ResidEval.h"

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
 * @attention This class currently does not have any test cases or examples. Its
 *     implementation may be incomplete, and future changes to Cantera may
 *     unexpectedly cause this class to stop working. If you use this class,
 *     please consider contributing examples or test cases. In the absence of
 *     new tests or examples, this class may be deprecated and removed in a
 *     future version of Cantera. See
 *     https://github.com/Cantera/cantera/issues/267 for additional information.
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
    MixtureFugacityTP();

    //! @}
    //! @name  Utilities
    //! @{

    virtual std::string type() const {
        return "MixtureFugacity";
    }

    virtual int standardStateConvention() const;

    //! Set the solution branch to force the ThermoPhase to exist on one branch
    //! or another
    /*!
     *  @param solnBranch  Branch that the solution is restricted to. the value
     *       -1 means gas. The value -2 means unrestricted. Values of zero or
     *       greater refer to species dominated condensed phases.
     */
    virtual void setForcedSolutionBranch(int solnBranch);

    //! Report the solution branch which the solution is restricted to
    /*!
     *  @return Branch that the solution is restricted to. the value -1 means
     *      gas. The value -2 means unrestricted. Values of zero or greater
     *      refer to species dominated condensed phases.
     */
    virtual int forcedSolutionBranch() const;

    //! Report the solution branch which the solution is actually on
    /*!
     *  @return Branch that the solution is restricted to. the value -1 means
     *      gas. The value -2 means superfluid.. Values of zero or greater refer
     *      to species dominated condensed phases.
     */
    virtual int reportSolnBranchActual() const;

    virtual void getdlnActCoeffdlnN_diag(doublereal* dlnActCoeffdlnN_diag) const {
        throw NotImplementedError("MixtureFugacityTP::getdlnActCoeffdlnN_diag");
    }

    //@}
    /// @name  Partial Molar Properties of the Solution
    //@{

    //! Get the array of non-dimensional species chemical potentials
    //! These are partial molar Gibbs free energies.
    /*!
     * \f$ \mu_k / \hat R T \f$.
     * Units: unitless
     *
     * We close the loop on this function, here, calling getChemPotentials() and
     * then dividing by RT. No need for child classes to handle.
     *
     * @param mu    Output vector of non-dimensional species chemical potentials
     *              Length: m_kk.
     */
    virtual void getChemPotentials_RT(doublereal* mu) const;

    //@}
    /*!
     * @name  Properties of the Standard State of the Species in the Solution
     *
     * Within MixtureFugacityTP, these properties are calculated via a common
     * routine, _updateStandardStateThermo(), which must be overloaded in
     * inherited objects. The values are cached within this object, and are not
     * recalculated unless the temperature or pressure changes.
     */
    //@{

    //! Get the array of chemical potentials at unit activity.
    /*!
     * These are the standard state chemical potentials \f$ \mu^0_k(T,P)
     * \f$. The values are evaluated at the current temperature and pressure.
     *
     * For all objects with the Mixture Fugacity approximation, we define the
     * standard state as an ideal gas at the current temperature and pressure
     * of the solution.
     *
     * @param mu   Output vector of standard state chemical potentials.
     *             length = m_kk. units are J / kmol.
     */
    virtual void getStandardChemPotentials(doublereal* mu) const;

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
    virtual void getEnthalpy_RT(doublereal* hrt) const;

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
    virtual void getEntropy_R(doublereal* sr) const;

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
    virtual void getGibbs_RT(doublereal* grt) const;

    //! Get the pure Gibbs free energies of each species. Species are assumed to
    //! be in their standard states.
    /*!
     * This is the same as getStandardChemPotentials().
     *
     * @param[out] gpure   Array of standard state Gibbs free energies. length =
     *     m_kk. units are J/kmol.
     */
    virtual void getPureGibbs(doublereal* gpure) const;

    //! Returns the vector of nondimensional internal Energies of the standard
    //! state at the current temperature and pressure of the solution for each
    //! species.
    /*!
     * For all objects with the Mixture Fugacity approximation, we define the
     * standard state as an ideal gas at the current temperature and pressure
     * of the solution.
     *
     * \f[
     *  u^{ss}_k(T,P) = h^{ss}_k(T)  - P * V^{ss}_k
     * \f]
     *
     * @param urt    Output vector of nondimensional standard state internal
     *               energies. length = m_kk.
     */
    virtual void getIntEnergy_RT(doublereal* urt) const;

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
    virtual void getCp_R(doublereal* cpr) const;

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
    virtual void getStandardVolumes(doublereal* vol) const;
    // @}

    //! Set the temperature of the phase
    /*!
     * Currently this passes down to setState_TP(). It does not make sense to
     * calculate the standard state without first setting T and P.
     *
     * @param temp  Temperature (kelvin)
     */
    virtual void setTemperature(const doublereal temp);

    //! Set the internally stored pressure (Pa) at constant temperature and
    //! composition
    /*!
     * Currently this passes down to setState_TP().  It does not make sense to
     * calculate the standard state without first setting T and P.
     *
     *  @param p input Pressure (Pa)
     */
    virtual void setPressure(doublereal p);

protected:
    /**
     * Calculate the density of the mixture using the partial molar volumes and
     * mole fractions as input
     *
     * The formula for this is
     *
     * \f[
     *     \rho = \frac{\sum_k{X_k W_k}}{\sum_k{X_k V_k}}
     * \f]
     *
     * where \f$X_k\f$ are the mole fractions, \f$W_k\f$ are the molecular
     * weights, and \f$V_k\f$ are the pure species molar volumes.
     *
     * Note, the basis behind this formula is that in an ideal solution the
     * partial molar volumes are equal to the pure species molar volumes. We
     * have additionally specified in this class that the pure species molar
     * volumes are independent of temperature and pressure.
     */
    virtual void calcDensity();

public:
    virtual void setState_TP(doublereal T, doublereal pres);
    virtual void setState_TR(doublereal T, doublereal rho);
    virtual void setState_TPX(doublereal t, doublereal p, const doublereal* x);

protected:
    virtual void compositionChanged();
    void setMoleFractions_NoState(const doublereal* const x);

protected:
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
public:

    /// @name Thermodynamic Values for the Species Reference States
    /*!
     * There are also temporary variables for holding the species reference-
     * state values of Cp, H, S, and V at the last temperature and reference
     * pressure called. These functions are not recalculated if a new call is
     * made using the previous temperature. All calculations are done within the
     * routine _updateRefStateThermo().
     */
    //@{

    virtual void getEnthalpy_RT_ref(doublereal* hrt) const;
    virtual void getGibbs_RT_ref(doublereal* grt) const;

protected:
    //! Returns the vector of nondimensional Gibbs free energies of the
    //! reference state at the current temperature of the solution and the
    //! reference pressure for the species.
    /*!
     * @return  Output vector contains the nondimensional Gibbs free energies
     *          of the reference state of the species
     *          length = m_kk, units = dimensionless.
     */
    const vector_fp& gibbs_RT_ref() const;

public:
    virtual void getGibbs_ref(doublereal* g) const;
    virtual void getEntropy_R_ref(doublereal* er) const;
    virtual void getCp_R_ref(doublereal* cprt) const;
    virtual void getStandardVolumes_ref(doublereal* vol) const;

    //@}
    //! @name Initialization Methods - For Internal use
    /*!
     * The following methods are used in the process of constructing
     * the phase and setting its parameters from a specification in an
     * input file. They are not normally used in application programs.
     * To see how they are used, see importPhase().
     */
    //@{

    virtual bool addSpecies(shared_ptr<Species> spec);
    virtual void setStateFromXML(const XML_Node& state);

protected:
    //! @name Special Functions for fugacity classes
    //! @{

    //!  Calculate the value of z
    /*!
     *  \f[
     *        z = \frac{P v}{R T}
     *  \f]
     *
     *  returns the value of z
     */
    doublereal z() const;

    //! Calculate the deviation terms for the total entropy of the mixture from
    //! the ideal gas mixture
    /*
     * Here we use the current state conditions
     *
     * @returns the change in entropy in units of J kmol-1 K-1.
     */
    virtual doublereal sresid() const;

    //! Calculate the deviation terms for the total enthalpy of the mixture from
    //! the ideal gas mixture
    /*
     * Here we use the current state conditions
     *
     * @returns the change in entropy in units of J kmol-1.
     */
    virtual doublereal hresid() const;

    //! Estimate for the saturation pressure
    /*!
     * Note: this is only used as a starting guess for later routines that
     * actually calculate an accurate value for the saturation pressure.
     *
     * @param TKelvin  temperature in kelvin
     * @return the estimated saturation pressure at the given temperature
     */
    virtual doublereal psatEst(doublereal TKelvin) const;

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
    virtual doublereal liquidVolEst(doublereal TKelvin, doublereal& pres) const;

    //! Calculates the density given the temperature and the pressure and a
    //! guess at the density.
    /*!
     * Note, below T_c, this is a multivalued function. We do not cross the
     * vapor dome in this. This is protected because it is called during
     * setState_TP() routines. Infinite loops would result if it were not
     * protected.
     *
     *  -> why is this not const?
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
    virtual doublereal densityCalc(doublereal TKelvin, doublereal pressure, int phaseRequested,
                                   doublereal rhoguess);

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
    int corr0(doublereal TKelvin, doublereal pres, doublereal& densLiq,
              doublereal& densGas, doublereal& liqGRT, doublereal& gasGRT);

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
    virtual doublereal densSpinodalLiquid() const;

    //! Return the value of the density at the gas spinodal point (on the gas
    //! side) for the current temperature.
    /*!
     * @returns the density with units of kg m-3
     */
    virtual doublereal densSpinodalGas() const;

public:
    //! Calculate the saturation pressure at the current mixture content for the
    //! given temperature
    /*!
     * @param TKelvin         (input) Temperature (Kelvin)
     * @param molarVolGas     (return) Molar volume of the gas
     * @param molarVolLiquid  (return) Molar volume of the liquid
     * @returns the saturation pressure at the given temperature
     */
    doublereal calculatePsat(doublereal TKelvin, doublereal& molarVolGas,
                             doublereal& molarVolLiquid);

public:
    //! Calculate the saturation pressure at the current mixture content for the
    //! given temperature
    /*!
     * @param TKelvin   Temperature (Kelvin)
     * @return          The saturation pressure at the given temperature
     */
    virtual doublereal satPressure(doublereal TKelvin);

protected:
    //! Calculate the pressure given the temperature and the molar volume
    /*!
     * @param   TKelvin   temperature in kelvin
     * @param   molarVol  molar volume ( m3/kmol)
     * @returns the pressure.
     */
    virtual doublereal pressureCalc(doublereal TKelvin, doublereal molarVol) const;

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
    virtual doublereal dpdVCalc(doublereal TKelvin, doublereal molarVol, doublereal& presCalc) const;

    virtual void updateMixingExpressions();

    //@}

protected:
    virtual void invalidateCache();

    //! Storage for the current values of the mole fractions of the species
    /*!
     * This vector is kept up-to-date when some the setState functions are called.
     */
    vector_fp moleFractions_;

    //! Current state of the fluid
    /*!
     *  There are three possible states of the fluid:
     *  - FLUID_GAS
     *  - FLUID_LIQUID
     *  - FLUID_SUPERCRIT
     */
    int iState_;

    //! Force the system to be on a particular side of the spinodal curve
    int forcedState_;

    //! The last temperature at which the reference state thermodynamic
    //! properties were calculated at.
    mutable doublereal m_Tlast_ref;

    //! Temporary storage for dimensionless reference state enthalpies
    mutable vector_fp m_h0_RT;

    //! Temporary storage for dimensionless reference state heat capacities
    mutable vector_fp m_cp0_R;

    //! Temporary storage for dimensionless reference state Gibbs energies
    mutable vector_fp m_g0_RT;

    //! Temporary storage for dimensionless reference state entropies
    mutable vector_fp m_s0_R;
};
}

#endif
