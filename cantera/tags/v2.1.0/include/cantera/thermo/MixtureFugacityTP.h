/**
 *  @file MixtureFugacityTP.h
 *    Header file for a derived class of ThermoPhase that handles
 *    non-ideal mixtures based on the fugacity models (see \ref thermoprops and
 *    class \link Cantera::MixtureFugacityTP MixtureFugacityTP\endlink).
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef CT_MIXTUREFUGACITYTP_H
#define CT_MIXTUREFUGACITYTP_H

#include "ThermoPhase.h"
#include "VPSSMgr.h"
#include "cantera/numerics/ResidEval.h"

namespace Cantera
{

class XML_Node;
class PDSS;

//! Various states of the Fugacity object. In general there can be multiple liquid
//! objects for a single phase identified with each species.

#define FLUID_UNSTABLE -4
#define FLUID_UNDEFINED -3
#define FLUID_SUPERCRIT -2
#define FLUID_GAS       -1
#define FLUID_LIQUID_0   0
#define FLUID_LIQUID_1   1
#define FLUID_LIQUID_2   2
#define FLUID_LIQUID_3   3
#define FLUID_LIQUID_4   4
#define FLUID_LIQUID_5   5
#define FLUID_LIQUID_6   6
#define FLUID_LIQUID_7   7
#define FLUID_LIQUID_8   8
#define FLUID_LIQUID_9   9

/**
 * @ingroup thermoprops
 *
 *  This is a filter class for ThermoPhase that implements some preparatory
 *  steps for efficiently handling mixture of gases that whose standard states
 *  are defined as ideal gases, but which describe also non-ideal solutions.
 *  In addition a multicomponent liquid phase below the critical temperature of the
 *  mixture is also allowed. The main subclass is currently a mixture Redlich-Kwong class.
 *
 *  Several concepts are introduced. The first concept is there are temporary
 *  variables for holding the species standard state values
 *  of Cp, H, S, G, and V at the last temperature and pressure called. These functions are not recalculated
 *  if a new call is made using the previous temperature and pressure.
 *
 *  The other concept is that the current state of the mixture is tracked.
 *  The state variable is either GAS, LIQUID, or SUPERCRIT fluid.  Additionally,
 *  the variable LiquidContent is used and may vary between 0 and 1.
 *
 *  To support the above functionality, pressure and temperature variables,
 *  m_Plast_ss and m_Tlast_ss, are kept which store the last pressure and temperature
 *  used in the evaluation of standard state properties.
 *
 *  Typically, only one liquid phase is allowed to be formed within these classes.
 *  Additionally, there is an inherent contradiction between three phase models and
 *  the ThermoPhase class. The ThermoPhase class is really only meant to represent a
 *  single instantiation of a phase. The three phase models may be in equilibrium with
 *  multiple phases of the fluid in equilibrium with each other. This has yet to be resolved.
 *
 *  This class is usually used for non-ideal gases.
 *
 *  @nosubgrouping
 */
class MixtureFugacityTP : public ThermoPhase
{
public:
    //! @name Constructors and Duplicators for %MixtureFugacityTP
    //! @{

    //! Constructor.
    MixtureFugacityTP();

    //! Copy Constructor.
    /*!
     *  @param b   Object to be copied
     */
    MixtureFugacityTP(const MixtureFugacityTP& b);

    //! Assignment operator
    /*!
     *  @param b   Object to be copied
     */
    MixtureFugacityTP& operator=(const MixtureFugacityTP& b);

    //! Duplication routine
    /*!
     *  @return  Returns a duplication
     */
    virtual ThermoPhase* duplMyselfAsThermoPhase() const;

    //! @}
    //! @name  Utilities
    //! @{
    /**
     * Equation of state type flag. The base class returns
     * zero. Subclasses should define this to return a unique
     * non-zero value. Constants defined for this purpose are
     * listed in mix_defs.h.
     */
    virtual int eosType() const {
        return 0;
    }

    //! This method returns the convention used in specification
    //! of the standard state, of which there are currently two,
    //! temperature based, and variable pressure based.
    /*!
     * Currently, there are two standard state conventions:
     *  - Temperature-based activities,
     *    `cSS_CONVENTION_TEMPERATURE 0` (default)
     *  - Variable Pressure and Temperature based activities,
     *    `cSS_CONVENTION_VPSS 1`
     */
    virtual int standardStateConvention() const;

    //! Set the solution branch to force the ThermoPhase to exist on one branch or another
    /*!
     *  @param solnBranch  Branch that the solution is restricted to.
     *                     the value -1 means gas. The value -2 means unrestricted.
     *                     Values of zero or greater refer to species dominated condensed phases.
     */
    virtual void setForcedSolutionBranch(int solnBranch);

    //! Report the solution branch which the solution is restricted to
    /*!
     *  @return            Branch that the solution is restricted to.
     *                     the value -1 means gas. The value -2 means unrestricted.
     *                     Values of zero or greater refer to species dominated condensed phases.
     */
    virtual int forcedSolutionBranch() const;

    //! Report the solution branch which the solution is actually on
    /*!
     *  @return            Branch that the solution is restricted to.
     *                     the value -1 means gas. The value -2 means superfluid..
     *                     Values of zero or greater refer to species dominated condensed phases.
     */
    virtual int reportSolnBranchActual() const;

    //! Get the array of log concentration-like derivatives of the
    //! log activity coefficients
    /*!
     * For ideal mixtures (unity activity coefficients), this can return zero.
     * Implementations should take the derivative of the logarithm of the
     * activity coefficient with respect to the logarithm of the
     * concentration-like variable (i.e. moles) that represents the standard
     * state.
     *
     * This quantity is to be used in conjunction with derivatives of
     * that concentration-like variable when the derivative of the chemical
     * potential is taken.
     *
     *  units = dimensionless
     *
     * @param dlnActCoeffdlnN_diag    Output vector of derivatives of the
     *                         log Activity Coefficients. length = m_kk
     */
    virtual void getdlnActCoeffdlnN_diag(doublereal* dlnActCoeffdlnN_diag) const {
        err("getdlnActCoeffdlnN_diag");
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
     * We close the loop on this function, here, calling
     * getChemPotentials() and then dividing by RT. No need for child
     * classes to handle.
     *
     * @param mu    Output vector of  non-dimensional species chemical potentials
     *              Length: m_kk.
     */
    void getChemPotentials_RT(doublereal* mu) const;

    //@}
    /*!
     * @name  Properties of the Standard State of the Species in the Solution
     *
     *  Within MixtureFugacityTP, these properties are calculated via a common routine,
     *  _updateStandardStateThermo(),
     *  which must be overloaded in inherited objects.
     *  The values are cached within this object, and are not recalculated unless
     *  the temperature or pressure changes.
     */
    //@{

    //!   Get the array of chemical potentials at unit activity.
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

    //! Get the nondimensional Enthalpy functions for the species
    //! at their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
    * For all objects with the Mixture Fugacity approximation, we define the
    * standard state as an ideal gas at the current temperature and pressure
    * of the solution.
    *
    * @param hrt     Output vector of standard state enthalpies.
    *                length = m_kk. units are unitless.
    */
    virtual void getEnthalpy_RT(doublereal* hrt) const;

    //! Get the array of nondimensional Enthalpy functions for the standard state species
    /*!
     * at the current <I>T</I> and <I>P</I> of the solution.
     * For all objects with the Mixture Fugacity approximation, we define the
     * standard state as an ideal gas at the current temperature and pressure
     * of the solution.
     *
     * @param sr     Output vector of nondimensional standard state
     *               entropies. length = m_kk.
     */
    virtual void getEntropy_R(doublereal* sr) const;

    //! Get the nondimensional Gibbs functions for the species
    //! at their standard states of solution at the current T and P of the solution.
    /*!
     * For all objects with the Mixture Fugacity approximation, we define the
     * standard state as an ideal gas at the current temperature and pressure
     * of the solution.
     *
     * @param grt    Output vector of nondimensional standard state
     *               Gibbs free energies. length = m_kk.
     */
    virtual void getGibbs_RT(doublereal* grt) const;

    //! Get the pure Gibbs free energies of each species.
    //! Species are assumed to be in their standard states. This is the same
    //! as getStandardChemPotentials().
    //! @param[out] gpure Array of standard state Gibbs free energies.
    //!        length = m_kk. units are J/kmol.
    void getPureGibbs(doublereal* gpure) const;

    //!  Returns the vector of nondimensional internal Energies of the standard state at the current temperature
    //!  and pressure of the solution for each species.
    /*!
     * For all objects with the Mixture Fugacity approximation, we define the
     * standard state as an ideal gas at the current temperature and pressure
     * of the solution.
     *
     * \f[
     *  u^{ss}_k(T,P) = h^{ss}_k(T)  - P * V^{ss}_k
     * \f]
     *
     * @param urt    Output vector of nondimensional standard state
     *               internal energies. length = m_kk.
     */
    virtual void getIntEnergy_RT(doublereal* urt) const;

    //! Get the nondimensional Heat Capacities at constant
    //! pressure for the standard state of the species  at the current T and P.
    /*!
     * For all objects with the Mixture Fugacity approximation, we define the
     * standard state as an ideal gas at the current temperature and pressure of the solution.
     *
     * @param cpr    Output vector containing the
     *               the nondimensional Heat Capacities at constant
     *               pressure for the standard state of the species.
     *               Length: m_kk.
     */
    virtual void getCp_R(doublereal* cpr) const;

    //! Get the molar volumes of each species in their standard
    //! states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * For all objects with the Mixture Fugacity approximation, we define the
     * standard state as an ideal gas at the current temperature and pressure of the solution.
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
     *    Currently this passes down to setState_TP(). It does not
     *    make sense to calculate the standard state without first
     *    setting T and P.
     *
     * @param temp  Temperature (kelvin)
     */
    virtual void setTemperature(const doublereal temp);

    //! Set the internally stored pressure (Pa) at constant
    //! temperature and composition
    /*!
     *  Currently this passes down to setState_TP().  It does not
     *    make sense to calculate the standard state without first
     *    setting T and P.
     *
     *  @param p input Pressure (Pa)
     */
    virtual void setPressure(doublereal p);

protected:
    /**
     * Calculate the density of the mixture using the partial
     * molar volumes and mole fractions as input
     *
     * The formula for this is
     *
     * \f[
     * \rho = \frac{\sum_k{X_k W_k}}{\sum_k{X_k V_k}}
     * \f]
     *
     * where \f$X_k\f$ are the mole fractions, \f$W_k\f$ are
     * the molecular weights, and \f$V_k\f$ are the pure species
     * molar volumes.
     *
     * Note, the basis behind this formula is that in an ideal
     * solution the partial molar volumes are equal to the pure
     * species molar volumes. We have additionally specified
     * in this class that the pure species molar volumes are
     * independent of temperature and pressure.
     */
    virtual void calcDensity();

public:
    //! Set the temperature and pressure at the same time
    /*!
     *  Note this function triggers a reevaluation of the standard
     *  state quantities.
     *
     *  @param T  temperature (kelvin)
     *  @param pres pressure (pascal)
     */
    virtual void setState_TP(doublereal T, doublereal pres);

    //! Set the internally stored temperature (K) and density (kg/m^3)
    /*!
     * @param T     Temperature in kelvin
     * @param rho   Density (kg/m^3)
     */
    virtual void setState_TR(doublereal T, doublereal rho);

    //! Set the temperature (K), pressure (Pa), and mole fractions.
    /*!
     * Note, the mole fractions are set first before the pressure is set.
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param t    Temperature (K)
     * @param p    Pressure (Pa)
     * @param x    Vector of mole fractions. Length is equal to m_kk.
     */
    virtual void setState_TPX(doublereal t, doublereal p, const doublereal* x);

    //! Set the mass fractions to the specified values, and then
    //! normalize them so that they sum to 1.0.
    /*!
     * @param y Array of unnormalized mass fraction values (input).
     *          Must have a length greater than or equal to the number of species.
     */
    virtual void setMassFractions(const doublereal* const y);

    //!Set the mass fractions to the specified values without normalizing.
    /*!
     * This is useful when the normalization
     * condition is being handled by some other means, for example
     * by a constraint equation as part of a larger set of
     * equations.
     *
     * @param y  Input vector of mass fractions. Length is m_kk.
     */
    virtual void setMassFractions_NoNorm(const doublereal* const y);

    //! Set the mole fractions to the specified values, and then
    //! normalize them so that they sum to 1.0.
    /*!
     * @param x Array of unnormalized mole fraction values (input).
     *          Must have a length greater than or equal to the number of species.
     */
    virtual void setMoleFractions(const doublereal* const x);

    //! Set the mole fractions to the specified values without normalizing.
    /*!
     * This is useful when the normalization
     * condition is being handled by some other means, for example
     * by a constraint equation as part of a larger set of equations.
     *
     * @param x  Input vector of mole fractions. Length is m_kk.
     */
    virtual void setMoleFractions_NoNorm(const doublereal* const x);

    //! Set the concentrations to the specified values within the phase.
    /*!
     * @param c The input vector to this routine is in dimensional
     *        units. For volumetric phases c[k] is the
     *        concentration of the kth species in kmol/m3.
     *        For surface phases, c[k] is the concentration
     *        in kmol/m2. The length of the vector is the number
     *        of species in the phase.
     */
    virtual void setConcentrations(const doublereal* const c);

protected:
    void setMoleFractions_NoState(const doublereal* const x);

public:
    //! Returns the current pressure of the phase
    /*!
     *  The pressure is an independent variable in this phase. Its current value
     *  is stored in the object MixtureFugacityTP.
     *
     * @return return the pressure in pascals.
     */
    doublereal pressure() const {
        return m_Pcurrent;
    }

protected:
    //! Updates the reference state thermodynamic functions at the current T  of the solution.
    /*!
     * This function must be called for every call to functions in this
     * class. It checks to see whether the temperature has changed and
     * thus the ss thermodynamics functions for all of the species
     * must be recalculated.
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

    /// @name Thermodynamic Values for the Species Reference States (MixtureFugacityTP)
    /*!
     *  There are also temporary
     *  variables for holding the species reference-state values of Cp, H, S, and V at the
     *  last temperature and reference pressure called. These functions are not recalculated
     *  if a new call is made using the previous temperature.
     *  All calculations are done within the routine  _updateRefStateThermo().
     */
    //@{

    //!  Returns the vector of nondimensional
    //!  enthalpies of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    /*!
     * @param hrt Output vector contains the nondimensional enthalpies
     *            of the reference state of the species
     *            length = m_kk, units = dimensionless.
     */
    virtual void getEnthalpy_RT_ref(doublereal* hrt) const;

#ifdef H298MODIFY_CAPABILITY
    //!  Modify the value of the 298 K Heat of Formation of the standard state of
    //!  one species in the phase (J kmol-1)
    /*!
     *   The 298K heat of formation is defined as the enthalpy change to create the standard state
     *   of the species from its constituent elements in their standard states at 298 K and 1 bar.
     *
     *   @param  k           Index of the species
     *   @param  Hf298New    Specify the new value of the Heat of Formation at 298K and 1 bar.
     *                       units = J/kmol.
     */
    void modifyOneHf298SS(const int k, const doublereal Hf298New);
#endif

    //!  Returns the vector of nondimensional
    //!  Gibbs free energies of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    /*!
     * @param grt Output vector contains the nondimensional Gibbs free energies
     *            of the reference state of the species
     *            length = m_kk, units = dimensionless.
     */
    virtual void getGibbs_RT_ref(doublereal* grt) const;

protected:
    //!  Returns the vector of nondimensional
    //!  Gibbs free energies of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    /*!
     * @return  Output vector contains the nondimensional Gibbs free energies
     *          of the reference state of the species
     *          length = m_kk, units = dimensionless.
     */
    const vector_fp& gibbs_RT_ref() const;

public:
    /*!
     *  Returns the vector of the
     *  gibbs function of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *  units = J/kmol
     *
     * @param g   Output vector contain the Gibbs free energies
     *            of the reference state of the species
     *            length = m_kk, units = J/kmol.
     */
    virtual void getGibbs_ref(doublereal* g) const;

    /*!
     *  Returns the vector of nondimensional
     *  entropies of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *
     * @param er  Output vector contain the nondimensional entropies
     *            of the species in their reference states
     *            length: m_kk, units: dimensionless.
     */
    virtual void getEntropy_R_ref(doublereal* er) const;

    /*!
     *  Returns the vector of nondimensional
     *  constant pressure heat capacities of the reference state
     *  at the current temperature of the solution
     *  and reference pressure for the species.
     *
     * @param cprt Output vector contains the nondimensional heat capacities
     *             of the species in their reference states
     *             length: m_kk, units: dimensionless.
     */
    virtual void getCp_R_ref(doublereal* cprt) const;

    //!  Get the molar volumes of the species reference states at the current
    //!  <I>T</I> and reference pressure of the solution.
    /*!
     * units = m^3 / kmol
     *
     * @param vol     Output vector containing the standard state volumes.
     *                Length: m_kk.
     */
    virtual void getStandardVolumes_ref(doublereal* vol) const;

    //@}
    //! @name Initialization Methods - For Internal use
    /*!
     * The following methods are used in the process of constructing
     * the phase and setting its parameters from a specification in an
     * input file. They are not normally used in application programs.
     * To see how they are used, see files importCTML.cpp and
     * ThermoFactory.cpp.
     */
    //@{

    //! Set the initial state of the phase to the conditions  specified in the state XML element.
    /*!
     * This method sets the temperature, pressure, and mole  fraction vector to a set default value.
     *
     * @param state An XML_Node object corresponding to
     *              the "state" entry for this phase in the input file.
     */
    virtual void setStateFromXML(const XML_Node& state);

    //! @internal Initialize the object
    /*!
     * This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species.  The base class implementation does nothing,
     * and subclasses that do not require initialization do not
     * need to overload this method.  When importing a CTML phase
     * description, this method is called after calling installSpecies()
     * for each species in the phase. It's called before calling
     * initThermoXML() for the phase. Therefore, it's the correct
     * place for initializing vectors which have lengths equal to the
     * number of species.
     *
     * @see importCTML.cpp
     */
    virtual void initThermo();

    //!   Initialize a ThermoPhase object, potentially reading activity
    //!   coefficient information from an XML database.
    /*!
     * This routine initializes the lengths in the current object and
     * then calls the parent routine.
     * This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species.  The base class implementation does nothing,
     * and subclasses that do not require initialization do not
     * need to overload this method.  When importing a CTML phase
     * description, this method is called just prior to returning
     * from function importPhase().
     *
     * @param phaseNode This object must be the phase node of a
     *             complete XML tree
     *             description of the phase, including all of the
     *             species data. In other words while "phase" must
     *             point to an XML phase object, it must have
     *             sibling nodes "speciesData" that describe
     *             the species in the phase.
     * @param id   ID of the phase. If nonnull, a check is done
     *             to see if phaseNode is pointing to the phase
     *             with the correct id.
     */
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

private:
    //!  @internal Initialize the internal lengths in this object.
    void initLengths();
    //@}

protected:
    //! @name Special Functions for fugacity classes
    //! @{

    //!  Calculate the value of z
    /*!
     *  \f[
     *        z = \frac{P v}{  R  T}
     *  \f]
     *
     *  returns the value of z
     */
    doublereal z() const;

    //! Calculate the deviation terms for the total entropy of the mixture from the
    //! ideal gas mixture
    /*
     *  Here we use the current state conditions
     *
     * @return  Returns the change in entropy in units of J kmol-1 K-1.
     */
    virtual doublereal sresid() const;

    //! Calculate the deviation terms for the total enthalpy of the mixture from the ideal gas mixture
    /*
     *  Here we use the current state conditions
     *
     * @return  Returns the change in entropy in units of J kmol-1.
     */
    virtual doublereal hresid() const;

    //! Estimate for the saturation pressure
    /*!
     *  Note: this is only used as a starting guess for later routines that actually calculate an
     *  accurate value for the saturation pressure.
     *
     *  @param TKelvin  temperature in kelvin
     *
     *  @return returns the estimated saturation pressure at the given temperature
     */
    virtual doublereal psatEst(doublereal TKelvin) const;

public:
    //! Estimate for the molar volume of the liquid
    /*!
     *   Note: this is only used as a starting guess for later routines that actually calculate an
     *  accurate value for the liquid molar volume.
     *  This routine doesn't change the state of the system.
     *
     *  @param TKelvin  temperature in kelvin
     *  @param pres     Pressure in Pa. This is used as an initial guess. If the routine
     *                  needs to change the pressure to find a stable liquid state, the
     *                  new pressure is returned in this variable.
     *
     *  @return Returns the estimate of the liquid volume.  If the liquid can't be found, this
     *          routine returns -1.
     */
    virtual doublereal liquidVolEst(doublereal TKelvin, doublereal& pres) const;

    //!  Calculates the density given the temperature and the pressure and a guess at the density.
    /*!
     * Note, below T_c, this is a multivalued function. We do not cross the vapor dome in this.
     * This is protected because it is called during setState_TP() routines. Infinite loops would result
     * if it were not protected.
     *
     *  -> why is this not const?
     *
     * parameters:
     *    @param TKelvin   Temperature in Kelvin
     *    @param pressure  Pressure in Pascals (Newton/m**2)
     *    @param phaseRequested     int representing the phase whose density we are requesting. If we put
     *                     a gas or liquid phase here, we will attempt to find a volume in that
     *                     part of the volume space, only, in this routine. A value of FLUID_UNDEFINED
     *                     means that we will accept anything.
     *
     *   @param rhoguess   Guessed density of the fluid. A value of -1.0 indicates that there
     *                     is no guessed density
     *
     *  @return   We return the density of the fluid at the requested phase. If we have not found any
     *            acceptable density we return a -1. If we have found an acceptable density at a
     *            different phase, we return a -2.
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
     * @param checkState If true, this function does a complete check to see where
     *        in parameters space we are
     *
     *  There are three values:
     *  - WATER_GAS   below the critical temperature but below the critical density
     *  - WATER_LIQUID  below the critical temperature but above the critical density
     *  - WATER_SUPERCRIT   above the critical temperature
     */
    int phaseState(bool checkState = false) const ;

    //! Return the value of the density at the liquid spinodal point (on the liquid side)
    //! for the current temperature.
    /*!
     * @return returns the density with units of kg m-3
     */
    virtual doublereal densSpinodalLiquid() const;

    //! Return the value of the density at the gas spinodal point (on the gas side)
    //! for the current temperature.
    /*!
     * @return returns the density with units of kg m-3
     */
    virtual doublereal densSpinodalGas() const;

public:
    //! Calculate the saturation pressure at the current mixture content for the given temperature
    /*!
     *   @param TKelvin         (input) Temperature (Kelvin)
     *   @param molarVolGas     (return) Molar volume of the gas
     *   @param molarVolLiquid  (return) Molar volume of the liquid
     *
     *   @return          Returns the saturation pressure at the given temperature
     */
    doublereal calculatePsat(doublereal TKelvin, doublereal& molarVolGas,
                             doublereal& molarVolLiquid);
    
public:
    //! Calculate the saturation pressure at the current mixture content for the given temperature
    /*!
     *   @param TKelvin         (input) Temperature (Kelvin)
     *   @param molarVolGas     (return) Molar volume of the gas
     *   @param molarVolLiquid  (return) Molar volume of the liquid
     *
     *   @return          Returns the saturation pressure at the given temperature
     */
    virtual doublereal satPressure(doublereal TKelvin);

protected:
    //! Calculate the pressure given the temperature and the molar volume
    /*!
     *  Calculate the pressure given the temperature and the molar volume
     *
     * @param   TKelvin   temperature in kelvin
     * @param   molarVol  molar volume ( m3/kmol)
     *
     * @return  Returns the pressure.
     */
    virtual doublereal pressureCalc(doublereal TKelvin, doublereal molarVol) const;

    //! Calculate the pressure and the pressure derivative given the temperature and the molar volume
    /*!
     *  Temperature and mole number are held constant
     *
     * @param   TKelvin   temperature in kelvin
     * @param   molarVol  molar volume ( m3/kmol)
     *
     * @param   presCalc  Returns the pressure.
     *
     *  @return  Returns the derivative of the pressure wrt the molar volume
     */
    virtual doublereal dpdVCalc(doublereal TKelvin, doublereal molarVol, doublereal& presCalc) const;

    virtual void updateMixingExpressions();

    //@}

    class spinodalFunc : public Cantera::ResidEval
    {
    public:
        spinodalFunc(MixtureFugacityTP* tp);
        virtual int evalSS(const doublereal t, const doublereal* const y, doublereal* const r);
        MixtureFugacityTP* m_tp;
    };

protected:
    //! Current value of the pressures
    /*!
     *  Because the pressure is now a calculation, we store the result of the calculation whenever
     *  it is recalculated.
     *
     *  units = Pascals
     */
    doublereal  m_Pcurrent;

    //! Storage for the current values of the mole fractions of the species
    /*!
     * This vector is kept up-to-date when some the setState functions are called.
     */
    std::vector<doublereal> moleFractions_;

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

    //! The last temperature at which the reference state thermodynamic properties were calculated at.
    mutable doublereal m_Tlast_ref;

    //! Temporary storage for log of p/rt
    mutable doublereal    m_logc0;

    //! Temporary storage for dimensionless reference state enthalpies
    mutable vector_fp      m_h0_RT;

    //! Temporary storage for dimensionless reference state heat capacities
    mutable vector_fp      m_cp0_R;

    //! Temporary storage for dimensionless reference state gibbs energies
    mutable vector_fp      m_g0_RT;

    //! Temporary storage for dimensionless reference state entropies
    mutable vector_fp      m_s0_R;

    spinodalFunc* fdpdv_;

private:
    //! MixtureFugacityTP has its own err routine
    /*!
     * @param msg  Error message string
     */
    doublereal err(const std::string& msg) const;
};
}

#endif
