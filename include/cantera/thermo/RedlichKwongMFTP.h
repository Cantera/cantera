/**
 *  @file RedlichKwongMFTP.h
 * Definition file for a derived class of ThermoPhase that assumes either
 * an ideal gas or ideal solution approximation and handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties (see \ref thermoprops and
 * class \link Cantera::RedlichKwongMFTP RedlichKwongMFTP\endlink).
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef CT_REDLICHKWONGMFTP_H
#define CT_REDLICHKWONGMFTP_H

#include "MixtureFugacityTP.h"
#include "cantera/base/Array.h"

namespace Cantera
{
/**
 * @ingroup thermoprops
 *
 *  This class can handle either an ideal solution or an ideal gas approximation
 *  of a phase.
 */
class RedlichKwongMFTP : public MixtureFugacityTP
{
public:
    //! @name Constructors and Duplicators
    //! @{

    //! Base constructor.
    RedlichKwongMFTP();

    //! Construct and initialize a RedlichKwongMFTP object directly from an
    //! ASCII input file
    /*!
     * @param infile    Name of the input file containing the phase XML data
     *                  to set up the object
     * @param id        ID of the phase in the input file. Defaults to the empty string.
     */
    RedlichKwongMFTP(const std::string& infile, std::string id="");

    //! Construct and initialize a RedlichKwongMFTP object directly from an
    //! XML database
    /*!
     *  @param phaseRef XML phase node containing the description of the phase
     *  @param id       id attribute containing the name of the phase.  (default is the empty string)
     */
    RedlichKwongMFTP(XML_Node& phaseRef, const std::string& id = "");

    //!  This is a special constructor, used to replicate test problems
    //!  during the initial verification of the object
    /*!
     *  test problems:
     *    1:     Pure CO2 problem
     *           input file = CO2_RedlickKwongMFTP.xml
     *
     * @param testProb Hard -coded test problem to instantiate.
     *                 Current valid values are 1.
     *  @deprecated To be removed after Cantera 2.2.
     */
    RedlichKwongMFTP(int testProb);

    //! Copy Constructor
    /*!
     * Copy constructor for the object. Constructed object will be a clone of this object, but will
     * also own all of its data. This is a wrapper around the assignment operator
     *
     * @param right Object to be copied.
     */
    RedlichKwongMFTP(const RedlichKwongMFTP& right);

    //! Assignment operator
    /*!
     * Assignment operator for the object. Constructed object will be a clone of this object, but will
     * also own all of its data.
     *
     * @param right Object to be copied.
     */
    RedlichKwongMFTP& operator=(const RedlichKwongMFTP& right);

    //! Duplicator from the ThermoPhase parent class
    /*!
     * Given a pointer to a ThermoPhase object, this function will
     * duplicate the ThermoPhase object and all underlying structures.
     * This is basically a wrapper around the copy constructor.
     *
     * @return returns a pointer to a ThermoPhase
     */
    virtual ThermoPhase* duplMyselfAsThermoPhase() const;

    /**
     * Equation of state type flag. The base class returns
     * zero. Subclasses should define this to return a unique
     * non-zero value. Constants defined for this purpose are
     * listed in mix_defs.h.
     */
    virtual int eosType() const;

    //! @name Molar Thermodynamic properties
    //! @{

    /// Molar enthalpy. Units: J/kmol.
    virtual doublereal enthalpy_mole() const;

    /// Molar entropy. Units: J/kmol/K.
    virtual doublereal entropy_mole() const;

    /// Molar heat capacity at constant pressure. Units: J/kmol/K.
    virtual doublereal cp_mole() const;

    /// Molar heat capacity at constant volume. Units: J/kmol/K.
    virtual doublereal cv_mole() const;

    //! @}
    //! @name Mechanical Properties
    //! @{

    //! Return the thermodynamic pressure (Pa).
    /*!
     *  Since the mass density, temperature, and mass fractions are stored,
     *  this method uses these values to implement the
     *  mechanical equation of state \f$ P(T, \rho, Y_1, \dots,  Y_K) \f$.
     *
     * \f[
     *    P = \frac{RT}{v-b_{mix}} - \frac{a_{mix}}{T^{0.5} v \left( v + b_{mix} \right) }
     * \f]
     */
    virtual doublereal pressure() const;

    // @}

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
     * solution the partial molar volumes are equal to the
     * species standard state molar volumes.
     * The species molar volumes may be functions
     * of temperature and pressure.
     */
    virtual void calcDensity();

protected:
    //! Set the temperature (K)
    /*!
     * This function sets the temperature, and makes sure that
     * the value propagates to underlying objects
     *
     * @param temp Temperature in kelvin
     */
    virtual void setTemperature(const doublereal temp);

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

public:
    //! This method returns an array of generalized concentrations
    /*!
     * \f$ C^a_k\f$ are defined such that \f$ a_k = C^a_k /
     * C^0_k, \f$ where \f$ C^0_k \f$ is a standard concentration
     * defined below and \f$ a_k \f$ are activities used in the
     * thermodynamic functions.  These activity (or generalized)
     * concentrations are used
     * by kinetics manager classes to compute the forward and
     * reverse rates of elementary reactions. Note that they may
     * or may not have units of concentration --- they might be
     * partial pressures, mole fractions, or surface coverages,
     * for example.
     *
     * @param c Output array of generalized concentrations. The
     *           units depend upon the implementation of the
     *           reaction rate expressions within the phase.
     */
    virtual void getActivityConcentrations(doublereal* c) const;

    //! Returns the standard concentration \f$ C^0_k \f$, which is used to normalize
    //! the generalized concentration.
    /*!
     * This is defined as the concentration  by which the generalized
     * concentration is normalized to produce the activity.
     * In many cases, this quantity will be the same for all species in a phase.
     * Since the activity for an ideal gas mixture is
     * simply the mole fraction, for an ideal gas \f$ C^0_k = P/\hat R T \f$.
     *
     * @param k Optional parameter indicating the species. The default
     *          is to assume this refers to species 0.
     * @return
     *   Returns the standard Concentration in units of m3 kmol-1.
     */
    virtual doublereal standardConcentration(size_t k=0) const;

    //! Returns the units of the standard and generalized concentrations.
    /*!
     * Note they have the same units, as their
     * ratio is defined to be equal to the activity of the kth
     * species in the solution, which is unitless.
     *
     * This routine is used in print out applications where the
     * units are needed. Usually, MKS units are assumed throughout
     * the program and in the XML input files.
     *
     * The base ThermoPhase class assigns the default quantities
     * of (kmol/m3) for all species.
     * Inherited classes are responsible for overriding the default
     * values if necessary.
     *
     * @param uA Output vector containing the units:
     *
     *     uA[0] = kmol units - default  = 1
     *     uA[1] = m    units - default  = -nDim(), the number of spatial
     *                                   dimensions in the Phase class.
     *     uA[2] = kg   units - default  = 0;
     *     uA[3] = Pa(pressure) units - default = 0;
     *     uA[4] = Temperature units - default = 0;
     *     uA[5] = time units - default = 0
     *
     * @param k species index. Defaults to 0.
     * @param sizeUA output int containing the size of the vector.
     *        Currently, this is equal to 6.
     * @deprecated To be removed after Cantera 2.2.
     */
    virtual void getUnitsStandardConc(double* uA, int k = 0, int sizeUA = 6) const;

    //! Get the array of non-dimensional activity coefficients at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     * For all objects with the Mixture Fugacity approximation, we define the
     * standard state as an ideal gas at the current temperature and pressure
     * of the solution. The activities are based on this standard state.
     *
     * @param ac Output vector of activity coefficients. Length: m_kk.
     */
    virtual void getActivityCoefficients(doublereal* ac) const;

    /// @name  Partial Molar Properties of the Solution
    //@{

    //! Get the array of non-dimensional species chemical potentials.
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

    //! Get the species chemical potentials. Units: J/kmol.
    /*!
     * This function returns a vector of chemical potentials of the
     * species in solution at the current temperature, pressure
     * and mole fraction of the solution.
     *
     * @param mu  Output vector of species chemical
     *            potentials. Length: m_kk. Units: J/kmol
     */
    virtual void getChemPotentials(doublereal* mu) const;

    //!  Get the species partial molar enthalpies. Units: J/kmol.
    /*!
     * @param hbar    Output vector of species partial molar enthalpies.
     *                Length: m_kk. units are J/kmol.
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;

    //! Get the species partial molar entropies. Units: J/kmol/K.
    /*!
     * @param sbar    Output vector of species partial molar entropies.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;

    //! Get the species partial molar enthalpies. Units: J/kmol.
    /*!
     * @param ubar    Output vector of species partial molar internal energies.
     *                Length = m_kk. units are J/kmol.
     */
    virtual void getPartialMolarIntEnergies(doublereal* ubar) const;

    //! Get the partial molar heat capacities Units: J/kmol/K
    /*!
     * @param cpbar   Output vector of species partial molar heat capacities
     *                at constant pressure.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const;

    //! Get the species partial molar volumes. Units: m^3/kmol.
    /*!
     *  @param vbar   Output vector of species partial molar volumes.
     *                Length = m_kk. units are m^3/kmol.
     */
    virtual void getPartialMolarVolumes(doublereal* vbar) const;

    //@}
    /// @name Critical State Properties.
    /// These methods are only implemented by some subclasses, and may
    /// be moved out of ThermoPhase at a later date.
    //@{

    /// Critical temperature (K).
    virtual doublereal critTemperature() const;

    /// Critical pressure (Pa).
    virtual doublereal critPressure() const;

    /// Critical volume (m3/kmol)
    virtual doublereal critVolume() const;

    // Critical compressibility (unitless)
    virtual doublereal critCompressibility() const;

    /// Critical density (kg/m3).
    virtual doublereal critDensity() const;

public:

    //@}
    //! @name Initialization Methods - For Internal use
    /*!
     * The following methods are used in the process of constructing
     * the phase and setting its parameters from a specification in an
     * input file. They are not normally used in application programs.
     * To see how they are used, see importPhase().
     */
    //@{

    //! Set equation of state parameter values from XML entries.
    /*!
     * This method is called by function importPhase() when processing a phase
     * definition in an input file. It should be overloaded in subclasses to set
     * any parameters that are specific to that particular phase model.
     *
     * @param thermoNode An XML_Node object corresponding to
     *                   the "thermo" entry for this phase in the input file.
     */
    virtual void setParametersFromXML(const XML_Node& thermoNode);

    //! @internal Initialize the object
    /*!
     * This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species.  The base class implementation does nothing,
     * and subclasses that do not require initialization do not
     * need to overload this method.  When importing a CTML phase
     * description, this method is called just prior to returning
     * from function importPhase().
     */
    virtual void initThermo();

    //!This method is used by the ChemEquil equilibrium solver.
    /*!
     * It sets the state such that the chemical potentials satisfy
     * \f[ \frac{\mu_k}{\hat R T} = \sum_m A_{k,m}
     * \left(\frac{\lambda_m} {\hat R T}\right) \f] where
     * \f$ \lambda_m \f$ is the element potential of element m. The
     * temperature is unchanged.  Any phase (ideal or not) that
     * implements this method can be equilibrated by ChemEquil.
     *
     * @param lambda_RT Input vector of dimensionless element potentials
     *                  The length is equal to nElements().
     */
    void setToEquilState(const doublereal* lambda_RT);

    //!   Initialize a ThermoPhase object, potentially reading activity
    //!   coefficient information from an XML database.
    /*!
     *
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
    //! Read the pure species RedlichKwong input parameters
    /*!
     *  @param pureFluidParam   XML_Node for the pure fluid parameters
     */
    void readXMLPureFluid(XML_Node& pureFluidParam);

    //! Apply mixing rules for a coefficients
    void applyStandardMixingRules();

    //! Read the cross species RedlichKwong input parameters
    /*!
     *  @param pureFluidParam   XML_Node for the cross fluid parameters
     */
    void readXMLCrossFluid(XML_Node& pureFluidParam);

private:
    //!  @internal Initialize the internal lengths in this object.
    /*!
     * Note this is not a virtual function and only handles this object
     */
    void initLengths();
    // @}

protected:
    // Special functions inherited from MixtureFugacityTP

    //! Calculate the deviation terms for the total entropy of the mixture from the
    //! ideal gas mixture
    /*!
     *  Here we use the current state conditions
     *
     * @return  Returns the change in entropy in units of J kmol-1 K-1.
     */
    virtual doublereal sresid() const;

    // Calculate the deviation terms for the total enthalpy of the mixture from the
    // ideal gas mixture
    /*
     *  Here we use the current state conditions
     *
     * @return  Returns the change in enthalpy in units of J kmol-1.
     */
    virtual doublereal hresid() const;
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
     *  @return Returns the estimate of the liquid volume.
     */
    virtual doublereal liquidVolEst(doublereal TKelvin, doublereal& pres) const;

public:
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
     *    @param phase     int representing the phase whose density we are requesting. If we put
     *                     a gas or liquid phase here, we will attempt to find a volume in that
     *                     part of the volume space, only, in this routine. A value of FLUID_UNDEFINED
     *                     means that we will accept anything.
     *
     *   @param rhoguess   Guessed density of the fluid. A value of -1.0 indicates that there
     *                     is no guessed density
     *
     *
     *  @return   We return the density of the fluid at the requested phase. If we have not found any
     *            acceptable density we return a -1. If we have found an acceptable density at a
     *            different phase, we return a -2.
     */
    virtual doublereal densityCalc(doublereal TKelvin, doublereal pressure, int phase, doublereal rhoguess);

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

    //! Calculate dpdV and dpdT at the current conditions
    /*!
     *  These are stored internally.
     */
    void pressureDerivatives() const;

    virtual void updateMixingExpressions();

    //! Update the a and b parameters
    /*!
     *  The a and the b parameters depend on the mole fraction and the temperature.
     *  This function updates the internal numbers based on the state of the object.
     */
    void updateAB();

    //!  Calculate the a and the b parameters given the temperature
    /*!
     *  This function doesn't change the internal state of the object, so it is a const
     *  function.  It does use the stored mole fractions in the object.
     *
     *  @param temp  Temperature (TKelvin)
     *
     *  @param aCalc (output)  Returns the a value
     *  @param bCalc (output)  Returns the b value.
     */
    void calculateAB(doublereal temp, doublereal& aCalc, doublereal& bCalc) const;

    // Special functions not inherited from MixtureFugacityTP

    doublereal da_dt() const;

    void calcCriticalConditions(doublereal a, doublereal b, doublereal a0_coeff, doublereal aT_coeff,
                                doublereal& pc, doublereal& tc, doublereal& vc) const;

    //! Solve the cubic equation of state
    /*!
     * The R-K equation of state may be solved via the following formula:
     *
     *     V**3 - V**2(RT/P)  - V(RTb/P - a/(P T**.5) + b*b) - (a b / (P T**.5)) = 0
     *
     * Returns the number of solutions found. If it only finds the liquid
     * branch solution, it will return a -1 or a -2 instead of 1 or 2.  If it
     * returns 0, then there is an error.
     */
    int NicholsSolve(double TKelvin, double pres, doublereal a, doublereal b,
                     doublereal Vroot[3]) const;

protected:
    //! boolean indicating whether standard mixing rules are applied
    /*!
     *  - 1 = Yes, there are standard cross terms in the a coefficient matrices.
     *  - 0 = No, there are nonstandard cross terms in the a coefficient matrices.
     */
    int m_standardMixingRules;

    //! Form of the temperature parameterization
    /*!
     *  - 0 = There is no temperature parameterization of a or b
     *  - 1 = The a_ij parameter is a linear function of the temperature
     */
    int m_formTempParam;

    //! Value of b in the equation of state
    /*!
     *  m_b is a function of the temperature and the mole fraction.
     */
    doublereal m_b_current;

    //! Value of a in the equation of state
    /*!
     *  a_b is a function of the temperature and the mole fraction.
     */
    doublereal m_a_current;

    vector_fp a_vec_Curr_;
    vector_fp b_vec_Curr_;

    Array2D  a_coeff_vec;

    vector_fp m_pc_Species;
    vector_fp m_tc_Species;
    vector_fp m_vc_Species;

    int NSolns_;

    doublereal Vroot_[3];

    //! Temporary storage - length = m_kk.
    mutable vector_fp m_pp;

    //! Temporary storage - length = m_kk.
    mutable vector_fp m_tmpV;

    // mutable vector_fp m_tmpV2;

    // Partial molar volumes of the species
    mutable vector_fp m_partialMolarVolumes;

    //! The derivative of the pressure wrt the volume
    /*!
     *  Calculated at the current conditions
     *  temperature and mole number kept constant
     */
    mutable doublereal dpdV_;

    //! The derivative of the pressure wrt the temperature
    /*!
     *  Calculated at the current conditions
     *  Total volume and mole number kept constant
     */
    mutable doublereal dpdT_;

    //! Vector of derivatives of pressure wrt mole number
    /*!
     *  Calculated at the current conditions
     *  Total volume, temperature and other mole number kept constant
     */
    mutable vector_fp dpdni_;

public:
    //! Omega constant for a -> value of a in terms of critical properties
    /*!
     *  this was calculated from a small nonlinear solve
     */
    static const doublereal omega_a;

    //! Omega constant for b
    static const doublereal omega_b;

    //! Omega constant for the critical molar volume
    static const doublereal omega_vc;
};
}

#endif
