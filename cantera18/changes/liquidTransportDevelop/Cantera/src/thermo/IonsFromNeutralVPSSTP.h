/**
 *  @file IonsFromNeutralVPSSTP.h
 *   Header for intermediate ThermoPhase object for phases which
 *   consist of ions whose thermodynamics is calculated from neutral molecule thermodynamics.
 *  (see \ref thermoprops 
 * and class \link Cantera::IonsFromNeutralVPSSTP IonsFromNeutralVPSSTP\endlink).
 *
 * Header file for a derived class of %ThermoPhase that handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties that are further based upon activities
 * based on the molality scale.  These include most of the methods for
 * calculating liquid electrolyte thermodynamics.
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 *  $Id$
 */

#ifndef CT_IONSFROMNEUTRALVPSSTP_H
#define CT_IONSFROMNEUTRALVPSSTP_H

#include "GibbsExcessVPSSTP.h"

namespace Cantera {

  //! enums for molten salt ion solution types
  /*!
   *  Types identify how complicated the solution is. If there
   *  is just mixing on one of the sublattices but not the other,
   *  then the math is considerably simpler.
   */
  enum IonSolnType_enumType {
    cIonSolnType_PASSTHROUGH = 2000 ,
    cIonSolnType_SINGLEANION ,
    cIonSolnType_SINGLECATION ,
    cIonSolnType_MULTICATIONANION
  };

  /*!
   * The IonsFromNeutralVPSSTP is a derived class of ThermoPhase
   * that handles the specification of the chemical potentials for
   * ionic species, given a specification of the chemical potentials
   * for the same phase expressed in terms of combinations of the
   * ionic species that represent neutral molecules. It's expected
   * that the neutral molecules will be represented in terms of 
   * an excess gibbs free energy approximation that is a derivative
   * of the GbbsExcessVPSSTP object. All of the e Excess Gibbs free
   *  energy formulations in this area employ
   *  symmetrical formulations. 
   *
   *  This class is used for molten salts.
   *
   *  This object actually employs 4 different mole fraction types.
   *  
   *  1) There is a mole fraction associated the the cations and
   *     anions and neutrals from this ThermoPhase object. This
   *     is the normal mole fraction vector for this object.
   *     Note, however, it isn't the appropriate mole fraction
   *     vector to use even for obtaining the correct ideal 
   *     free energies of mixing. 
   *  2) There is a mole fraction vector associated with the
   *     neutral molecule ThermoPhase object. 
   *  3) There is a mole fraction vector associated with the
   *     cation lattice. 
   *  4) There is a mole fraction vector associated with the
   *     anion lattice
   *
   *  This object can translate between any of the four mole
   *  fraction representations.
   *
   *
   */
  class IonsFromNeutralVPSSTP : public GibbsExcessVPSSTP  {

  public:
        
    /// Constructors 
    /*!
     * Default constructor
     */
    IonsFromNeutralVPSSTP();

    //! Construct and initialize an IonsFromNeutralVPSSTP object
    //! directly from an asci input file
    /*!
     * Working constructors
     *
     *  The two constructors below are the normal way
     *  the phase initializes itself. They are shells that call
     *  the routine initThermo(), with a reference to the
     *  XML database to get the info for the phase.
     *
     * @param inputFile Name of the input file containing the phase XML data
     *                  to set up the object
     * @param id        ID of the phase in the input file. Defaults to the
     *                  empty string.
     * @param neutralPhase   The object takes a neutralPhase ThermoPhase
     *                       object as input. It can either take a pointer
     *                       to an existing object in the parameter list,
     *                       in which case it does not own the object, or
     *                       it can construct a neutral Phase as a slave
     *                       object, in which case, it does own the slave
     *                       object, for purposes of who gets to destroy
     *                       the object.
     *                       If this parameter is zero, then a slave
     *                       neutral phase object is created and used.
     */
    IonsFromNeutralVPSSTP(std::string inputFile, std::string id = "", 
			    ThermoPhase *neutralPhase = 0);
  
    //! Construct and initialize an IonsFromNeutralVPSSTP object
    //! directly from an XML database
    /*!
     *  @param phaseRoot XML phase node containing the description of the phase
     *  @param id     id attribute containing the name of the phase.
     *                (default is the empty string)
     * @param neutralPhase   The object takes a neutralPhase ThermoPhase
     *                       object as input. It can either take a pointer
     *                       to an existing object in the parameter list,
     *                       in which case it does not own the object, or
     *                       it can construct a neutral Phase as a slave
     *                       object, in which case, it does own the slave
     *                       object, for purposes of who gets to destroy
     *                       the object.
     *                       If this parameter is zero, then a slave
     *                       neutral phase object is created and used.
     */
    IonsFromNeutralVPSSTP(XML_Node& phaseRoot,  std::string id = "", 
			  ThermoPhase *neutralPhase = 0);


    //! Copy constructor
    /*!
     *  Note this stuff will not work until the underlying phase
     *  has a working copy constructor
     *
     * @param b class to be copied
     */
    IonsFromNeutralVPSSTP(const IonsFromNeutralVPSSTP &b);

    /// Assignment operator
    /*!
     *
     * @param b class to be copied.
     */
    IonsFromNeutralVPSSTP & operator=(const IonsFromNeutralVPSSTP &b);

    /// Destructor. 
    virtual ~IonsFromNeutralVPSSTP();

    //! Duplication routine for objects which inherit from ThermoPhase.
    /*!
     *  This virtual routine can be used to duplicate ThermoPhase objects
     *  inherited from ThermoPhase even if the application only has
     *  a pointer to ThermoPhase to work with.
     */
    virtual ThermoPhase *duplMyselfAsThermoPhase() const;
    
    /**
     *   
     * @name  Utilities  
     * @{
     */

   
    //! Equation of state type flag.
    /*!
     * The ThermoPhase base class returns
     * zero. Subclasses should define this to return a unique
     * non-zero value. Known constants defined for this purpose are
     * listed in mix_defs.h. The MolalityVPSSTP class also returns
     * zero, as it is a non-complete class.
     */
    virtual int eosType() const;

   

    /**
     * @} 
     * @name  Molar Thermodynamic Properties 
     * @{
     */
    
    //! Return the Molar enthalpy. Units: J/kmol.
    /*!
     * This is calculated from the partial molar enthalpies of the species
     */
    virtual doublereal enthalpy_mole() const;

    /**
     * Molar internal energy. J/kmol. 
     *  *
     * This is calculated from the soln enthalpy and then
     * subtracting pV.
     */
    virtual doublereal intEnergy_mole() const;

    /**
     * Molar entropy. Units: J/kmol/K.
     * 
     * 
     */
    virtual doublereal entropy_mole() const;

    /**
     * Molar Gibbs free Energy for an ideal gas.
     * Units =  J/kmol.
     */
    virtual doublereal gibbs_mole() const;

    /**
     * Molar heat capacity at constant pressure. Units: J/kmol/K.
     * For an ideal gas mixture, 
     *
     */
    virtual doublereal cp_mole() const;

    /**
     * Molar heat capacity at constant volume. Units: J/kmol/K.
     *
     */
    virtual doublereal cv_mole() const;


    /**
     * @}
     * @name Utilities 
     * @{
     */


 

    /**
     * @}
     * @name Mechanical Properties
     * @{
     */

    /**
     * @} 
     * @name Potential Energy
     * 
     * Species may have an additional potential energy due to the
     * presence of external gravitation or electric fields. These
     * methods allow specifying a potential energy for individual
     * species.
     * @{
     */

    /**
     * @}
     * @name Activities, Standard States, and Activity Concentrations
     *
     * The activity \f$a_k\f$ of a species in solution is
     * related to the chemical potential by \f[ \mu_k = \mu_k^0(T)
     * + \hat R T \log a_k. \f] The quantity \f$\mu_k^0(T,P)\f$ is
     * the chemical potential at unit activity, which depends only
     * on temperature and pressure.
     * @{
     */

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
  
    //! Return the standard concentration for the kth species
    /*!
     * The standard concentration \f$ C^0_k \f$ used to normalize
     * the activity (i.e., generalized) concentration. In many cases, this quantity
     * will be the same for all species in a phase - for example,
     * for an ideal gas \f$ C^0_k = P/\hat R T \f$. For this
     * reason, this method returns a single value, instead of an
     * array.  However, for phases in which the standard
     * concentration is species-specific (e.g. surface species of
     * different sizes), this method may be called with an
     * optional parameter indicating the species.
     *
     *  Here we define the standard concentration as being equal to 1.0.
     *  Therefore, the kinetics operators will be dealing in unitless
     *  activities for all kinetics expressions involving the molten
     *  salts. This assignment is subject to further assessment.
     *
     * @param k Optional parameter indicating the species. The default
     *          is to assume this refers to species 0.
     * @return 
     *   Returns the standard concentration. The units are by definition
     *   dependent on the ThermoPhase and kinetics manager representation.
     */
    virtual doublereal standardConcentration(int k=0) const;


    //! Natural logarithm of the standard concentration of the kth species.
    /*!
     * @param k    index of the species (defaults to zero)
     */
    virtual doublereal logStandardConc(int k=0) const; 

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
     * The base %ThermoPhase class assigns the default quantities
     * of (kmol/m3) for all species.
     * Inherited classes are responsible for overriding the default 
     * values if necessary.
     *
     * @param uA Output vector containing the units
     *  uA[0] = kmol units - default  = 1
     *  uA[1] = m    units - default  = -nDim(), the number of spatial
     *                                dimensions in the Phase class.
     *  uA[2] = kg   units - default  = 0;
     *  uA[3] = Pa(pressure) units - default = 0;
     *  uA[4] = Temperature units - default = 0;
     *  uA[5] = time units - default = 0
     * @param k species index. Defaults to 0.
     * @param sizeUA output int containing the size of the vector.
     *        Currently, this is equal to 6.
     */
    virtual void getUnitsStandardConc(double *uA, int k = 0,
				      int sizeUA = 6) const;

   
    //! Get the array of non-dimensional molar-based activity coefficients at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     * @param ac Output vector of activity coefficients. Length: m_kk.
     */
    virtual void getActivityCoefficients(doublereal* ac) const;


 
    //@}
    /// @name  Partial Molar Properties of the Solution 
    //@{

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

    
    //! Returns an array of partial molar enthalpies for the species
    //! in the mixture.
    /*!
     * Units (J/kmol)
     *
     * For this phase, the partial molar enthalpies are equal to the
     * standard state enthalpies modified by the derivative of the
     * molality-based activity coefficent wrt temperature
     *
     *  \f[
     * \bar h_k(T,P) = h^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
     * \f]
     *
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;

    //! Returns an array of partial molar entropies for the species
    //! in the mixture.
    /*!
     * Units (J/kmol)
     *
     * For this phase, the partial molar enthalpies are equal to the
     * standard state enthalpies modified by the derivative of the
     * activity coefficent wrt temperature
     *
     *  \f[
     *   \bar s_k(T,P) = s^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
     *                              - R \ln( \gamma_k X_k)
     *                              - R T \frac{d \ln(\gamma_k) }{dT}
     *  \f]
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;

    //! Get the array of change in the log activity coefficients w.r.t. change in state (change temp, change mole fractions)
    /*!
     * This function is a virtual class, but it first appears in GibbsExcessVPSSTP
     * class and derived classes from GibbsExcessVPSSTP.
     *
     * This function is a virtual method.  For ideal mixtures 
     * (unity activity coefficients), this can gradX/X.  
     *
     * @param dT    Input of temperature change
     * @param dX    Input vector of changes in mole fraction. length = m_kk
     * @param dlnActCoeff    Output vector of derivatives of the 
     *                         log Activity Coefficients. length = m_kk
     */
     virtual void getdlnActCoeff(const doublereal dT, const doublereal * const dX, doublereal *dlnActCoeff) const;

    //! Get the array of log concentration-like derivatives of the 
    //! log activity coefficients
    /*!
     * This function is a virtual method.  For ideal mixtures 
     * (unity activity coefficients), this can return zero.  
     * Implementations should take the derivative of the 
     * logarithm of the activity coefficient with respect to the 
     * logarithm of the concentration-like variable (i.e. mole fraction) 
     * that represents the standard state.  
     * This quantity is to be used in conjunction with derivatives of 
     * that concentration-like variable when the derivative of the chemical 
     * potential is taken.  
     *
     *  units = dimensionless
     *
     * @param dlnActCoeffdlnX    Output vector of log(mole fraction)  
     *                 derivatives of the log Activity Coefficients.
     *                 length = m_kk
     */
    virtual void getdlnActCoeffdlnX(doublereal *dlnActCoeffdlnX) const;

    //! Get the array of log concentration-like derivatives of the 
    //! log activity coefficients
    /*!
     * This function is a virtual method.  For ideal mixtures 
     * (unity activity coefficients), this can return zero.  
     * Implementations should take the derivative of the 
     * logarithm of the activity coefficient with respect to the 
     * logarithm of the concentration-like variable (i.e. number of moles)
     * that represents the standard state.  
     * This quantity is to be used in conjunction with derivatives of 
     * that concentration-like variable when the derivative of the chemical 
     * potential is taken.  
     *
     *  units = dimensionless
     *
     * @param dlnActCoeffdlnN    Output vector of log(mole fraction)  
     *                 derivatives of the log Activity Coefficients.
     *                 length = m_kk
     */
    virtual void getdlnActCoeffdlnN(doublereal *dlnActCoeffdlnN) const;

    //! Get the Salt Dissociation Coefficients
    //! Returns the vector of dissociation coefficients and vector of charges
    virtual void getDissociationCoeffs(vector_fp& coeffs, vector_fp& charges, std::vector<int>& neutMolIndex);

    virtual void getNeutralMolecMoleFractions(vector_fp& fracs){fracs=NeutralMolecMoleFractions_;}

    //! Calculate neutral molecule mole fractions
    /*!
     *  This routine calculates the neutral molecule mole
     *  fraction given the vector of ion mole fractions,
     *  i.e., the mole fractions from this ThermoPhase.
     *  Note, this routine basically assumes that there
     *  is charge neutrality. If there isn't, then it wouldn't
     *  make much sense. 
     *
     *  for the case of  cIonSolnType_SINGLEANION, some slough
     *  in the charge neutrality is allowed. The cation number
     *  is followed, while the difference in charge neutrality
     *  is dumped into the anion mole number to fix the imbalance.
     */
    virtual void getNeutralMoleculeMoleGrads(const doublereal * const x, doublereal *y) const;

    virtual void getCationList(std::vector<int>& cation){cation=cationList_;}
    virtual void getAnionList(std::vector<int>& anion){anion=anionList_;}
    virtual void getSpeciesNames(std::vector<std::string>& names){names=m_speciesNames;}


    //@}
    /// @name  Properties of the Standard State of the Species in the Solution
    //@{

     

    //@}
    /// @name Thermodynamic Values for the Species Reference States
    //@{


    ///////////////////////////////////////////////////////
    //
    //  The methods below are not virtual, and should not
    //  be overloaded.
    //
    //////////////////////////////////////////////////////

    /**
     * @name Specific Properties
     * @{
     */


    /**
     * @name Setting the State
     *
     * These methods set all or part of the thermodynamic
     * state.
     * @{
     */

    virtual void setTemperature(const doublereal t);
    virtual void setPressure(doublereal p);

    //! Set the temperature (K) and pressure (Pa)
    /*!
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param t    Temperature (K)
     * @param p    Pressure (Pa)
     */
    virtual void setState_TP(doublereal t, doublereal p);


    //! Calculate ion mole fractions from neutral molecule 
    //! mole fractions.
    /*!
     *  @param mf Dump the mole fractions into this vector.
     */
    virtual void calcIonMoleFractions(doublereal *const mf) const;

    //! Calculate neutral molecule mole fractions
    /*!
     *  This routine calculates the neutral molecule mole
     *  fraction given the vector of ion mole fractions,
     *  i.e., the mole fractions from this ThermoPhase.
     *  Note, this routine basically assumes that there
     *  is charge neutrality. If there isn't, then it wouldn't
     *  make much sense. 
     *
     *  for the case of  cIonSolnType_SINGLEANION, some slough
     *  in the charge neutrality is allowed. The cation number
     *  is followed, while the difference in charge neutrality
     *  is dumped into the anion mole number to fix the imbalance.
     */
    virtual void calcNeutralMoleculeMoleFractions() const;

   /**
     * Set the mass fractions to the specified values, and then
     * normalize them so that they sum to 1.0.
     * @param y Array of unnormalized mass fraction values (input).
     * Must have a length greater than or equal to the number of
     * species.
     *
     * @param y  Input vector of mass fractions.
     *           Length is m_kk.
     */
    virtual void setMassFractions(const doublereal* const y);

    /**
     * Set the mass fractions to the specified values without
     * normalizing. This is useful when the normalization
     * condition is being handled by some other means, for example
     * by a constraint equation as part of a larger set of
     * equations.
     *
     * @param y  Input vector of mass fractions.
     *           Length is m_kk.
     */
    virtual void setMassFractions_NoNorm(const doublereal* const y);

    /**
     * Set the mole fractions to the specified values, and then
     * normalize them so that they sum to 1.0.
     * @param x Array of unnormalized mole fraction values (input).
     * Must have a length greater than or equal to the number of
     * species.
     *
     * @param x  Input vector of mole fractions.
     *           Length is m_kk.
     */
    virtual void setMoleFractions(const doublereal* const x);

    /**
     * Set the mole fractions to the specified values without
     * normalizing. This is useful when the normalization
     * condition is being handled by some other means, for example
     * by a constraint equation as part of a larger set of
     * equations.
     *
     * @param x  Input vector of mole fractions.
     *           Length is m_kk.
     */
    virtual void setMoleFractions_NoNorm(const doublereal* const x);

    /**
     * Set the concentrations to the specified values within the
     * phase.
     *
     * @param c The input vector to this routine is in dimensional
     *        units. For volumetric phases c[k] is the
     *        concentration of the kth species in kmol/m3.
     *        For surface phases, c[k] is the concentration
     *        in kmol/m2. The length of the vector is the number
     *        of species in the phase.
     */
    virtual void setConcentrations(const doublereal* const c);

    //@}

    /**
     * @name Chemical Equilibrium
     * Routines that implement the Chemical equilibrium capability
     * for a single phase, based on the element-potential method.
     * @{
     */

   

    //@}



    /// The following methods are used in the process of constructing
    /// the phase and setting its parameters from a specification in an 
    /// input file. They are not normally used in application programs.
    /// To see how they are used, see files importCTML.cpp and 
    /// ThermoFactory.cpp.

    //! Initialization of an IonsFromNeutralVPSSTP phase using an xml file
    /*!
     * This routine is a precursor to initThermo(XML_Node*)
     * routine, which does most of the work.
     *
     * @param inputFile XML file containing the description of the
     *        phase
     *
     * @param id  Optional parameter identifying the name of the
     *            phase. If none is given, the first XML
     *            phase element will be used.
     */
    void constructPhaseFile(std::string inputFile, std::string id);

    //!   Import and initialize an IonsFromNeutralVPSSTP phase
    //!   specification in an XML tree into the current object.
    /*!
     *   Here we read an XML description of the phase.
     *   We import descriptions of the elements that make up the
     *   species in a phase.
     *   We import information about the species, including their
     *   reference state thermodynamic polynomials. We then freeze
     *   the state of the species.
     *
     *   Then, we read the species molar volumes from the xml
     *   tree to finish the initialization.
     *
     * @param phaseNode This object must be the phase node of a
     *             complete XML tree
     *             description of the phase, including all of the
     *             species data. In other words while "phase" must
     *             point to an XML phase object, it must have
     *             sibling nodes "speciesData" that describe
     *             the species in the phase.
     *
     * @param id   ID of the phase. If nonnull, a check is done
     *             to see if phaseNode is pointing to the phase
     *             with the correct id.
     */
    void constructPhaseXML(XML_Node& phaseNode, std::string id);


    /*!
     * @internal Initialize. This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species.  The base class implementation does nothing,
     * and subclasses that do not require initialization do not
     * need to overload this method.  When importing a CTML phase
     * description, this method is called just prior to returning
     * from function importPhase.
     *
     * @see importCTML.cpp
     */
    virtual void initThermo();


    /**
     *   Import and initialize a ThermoPhase object
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
    void initThermoXML(XML_Node& phaseNode, std::string id);

 
    //! returns a summary of the state of the phase as a string
    /*!
     * @param show_thermo If true, extra information is printed out
     *                    about the thermodynamic state of the system.
     */
    virtual std::string report(bool show_thermo = true) const;

    //! returns a summary of the state of the phase to specified 
    //! comma separated files
    /*!
     * @param textFile    ofstream file to print textual variable names that
     *                    will correspod to data in comma separated file (csv).
     * @param csvFile     ofstream file to print comma separated data for
     *                    the phase
     * @param show_thermo If true, extra information is printed out
     *                    about the thermodynamic state of the system.
     */
    virtual void reportCSV(std::ofstream& textFile, std::ofstream& csvFile, bool show_thermo = true) const;

  private:
  

    //! Initialize lengths of local variables after all species have
    //! been identified.
    void initLengths();
            
    //! Update the activity coefficients
    /*!
     * This function will be called to update the internally storred
     * natural logarithm of the activity coefficients
     */
    void s_update_lnActCoeff() const;

    //! Update the temperature derivative of the ln activity coefficients
    /*!
     * This function will be called to update the internally storred
     * temperature derivative of the natural logarithm of the activity coefficients
     */
    void s_update_dlnActCoeffdT() const;

    //! Update the change in the ln activity coefficients
    /*!
     * This function will be called to update the internally storred
     * change of the natural logarithm of the activity coefficients
     * w.r.t a change in state (temp, mole fraction, etc)
     */
    void s_update_dlnActCoeff() const;

    //! Update the derivative of the log of the activity coefficients
    //!  wrt log(mole fraction)
    /*!
     * This function will be called to update the internally storred
     * derivative of the natural logarithm of the activity coefficients
     * wrt logarithm of the mole fractions.
     */
    void s_update_dlnActCoeff_dlnX() const;

    //! Update the derivative of the log of the activity coefficients
    //!  wrt log(number of moles)
    /*!
     * This function will be called to update the internally storred
     * derivative of the natural logarithm of the activity coefficients
     * wrt logarithm of the number of moles of given species.
     */
    void s_update_dlnActCoeff_dlnN() const;


  private:
    //! Error function
    /*!
     *  Print an error string and exit
     *
     * @param msg  Message to be printed
     */
    doublereal err(std::string msg) const;

  protected:

    //! Ion solution type
    /*!
     *  There is either mixing on the anion, cation, or both lattices.
     *  There is also a passthrough option
     * 
     *  Defaults to cIonSolnType_SINGLEANION, so that LiKCl can be hardwired
     */
    IonSolnType_enumType ionSolnType_;

    //! Number of neutral molecule species
    /*!
     *  This is equal to the number of species in the 
     *  neutralMoleculePhase_ ThermoPhase.
     */
    int numNeutralMoleculeSpecies_;

    //! Index of special species
    int indexSpecialSpecies_;

    //! Index of special species
    int indexSecondSpecialSpecies_;
    
    //! Formula Matrix for composition of neutral molecules
    //! in terms of the molecules in this ThermoPhase
    /*!
     *       fm_neutralMolec_ions[ i + jNeut * m_kk                ]
     *
     *             This is the number of ions of type i in the neutral
     *             molecule jNeut.
     */
    std::vector<double> fm_neutralMolec_ions_;

    //! Mapping between ion species and neutral molecule for quick invert.
    /*!
     *
     * fm_invert_ionForNeutral returns vector of int. Each element represents 
     * an ionic species and stores the value of the corresponding neutral 
     * molecule
     *
     *  For the case of fm_invert_simple_ = true, we assume that there
     *  is a quick way to invert the formula matrix so that we can
     *  quickly calculate the neutral molecule mole fraction
     *  given the ion mole fraction vector. 
     *
     *  We assume that for a selected set of ion species, that that
     *  ion is only in the neutral molecule, jNeut.
     *
     *     therefore,    
     *
     *      NeutralMolecMoleFractions_[jNeut] += moleFractions_[i_ion] / fmij;
     *
     *   where fmij is the number of ions in neutral molecule jNeut.
     *
     *  Thus, we formulate the neutral molecule mole fraction NeutralMolecMoleFractions_[]
     *  vector from this association. We further assume that there are
     *  no other associations.  If  fm_invert_simple_ is not true,
     *  then we need to do a formal inversion which takes a great
     *  deal of time and is not currently implemented.
     */
    std::vector<int> fm_invert_ionForNeutral;

    //! Mole fractions using the Neutral Molecule Mole fraction basis
    mutable std::vector<doublereal> NeutralMolecMoleFractions_;


    //! List of the species in this ThermoPhase which are cation species
    std::vector<int> cationList_;

    //! Number of cation species
    int numCationSpecies_;

    //! List of the species in this ThermoPhase which are anion species
    std::vector<int> anionList_;

    //! Number of anion species
    int numAnionSpecies_;

    //! List of the species in this ThermoPhase which are passed
    //! through to the neutralMoleculePhase ThermoPhase.
    /*!
     *  These have neutral charges.
     */ 
    std::vector<int> passThroughList_;

    //! Number of the species in this ThermoPhase which are passed
    //! through to the neutralMoleculePhase ThermoPhase
    int numPassThroughSpecies_;

  public:
    //! This is a pointer to the neutral Molecule Phase
    /*!
     *  If the variable, IOwnNThermoPhase_ is true, then we own
     *  the pointer. If not, then this is considered a shallow pointer.
     */
    ThermoPhase *neutralMoleculePhase_;

  private:

    //! If true then we own the underlying neutral Molecule Phase
    /*!
     *  If this is false, then the neutral molecule phase is considered
     *  as a shallow pointer.
     */
    bool IOwnNThermoPhase_;

    //! ThermoPhase for the cation lattice
    /*!
     *  Currently this is unimplemented and may be deleted
     */
    // ThermoPhase *cationPhase_;

    //! ThermoPhase for the anion lattice
    /*!
     *  Currently this is unimplemented and may be deleted
     */
    //ThermoPhase *anionPhase_;
   
    //! Temporary mole fraction vector
    mutable std::vector<doublereal> moleFractionsTmp_;

    mutable std::vector<doublereal> muNeutralMolecule_;
    mutable std::vector<doublereal> gammaNeutralMolecule_;
    mutable std::vector<doublereal> dlnActCoeff_NeutralMolecule_;
    mutable std::vector<doublereal> dlnActCoeffdT_NeutralMolecule_;
    mutable std::vector<doublereal> dlnActCoeffdlnX_NeutralMolecule_;
    mutable std::vector<doublereal> dlnActCoeffdlnN_NeutralMolecule_;

  };





}
        
#endif
