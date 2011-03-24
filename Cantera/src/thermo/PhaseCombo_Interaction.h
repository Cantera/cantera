/**
 *  @file 
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 *  $Id: MargulesVPSSTP.h 641 2010-11-12 21:37:41Z hkmoffa $
 */

#ifndef CT_LIFES_X_ONEPHASE_VPSSTP_H
#define CT_LIFES_X_ONEPHASE_VPSSTP_H

#include "PseudoBinaryVPSSTP.h"
#include "GibbsExcessVPSSTP.h"

namespace Cantera {

  /**
   * @ingroup thermoprops
   */

  
  //!  MargulesVPSSTP is a derived class of GibbsExcessVPSSTP that employs
  //!  the Margules approximation for the excess gibbs free energy

  class PhaseCombo_Interaction : public GibbsExcessVPSSTP {

  public:
        
    //! Constructor 
    /*!
     * This doesn't do much more than initialize constants with
     * default values for water at 25C. Water molecular weight 
     * comes from the default elements.xml file. It actually
     * differs slightly from the IAPWS95 value of 18.015268. However,
     * density conservation and therefore element conservation
     * is the more important principle to follow.
     */
    PhaseCombo_Interaction();

    //! Construct and initialize a PhaseCombo_Interaction ThermoPhase object
    //! directly from an xml input file
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
     */
    PhaseCombo_Interaction(std::string inputFile, std::string id = "");

    //! Construct and initialize a PhaseCombo_Interaction ThermoPhase object
    //! directly from an XML database
    /*!
     *  @param phaseRef XML phase node containing the description of the phase
     *  @param id     id attribute containing the name of the phase.
     *                (default is the empty string)
     */
    PhaseCombo_Interaction(XML_Node& phaseRef, std::string id = "");


    //! Special constructor for a hard-coded problem
    /*!
     * 
     *  @param testProb Hard-coded value. Only the value of 1 is
     *                  used. It's for 
     *                  a LiKCl system
     *                  -> test to predict the eutectic and liquidus correctly.
     */
    PhaseCombo_Interaction(int testProb);

    //! Copy constructor
    /*!
     *  Note this stuff will not work until the underlying phase
     *  has a working copy constructor
     *
     * @param b class to be copied
     */
    PhaseCombo_Interaction(const PhaseCombo_Interaction& b);

    //! Assignment operator
    /*!
     *
     * @param b class to be copied.
     */
    PhaseCombo_Interaction& operator=(const PhaseCombo_Interaction &b);

    //! Destructor
    virtual ~PhaseCombo_Interaction();

    //! Duplication routine for objects which inherit from  ThermoPhase.
    /*!
     *  This virtual routine can be used to duplicate thermophase objects
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

    //! Initialization of a phase using an xml file
    /*!
     * This routine is a precursor to 
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

    //!   Import and initialize a phase 
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

    /**
     * @} 
     * @name  Molar Thermodynamic Properties 
     * @{
     */


    /**
     * @}
     * @name Utilities for Solvent ID and Molality
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
 

    /**
     * The standard concentration \f$ C^0_k \f$ used to normalize
     * the generalized concentration. In many cases, this quantity
     * will be the same for all species in a phase - for example,
     * for an ideal gas \f$ C^0_k = P/\hat R T \f$. For this
     * reason, this method returns a single value, instead of an
     * array.  However, for phases in which the standard
     * concentration is species-specific (e.g. surface species of
     * different sizes), this method may be called with an
     * optional parameter indicating the species.
     *
     * @param k species index. Defaults to zero.
     */
    virtual doublereal standardConcentration(int k=0) const;

    /**
     * Returns the natural logarithm of the standard 
     * concentration of the kth species
     *
     * @param k  species index
     */
    virtual doublereal logStandardConc(int k=0) const;

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

      /// Molar enthalpy. Units: J/kmol. 
    virtual doublereal enthalpy_mole() const;

      /// Molar entropy. Units: J/kmol. 
    virtual doublereal entropy_mole() const;

    /// Molar heat capacity at constant pressure. Units: J/kmol/K. 
    virtual doublereal cp_mole() const;

    /// Molar heat capacity at constant volume. Units: J/kmol/K. 
    virtual doublereal cv_mole() const;

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
     *   \bar h_k(T,P) = h^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
     *  \f]
     *
     * @param hbar  Vector of returned partial molar enthalpies
     *              (length m_kk, units = J/kmol)
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
     *
     * @param sbar  Vector of returned partial molar entropies
     *              (length m_kk, units = J/kmol/K)
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;

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
     *   ???????????????
     *   \bar s_k(T,P) = s^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
     *                              - R \ln( \gamma_k X_k)
     *                              - R T \frac{d \ln(\gamma_k) }{dT}
     *   ???????????????
     *  \f]
     *
     * @param cpbar  Vector of returned partial molar heat capacities
     *              (length m_kk, units = J/kmol/K)
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const;

    
    //! Return an array of partial molar volumes for the
    //! species in the mixture. Units: m^3/kmol.
    /*!
     *  Frequently, for this class of thermodynamics representations,
     *  the excess Volume due to mixing is zero. Here, we set it as
     *  a default. It may be overriden in derived classes.
     *
     *  @param vbar   Output vector of speciar partial molar volumes.
     *                Length = m_kk. units are m^3/kmol.
     */
    virtual void getPartialMolarVolumes(doublereal* vbar) const;

    //! Get the species electrochemical potentials.
    /*!
     * These are partial molar quantities.
     * This method adds a term \f$ Fz_k \phi_k \f$ to the 
     * to each chemical potential.
     *
     * Units: J/kmol
     *
     * @param mu     output vector containing the species electrochemical potentials.
     *               Length: m_kk., units = J/kmol
     */
    void getElectrochemPotentials(doublereal* mu) const;

    //! Get the array of temperature second derivatives of the log activity coefficients
    /*!
     * This function is a virtual class, but it first appears in GibbsExcessVPSSTP
     * class and derived classes from GibbsExcessVPSSTP.
     *
     *  units = 1/Kelvin
     *
     * @param d2lnActCoeffdT2  Output vector of temperature 2nd derivatives of the 
     *                         log Activity Coefficients. length = m_kk
     *
     */
    virtual void getd2lnActCoeffdT2(doublereal *d2lnActCoeffdT2) const;

    //! Get the array of temperature derivatives of the log activity coefficients
    /*!
     * This function is a virtual class, but it first appears in GibbsExcessVPSSTP
     * class and derived classes from GibbsExcessVPSSTP.
     *
     *  units = 1/Kelvin
     *
     * @param dlnActCoeffdT    Output vector of temperature derivatives of the 
     *                         log Activity Coefficients. length = m_kk
     *
     */
    virtual void getdlnActCoeffdT(doublereal *dlnActCoeffdT) const;


 
    //@}
    /// @name  Properties of the Standard State of the Species in the Solution
    //@{

     

    //@}
    /// @name Thermodynamic Values for the Species Reference States
    //@{


   


    /// The following methods are used in the process of constructing
    /// the phase and setting its parameters from a specification in an 
    /// input file. They are not normally used in application programs.
    /// To see how they are used, see files importCTML.cpp and 
    /// ThermoFactory.cpp.


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

    /**
     * @} 
     * @name  Derivatives of Thermodynamic Variables needed for Applications
     * @{
     */

    //! Get the change in activity coefficients w.r.t. change in state (temp, mole fraction, etc.) along
    //! a line in parameter space or along a line in physical space
    /*!
     *
     * @param dTds           Input of temperature change along the path
     * @param dXds           Input vector of changes in mole fraction along the path. length = m_kk
     *                       Along the path length it must be the case that the mole fractions sum to one.
     * @param dlnActCoeffds  Output vector of the directional derivatives of the 
     *                       log Activity Coefficients along the path. length = m_kk
     *  units are 1/units(s). if s is a physical coordinate then the units are 1/m.
     */
    virtual void getdlnActCoeffds(const doublereal dTds, const doublereal * const dXds, doublereal *dlnActCoeffds) const;
 
    //! Get the array of log concentration-like derivatives of the 
    //! log activity coefficients - diagonal component
    /*!
     * This function is a virtual method.  For ideal mixtures 
     * (unity activity coefficients), this can return zero.  
     * Implementations should take the derivative of the 
     * logarithm of the activity coefficient with respect to the 
     * logarithm of the mole fraction.
     *
     *  units = dimensionless
     *
     * @param dlnActCoeffdlnX_diag    Output vector of the diagonal component of the log(mole fraction)  
     *                 derivatives of the log Activity Coefficients.
     *                 length = m_kk
     */
    virtual void getdlnActCoeffdlnX_diag(doublereal *dlnActCoeffdlnX_diag) const;

    //! Get the array of  derivatives of the log activity coefficients wrt mole numbers - diagonal only
    /*!
     * This function is a virtual method.  For ideal mixtures 
     * (unity activity coefficients), this can return zero.  
     * Implementations should take the derivative of the 
     * logarithm of the activity coefficient with respect to the 
     * logarithm of the concentration-like variable (i.e. mole fraction,
     * molality, etc.) that represents the standard state.  
     *
     *  units = dimensionless
     *
     * @param dlnActCoeffdlnN_diag    Output vector of the diagonal entries for the log(mole fraction)  
     *                 derivatives of the log Activity Coefficients.
     *                 length = m_kk
     */
    virtual void getdlnActCoeffdlnN_diag(doublereal *dlnActCoeffdlnN_diag) const;


    //! Get the array of derivatives of the log activity coefficients with respect to the ln species mole numbers
    /*!
     * Implementations should take the derivative of the logarithm of the activity coefficient with respect to a
     * log of a species mole number (with all other species mole numbers held constant)
     * 
     *  units = 1 / kmol
     *
     *  dlnActCoeffdlnN[ ld * k  + m]  will contain the derivative of log act_coeff for the <I>m</I><SUP>th</SUP> 
     *                                 species with respect to the number of moles of the <I>k</I><SUP>th</SUP> species.
     *
     * \f[
     *        \frac{d \ln(\gamma_m) }{d \ln( n_k ) }\Bigg|_{n_i}
     * \f]
     *
     * @param ld               Number of rows in the matrix
     * @param dlnActCoeffdlnN    Output vector of derivatives of the 
     *                         log Activity Coefficients. length = m_kk * m_kk        
     */
    virtual void getdlnActCoeffdlnN(const int ld, doublereal * const dlnActCoeffdlnN) const;

   //@}

  private:
  
    //! Process an XML node called "binaryNeutralSpeciesParameters"
    /*!
     * This node contains all of the parameters necessary to describe
     * the Margules model for a particular binary interaction.
     * This function reads the XML file and writes the coefficients
     * it finds to an internal data structures.
     *
     * @param xmlBinarySpecies  Reference to the XML_Node named "binaryNeutralSpeciesParameters"
     *                          containing the binary interaction
     */
    void readXMLBinarySpecies(XML_Node &xmlBinarySpecies);

    //! Resize internal arrays within the object that depend upon the number
    //! of binary Margules interaction terms
    /*!
     *  @param num Number of binary Margules interaction terms
     */
    void resizeNumInteractions(const int num);


    //! Initialize lengths of local variables after all species have
    //! been identified.
    void initLengths();

    //! Update the activity coefficients
    /*!
     * This function will be called to update the internally storred
     * natural logarithm of the activity coefficients
     */
    void s_update_lnActCoeff() const;

    //! Update the derivative of the log of the activity coefficients wrt T
    /*!
     * This function will be called to update the internally storred
     * derivative of the natural logarithm of the activity coefficients
     * wrt temperature.
     */
    void s_update_dlnActCoeff_dT() const;

    //! Update the derivative of the log of the activity coefficients
    //!  wrt log(mole fraction)
    /*!
     * This function will be called to update the internally storred
     * derivative of the natural logarithm of the activity coefficients
     * wrt logarithm of the mole fractions.
     */
    void s_update_dlnActCoeff_dlnX_diag() const;

    //! Update the derivative of the log of the activity coefficients
    //!  wrt log(moles) - diagonal only
    /*!
     * This function will be called to update the internally storred diagonal entries for the
     * derivative of the natural logarithm of the activity coefficients
     * wrt logarithm of the moles.
     */
    void s_update_dlnActCoeff_dlnN_diag() const;

    //! Update the derivative of the log of the activity coefficients  wrt log(moles_m)
    /*!
     * This function will be called to update the internally storred
     * derivative of the natural logarithm of the activity coefficients
     * wrt logarithm of the mole number of species
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


    //! number of binary interaction expressions
    int numBinaryInteractions_;

    //! Enthalpy term for the binary mole fraction interaction of the
    //! excess gibbs free energy expression
    mutable vector_fp m_HE_b_ij;

    //! Enthalpy term for the ternary mole fraction interaction of the
    //! excess gibbs free energy expression
    mutable vector_fp m_HE_c_ij;

    //! Enthalpy term for the quaternary mole fraction interaction of the
    //! excess gibbs free energy expression
    mutable vector_fp m_HE_d_ij;

    //! Entropy term for the binary mole fraction interaction of the
    //! excess gibbs free energy expression
    mutable vector_fp m_SE_b_ij;

    //! Entropy term for the ternary mole fraction interaction of the
    //! excess gibbs free energy expression
    mutable vector_fp m_SE_c_ij;

    //! Entropy term for the quaternary mole fraction interaction of the
    //! excess gibbs free energy expression
    mutable vector_fp m_SE_d_ij;

    //! Enthalpy term for the binary mole fraction interaction of the
    //! excess gibbs free energy expression
    mutable vector_fp m_VHE_b_ij;

    //! Enthalpy term for the ternary mole fraction interaction of the
    //! excess gibbs free energy expression
    mutable vector_fp m_VHE_c_ij;

    //! Enthalpy term for the quaternary mole fraction interaction of the
    //! excess gibbs free energy expression
    mutable vector_fp m_VHE_d_ij;

    //! Entropy term for the binary mole fraction interaction of the
    //! excess gibbs free energy expression
    mutable vector_fp m_VSE_b_ij;

    //! Entropy term for the ternary mole fraction interaction of the
    //! excess gibbs free energy expression
    mutable vector_fp m_VSE_c_ij;

    //! Entropy term for the quaternary mole fraction interaction of the
    //! excess gibbs free energy expression
    mutable vector_fp m_VSE_d_ij;


    
    //! vector of species indices representing species A in the interaction
    /*!
     *  Each Margules excess Gibbs free energy term involves two species, A and B.
     *  This vector identifies species A.
     */
    vector_int m_pSpecies_A_ij;

    //! vector of species indices representing species B in the interaction
    /*!
     *  Each Margules excess Gibbs free energy term involves two species, A and B.
     *  This vector identifies species B.
     */
    vector_int m_pSpecies_B_ij;

    //! form of the Margules interaction expression
    /*!
     *  Currently there is only one form.
     */
    int formMargules_;

    //! form of the temperatuer dependence of the Margules interaction expression
    /*!
     *  Currently there is only one form -> constant wrt temperature.
     */
    int formTempModel_;

  
  };



}
        
#endif





