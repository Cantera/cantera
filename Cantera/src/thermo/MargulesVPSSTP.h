/**
 *  @file  MargulesVPSSTP.h
 *   Header for intermediate ThermoPhase object for phases which
 *   employ gibbs excess free energy based formulations
 *  (see \ref thermoprops 
 * and class \link Cantera::MargulesVPSSTP MargulesVPSSTP\endlink).
 *
 * Header file for a derived class of ThermoPhase that handles
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
 *  $Id: MargulesVPSSTP.h,v 1.1 2009/03/03 21:08:31 hkmoffa Exp $
 */

#ifndef CT_MARGULESVPSSTP_H
#define CT_MARGULESVPSSTP_H

#include "PseudoBinaryVPSSTP.h"

namespace Cantera {

  /**
   * @ingroup thermoprops
   */

  
  //!  MargulesVPSSTP is a derived class of GibbsExcessVPSSTP that employs
  //!  the Margules approximation for the excess gibbs free energy
  /*!
   *  
   * %MargulesVPSSTP derives from class GibbsExcessVPSSTP which is derived
   * from VPStandardStateTP,
   * and overloads the virtual methods defined there with ones that
   * use expressions appropriate for the Margules Excess gibbs free energy
   * approximation.
   *
   * The independent unknowns are pressure, temperature, and mass fraction.
   *
   *  Several concepts are introduced. The first concept is there are temporary
   *  variables for holding the species standard state values
   *  of Cp, H, S, G, and V at the
   *  last temperature and pressure called. These functions are not recalculated
   *  if a new call is made using the previous temperature and pressure. Currently,
   *  these variables and the calculation method are handled by the VPSSMgr class,
   *  for which VPStandardStateTP owns a pointer to.
   *
   *  To support the above functionality, pressure and temperature variables,
   *  m_plast_ss and m_tlast_ss, are kept which store the last pressure and temperature
   *  used in the evaluation of standard state properties.
   *
   *  This class is usually used for nearly incompressible phases. For those phases, it
   *  makes sense to change the equation of state independent variable from
   *  density to pressure. The variable m_Pcurrent contains the current value of the
   *  pressure within the phase.
   *
   *
   * <HR>
   * <H2> Specification of Species Standard %State Properties </H2>
   * <HR>
   *
   *  All species are defined to have standard states that depend upon both
   *  the temperature and the pressure. The Margules approximation assumes
   *  symmetric standard states, where all of the standard state assume
   *  that the species are in pure component states at the temperatue
   *  and pressure of the solution.  I don't think it prevents, however,
   *  some species from being dilute in the solution.
   *
   *
   * <HR>
   * <H2> Specification of Solution Thermodynamic Properties </H2>
   * <HR>
   *
   * The excess Gibbs free energy
   *
   *      \f[
   *         G^E = \sum_i \left(  H_{Ei} - T S_{Ei} \right)
   *      \f]
   *      \f[
   *         H^E_i = X_{Ai} X_{Bi} \left( h_{o,i} +  h_{1,i} X_{Bi} \right)
   *      \f]
   *      \f[
   *         S^E_i = X_{Ai} X_{Bi} \left( s_{o,i} +  s_{1,i} X_{Bi} \right)
   *      \f]
   *      
   *
   * The activity of a species defined in the phase is given by an excess 
   * Gibbs free energy formulation.
   *
   *       \f[
   *            a_k = \gamma_k  X_k
   *       \f]
   *
   * where \f$ X_k \f$ is the mole fraction of species <I>k</I>.
   * The chemical potential for species <I>k</I> is equal to
   *
   *       \f[
   *            \mu_k(T,P) = \mu^o_k(T, P) + R T \log(\gamma_k X_k)
   *       \f]
   *
   * In terms of the reference state, the above can be rewritten
   *
   *
   *       \f[
   *            \mu_k(T,P) = \mu^{ref}_k(T, P) + R T \log(\frac{P X_k}{P_{ref}})
   *       \f]
   *
   * The partial molar entropy for species <I>k</I> is given by the following relation,
   *
   *       \f[
   *            \tilde{s}_k(T,P) = s^o_k(T,P) - R \log(X_k) = s^{ref}_k(T) - R \log(\frac{P X_k}{P_{ref}})
   *       \f]
   *
   * The partial molar enthalpy for species <I>k</I> is
   *
   *       \f[
   *            \tilde{h}_k(T,P) = h^o_k(T,P) = h^{ref}_k(T)
   *       \f]
   *
   * The partial molar Internal Energy for species <I>k</I> is
   *
   *       \f[
   *            \tilde{u}_k(T,P) = u^o_k(T,P) = u^{ref}_k(T)
   *       \f]
   *
   * The partial molar Heat Capacity for species <I>k</I> is
   *
   *       \f[
   *            \tilde{Cp}_k(T,P) = Cp^o_k(T,P) = Cp^{ref}_k(T)
   *       \f]
   *
   * <HR>
   * <H2> %Application within %Kinetics Managers </H2>
   * <HR>
   *
   *   \f$ C^a_k\f$ are defined such that \f$ a_k = C^a_k /
   *   C^s_k, \f$ where \f$ C^s_k \f$ is a standard concentration
   *   defined below and \f$ a_k \f$ are activities used in the
   *   thermodynamic functions.  These activity (or generalized)
   *   concentrations are used
   *   by kinetics manager classes to compute the forward and
   *   reverse rates of elementary reactions.
   *   The activity concentration,\f$  C^a_k \f$,is given by the following expression.
   *
   *       \f[
   *            C^a_k = C^s_k  X_k  = \frac{P}{R T} X_k
   *       \f]
   *
   * The standard concentration for species <I>k</I> is independent of <I>k</I> and equal to
   *
   *        \f[
   *            C^s_k =  C^s = \frac{P}{R T}
   *        \f]
   *
   * For example, a bulk-phase binary gas reaction between species j and k, producing
   * a new gas species l would have the
   * following equation for its rate of progress variable, \f$ R^1 \f$, which has
   * units of kmol m-3 s-1.
   *
   *   \f[
   *    R^1 = k^1 C_j^a C_k^a =  k^1 (C^s a_j) (C^s a_k)
   *   \f]
   *  where
   *   \f[
   *      C_j^a = C^s a_j \mbox{\quad and \quad} C_k^a = C^s a_k
   *   \f]
   *
   *
   *  \f$ C_j^a \f$ is the activity concentration of species j, and
   *  \f$ C_k^a \f$ is the activity concentration of species k. \f$ C^s \f$
   *  is the standard concentration. \f$ a_j \f$ is
   *  the activity of species j which is equal to the mole fraction of j.
   *
   *  The reverse rate constant can then be obtained from the law of microscopic reversibility
   *  and the equilibrium expression for the system.
   *
   *   \f[
   *         \frac{a_j a_k}{ a_l} = K_a^{o,1} = \exp(\frac{\mu^o_l - \mu^o_j - \mu^o_k}{R T} )
   *   \f]
   *
   *  \f$  K_a^{o,1} \f$ is the dimensionless form of the equilibrium constant, associated with
   *  the pressure dependent standard states \f$ \mu^o_l(T,P) \f$ and their associated activities,
   *  \f$ a_l \f$, repeated here:
   *
   *       \f[
   *            \mu_l(T,P) = \mu^o_l(T, P) + R T \log(a_l)
   *       \f]
   *
   *  We can switch over to expressing the equilibrium constant in terms of the reference
   *  state chemical potentials
   *
   *   \f[
   *       K_a^{o,1} = \exp(\frac{\mu^{ref}_l - \mu^{ref}_j - \mu^{ref}_k}{R T} ) * \frac{P_{ref}}{P}
   *   \f]
   *
   *   The concentration equilibrium constant, \f$ K_c \f$, may be obtained by changing over
   *   to activity concentrations. When this is done:
   *
   *   \f[
   *         \frac{C^a_j C^a_k}{ C^a_l} = C^o K_a^{o,1} = K_c^1 =
   *             \exp(\frac{\mu^{ref}_l - \mu^{ref}_j - \mu^{ref}_k}{R T} ) * \frac{P_{ref}}{RT}
   *   \f]
   *
   *    %Kinetics managers will calculate the concentration equilibrium constant, \f$ K_c \f$,
   *    using the second and third part of the above expression as a definition for the concentration
   *    equilibrium constant.
   *
   *    For completeness, the pressure equilibrium constant may be obtained as well
   *
   *   \f[
   *         \frac{P_j P_k}{ P_l P_{ref}} = K_p^1 = \exp(\frac{\mu^{ref}_l - \mu^{ref}_j - \mu^{ref}_k}{R T} )
   *   \f]
   *
   *   \f$ K_p \f$ is the simplest form of the equilibrium constant for ideal gases. However, it isn't
   *   necessarily the simplest form of the equilibrium constant for other types of phases; \f$ K_c \f$ is
   *   used instead because it is completely general.
   *
   *   The reverse rate of progress may be written down as
   *   \f[
   *    R^{-1} = k^{-1} C_l^a =  k^{-1} (C^o a_l)
   *   \f]
   *
   *  where we can use the concept of microscopic reversibility to
   *  write the reverse rate constant in terms of the
   *  forward reate constant and the concentration equilibrium
   *   constant, \f$ K_c \f$.
   *
   *    \f[
   *       k^{-1} =  k^1 K^1_c
   *    \f]
   *
   *  \f$k^{-1} \f$ has units of s-1.
   *
   *
   * <HR>
   * <H2> Instantiation of the Class </H2>
   * <HR>
   *
   *
   * The constructor for this phase is located in the default ThermoFactory
   * for %Cantera. A new %IdealGasPhase may be created by the following code
   * snippet:
   *
   * @code
   *    XML_Node *xc = get_XML_File("silane.xml");
   *    XML_Node * const xs = xc->findNameID("phase", "silane");
   *    ThermoPhase *silane_tp = newPhase(*xs);
   *    IdealGasPhase *silaneGas = dynamic_cast <IdealGasPhase *>(silane_tp);
   * @endcode
   *
   * or by the following constructor:
   *
   * @code
   *    XML_Node *xc = get_XML_File("silane.xml");
   *    XML_Node * const xs = xc->findNameID("phase", "silane");
   *    IdealGasPhase *silaneGas = new IdealGasPhase(*xs);
   * @endcode
   *
   * <HR>
   * <H2> XML Example </H2>
   * <HR>
   *   An example of an XML Element named phase setting up a IdealGasPhase
   *   object named silane is given below.
   *
   *
   * @verbatim
   <!--     phase silane      -->
   <phase dim="3" id="silane">
   <elementArray datasrc="elements.xml"> Si  H  He </elementArray>
   <speciesArray datasrc="#species_data">
   H2  H  HE  SIH4  SI  SIH  SIH2  SIH3  H3SISIH  SI2H6
   H2SISIH2  SI3H8  SI2  SI3
   </speciesArray>
   <reactionArray datasrc="#reaction_data"/>
   <thermo model="IdealGas"/>
   <kinetics model="GasKinetics"/>
   <transport model="None"/>
   </phase>
   @endverbatim
   *
   *   The model attribute "IdealGas" of the thermo XML element identifies the phase as
   *   being of the type handled by the IdealGasPhase object.
   *
   *    @ingroup thermoprops
   *

  */
  class MargulesVPSSTP : public GibbsExcessVPSSTP {

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
    MargulesVPSSTP();

    //! Construct and initialize a MargulesVPSSTP ThermoPhase object
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
    MargulesVPSSTP(std::string inputFile, std::string id = "");

    //! Construct and initialize a MargulesVPSSTP ThermoPhase object
    //! directly from an XML database
    /*!
     *  @param phaseRef XML phase node containing the description of the phase
     *  @param id     id attribute containing the name of the phase.
     *                (default is the empty string)
     */
    MargulesVPSSTP(XML_Node& phaseRef, std::string id = "");


    //! Special constructor for a hard-coded problem
    /*!
     *
     *   LiKCl treating the PseudoBinary layer as passthrough.
     *   -> test to predict the eutectic and liquidus correctly.
     *
     */
    MargulesVPSSTP(int testProb);

    //! Copy constructor
    /*!
     *  Note this stuff will not work until the underlying phase
     *  has a working copy constructor
     *
     * @param b class to be copied
     */
    MargulesVPSSTP(const MargulesVPSSTP&b);

    //! Assignment operator
    /*!
     *
     * @param b class to be copied.
     */
    MargulesVPSSTP& operator=(const MargulesVPSSTP &b);

    //! Destructor
    virtual ~MargulesVPSSTP();

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

    
    //! Get the species electrochemical potentials.
    /*!
     * These are partial molar quantities.
     * This method adds a term \f$ Fz_k \phi_k \f$ to the 
     * to each chemical potential.
     *
     * Units: J/kmol
     *
     * @param mu     output vector containing the species electrochemical potentials.
     *               Length: m_kk.
     */
    void getElectrochemPotentials(doublereal* mu) const;

 
    //! Get the array of temperature derivatives of the log activity coefficients
    /*!
     * This function is a virtual class, but it first appears in GibbsExcessVPSSTP
     * class and derived classes from GibbsExcessVPSSTP.
     *
     *  units = 1/Kelvin
     *
     * @param dlnActCoeffdT    Output vector of temperature derivatives of the 
     *                         log Activity Coefficients. length = m_kk
     */
    virtual void getdlnActCoeffdT(doublereal *dlnActCoeffdT) const;

 
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


  private:
  
    //! Process an XML node called "binaryNeutralSpeciesParameters"
    /*!
     * This node contains all of the parameters necessary to describe
     * the Margules model for a particular binary interaction.
     * This function reads the XML file and writes the coefficients
     * it finds to an internal data structures.
     *
     * @param BinSalt  reference to the XML_Node named "binaryNeutralSpeciesParameters"
     *                 containing the binary interaction
     */
    void readXMLBinarySpecies(XML_Node &xmLBinarySpecies);

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

    // Update the derivative of the log of the activity coefficients wrt T
    /*
     * This function will be called to update the internally storred
     * natural logarithm of the activity coefficients
     *
     */
    void s_update_dlnActCoeff_dT() const;


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

    mutable vector_fp m_HE_b_ij;

    mutable vector_fp m_HE_c_ij;

    mutable vector_fp m_HE_d_ij;


    mutable vector_fp m_SE_b_ij;

    mutable vector_fp m_SE_c_ij;

    mutable vector_fp m_SE_d_ij;
    
    vector_int m_pSpecies_A_ij;
    vector_int m_pSpecies_B_ij;


    int formMargules_;
    int formTempModel_;

  private:

  
  };



}
        
#endif





