/**
 *  @file gibbsExcessVPSSTP.h
 *   Header for intermediate ThermoPhase object for phases which
 *   employ gibbs excess free energy based formulations
 *  (see \ref thermoprops 
 * and class \link Cantera::gibbsExcessVPSSTP gibbsExcessVPSSTP\endlink).
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
 *  $Id: GibbsExcessVPSSTP.h,v 1.3 2009/03/03 21:08:31 hkmoffa Exp $
 */

#ifndef CT_GIBBSEXCESSVPSSTP_H
#define CT_GIBBSEXCESSVPSSTP_H

#include "VPStandardStateTP.h"

namespace Cantera {

  /**
   * @ingroup thermoprops
   */

  /*!
   * GibbsExcessVPSSTP is a derived class of ThermoPhase that handles
   * variable pressure standard state methods for calculating
   * thermodynamic properties that are further based on
   * expressing the Excess Gibbs free energy as a function of
   * the mole fractions (or pseudo mole fractions) of consitituents.
   * This category is the workhorse for describing molten salts, 
   * solid-phase mixtures of semiconductors, and mixtures of miscible
   * and semi-miscible compounds.
   *
   * It includes 
   *   . regular solutions
   *   . Margueles expansions
   *   . NTRL equation
   *   . Wilson's equation
   *   . UNIQUAC equation of state.
   *
   * This class adds additional functions onto the %ThermoPhase interface
   * that handles the calculation of the excess Gibbs free energy. The %ThermoPhase
   * class includes a member function, ThermoPhase::activityConvention() 
   * that indicates which convention the activities are based on. The
   * default is to assume activities are based on the molar convention.
   * That default is used here. 
   *
   * All of the Excess Gibbs free energy formulations in this area employ
   * symmetrical formulations. 
   *
   *
   *  Chemical potentials
   * of species k, \f$ \mu_o \f$, has the following general format:
   *
   * \f[
   *    \mu_k = \mu^o_k(T,P) + R T ln( \gamma_k X_k ) 
   * \f]
   *
   *
   * where \f$ \gamma_k^{\triangle} \f$ is a molar based activity coefficient for species
   * \f$k\f$.
   * 
   * GibbsExcessVPSSTP contains an internal vector with the current mole
   * fraction vector. That's one of its primary usages.
   *
   *  <H3> SetState Strategy  </H3>
   *
   *   The gibbsExcessVPSSTP object does not have a setState strategy.
   *   It's strictly an interfacial layer that writes the current mole fractions to the
   *   State object.
   *
   *
   */
  class GibbsExcessVPSSTP : public VPStandardStateTP  {

  public:
        
    /// Constructors 
    /*!
     * This doesn't do much more than initialize constants with
     * default values for water at 25C. Water molecular weight 
     * comes from the default elements.xml file. It actually
     * differs slightly from the IAPWS95 value of 18.015268. However,
     * density conservation and therefore element conservation
     * is the more important principle to follow.
     */
    GibbsExcessVPSSTP();

    //! Copy constructor
    /*!
     *  Note this stuff will not work until the underlying phase
     *  has a working copy constructor
     *
     * @param b class to be copied
     */
    GibbsExcessVPSSTP(const GibbsExcessVPSSTP &b);

    /// Assignment operator
    /*!
     *
     * @param b class to be copied.
     */
    GibbsExcessVPSSTP& operator=(const GibbsExcessVPSSTP &b);

    /// Destructor. 
    virtual ~GibbsExcessVPSSTP();

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

   

    /**
     * @} 
     * @name  Molar Thermodynamic Properties 
     * @{
     */


   
 

    /**
     * @}
     * @name Mechanical Properties
     * @{
     */

    //! Set the internally storred pressure (Pa) at constant
    //! temperature and composition
    /*!
     *  This method sets the pressure within the object.
     *  The water model is a completely compressible model.
     *  Also, the dielectric constant is pressure dependent.
     *
     *  @param p input Pressure (Pa)
     *
     * @todo Implement a variable pressure capability
     */
    virtual void setPressure(doublereal p);

  private:
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
     *
     * NOTE: This is a non-virtual function, which is not a 
     *       member of the ThermoPhase base class. 
     */
    void calcDensity();

  public:
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

    /**
     * Returns the units of the standard and generalized
     * concentrations Note they have the same units, as their
     * ratio is defined to be equal to the activity of the kth
     * species in the solution, which is unitless.
     *
     * This routine is used in print out applications where the
     * units are needed. Usually, MKS units are assumed throughout
     * the program and in the XML input files.
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

    
    //! Get the array of non-dimensional activities (molality
    //! based for this class and classes that derive from it) at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     * All standard state properties for molality-based phases are
     * evaluated consistent with the molality scale. Therefore, this function
     * must return molality-based activities.
     *
     * \f[
     *  a_i^\triangle = \gamma_k^{\triangle} \frac{m_k}{m^\triangle}
     * \f]
     *
     * This function must be implemented in derived classes.
     *
     * @param ac     Output vector of molality-based activities. Length: m_kk.
     */
    virtual void getActivities(doublereal* ac) const;

   
 
    //@}
    /// @name  Partial Molar Properties of the Solution 
    //@{


    /**
     * Get the species electrochemical potentials. 
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
  

    //! Initialize lengths of local variables after all species have
    //! been identified.
    void initLengths();
            


  private:
    //! Error function
    /*!
     *  Print an error string and exit
     *
     * @param msg  Message to be printed
     */
    doublereal err(std::string msg) const;

  protected:

    double checkMFSum(const doublereal * const x) const;

  protected:

    //! Storage for the current values of the mole fractions of the species
    mutable std::vector<doublereal> moleFractions_;

    //! Storage for the current values of the activity coefficients of the
    //! species, divided by RT
    mutable std::vector<doublereal> lnActCoeff_Scaled_;

    mutable std::vector<doublereal> m_pp;

  };


}
        
#endif





