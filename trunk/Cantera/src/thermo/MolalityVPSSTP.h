/**
 *  @file MolalityVPSSTP.h
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
 *  $Id: MolalityVPSSTP.h,v 1.6 2006/07/06 15:39:34 hkmoffa Exp $
 */

#ifndef CT_MOLALITYVPSSTP_H
#define CT_MOLALITYVPSSTP_H

#include "VPStandardStateTP.h"

namespace Cantera {

  /**
   * @ingroup thermoprops
   */

  /**
   * MolalityVPSSTP is a derived class of ThermoPhase that handles
   * variable pressure standard state methods for calculating
   * thermodynamic properties that are further based on
   * molality-scaled activities.
   * These include most of the methods
   * for calculating liquid electrolyte thermodynamics.
   */
  class MolalityVPSSTP : public VPStandardStateTP  {

  public:
        
    /// Constructors 
    MolalityVPSSTP();
    MolalityVPSSTP(const MolalityVPSSTP &);
    /// Assignment operator
    MolalityVPSSTP& operator=(const	MolalityVPSSTP&);

    /// Destructor. 
    virtual ~MolalityVPSSTP();

    /**
     * Duplication routine for objects which inherit from 
     * ThermoPhase.
     *
     *  This virtual routine can be used to duplicate thermophase objects
     *  inherited from ThermoPhase even if the application only has
     *  a pointer to ThermoPhase to work with.
     */
    virtual ThermoPhase *duplMyselfAsThermoPhase();

    /**
     *   
     * @name  Utilities  
     * @{
     */

    /** 
     * Equation of state type flag. The ThermoPhase base class returns
     * zero. Subclasses should define this to return a unique
     * non-zero value. Known constants defined for this purpose are
     * listed in mix_defs.h. The MolalityVPSSTP class also returns
     * zero, as it is a non-complete class.
     */
    virtual int eosType() const { return 0; }


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
     * This routine sets the index number of the solvent for
     * the phase.
     *  
     *  Note, having a solvent 
     *   is a precursor to many things having to do with molality.
     *
     * @param k the solvent index number
     */
    void setSolvent(int k);

    /**
     * Sets the minimum mole fraction in the molality formulation.
     * Note the molality formulation is singular in the limit that
     * the solvent mole fraction goes to zero. Numerically, how
     * this limit is treated and resolved is an ongoing issue within
     * Cantera.
     */
    void setMoleFSolventMin(doublereal xmolSolventMIN);

    /**
     * Returns the solvent index.
     */
    int solventIndex() const;

    /**
     * Returns the minimum mole fraction in the molality 
     * formulation.
     */
    doublereal moleFSolventMin() const;

    /**
     * Calculates the molality of all species and
     * stores the result internally.
     */
    void calcMolalities() const;

    /**
     * getMolalities()
     *   This function will return the molalities of the
     *   species.
     *
     */
    void getMolalities(doublereal * const molal) const;


    void setMolalities(const doublereal * const molal);
    void setMolalitiesByName(compositionMap& xMap);
    void setMolalitiesByName(const string &);

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
     * This method returns the activity convention.
     * Currently, there are two activity conventions
     *  Molar-based activities
     *       Unit activity of species at either a hypothetical pure
     *       solution of the species or at a hypothetical
     *       pure ideal solution at infinite dilution
     *   cAC_CONVENTION_MOLAR 0
     *      - default
     *  
     *  Molality based acvtivities
     *       (unit activity of solutes at a hypothetical 1 molal
     *        solution referenced to infinite dilution at all
     *        pressures and temperatures).
     *   cAC_CONVENTION_MOLALITY 1
     *
     *  We set the convention to molality here.
     */
    int activityConvention() const;

    /**
     * This method returns an array of generalized concentrations
     * \f$ C_k\f$ that are defined such that 
     * \f$ a_k = C_k / C^0_k, \f$ where \f$ C^0_k \f$ 
     * is a standard concentration
     * defined below.  These generalized concentrations are used
     * by kinetics manager classes to compute the forward and
     * reverse rates of elementary reactions. 
     *
     * @param c Array of generalized concentrations. The 
     *           units depend upon the implementation of the
     *           reaction rate expressions within the phase.
     */
    virtual void getActivityConcentrations(doublereal* c) const {
      err("getActivityConcentrations");
    }

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
     */
    virtual doublereal standardConcentration(int k=0) const {
      err("standardConcentration");
      return -1.0;
    }

    /**
     * Returns the natural logarithm of the standard 
     * concentration of the kth species
     */
    virtual doublereal logStandardConc(int k=0) const {
      err("logStandardConc");
      return -1.0;
    }

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
     *  uA[0] = kmol units - default  = 1
     *  uA[1] = m    units - default  = -nDim(), the number of spatial
     *                                dimensions in the Phase class.
     *  uA[2] = kg   units - default  = 0;
     *  uA[3] = Pa(pressure) units - default = 0;
     *  uA[4] = Temperature units - default = 0;
     *  uA[5] = time units - default = 0
     */
    virtual void getUnitsStandardConc(double *uA, int k = 0,
				      int sizeUA = 6);

    /**
     * Get the array of non-dimensional activities (molality
     * based for this class and classes that derive from it) at
     * the current solution temperature, pressure, and
     * solution concentration.
     */
    virtual void getActivities(doublereal* ac) const {
      err("getActivities");
    }

    /**
     * Get the array of non-dimensional activity coefficients at
     * the current solution temperature, pressure, and
     * solution concentration.
     * These are mole fraction based activity coefficients. In this
     * object, their calculation is based on translating the values
     * of Molality based activity coefficients.
     *  See Denbigh p. 278 for a thorough discussion
     */
    void getActivityCoefficients(doublereal* ac) const;

    /**
     * Get the array of non-dimensional molality based 
     * activity coefficients at the current solution temperature, 
     * pressure, and  solution concentration.
     *  See Denbigh p. 278 for a thorough discussion
     */
    virtual void getMolalityActivityCoefficients(doublereal *acMolality)
      const {
      err("getMolalityActivityCoefficients");
    }

    /**
     * Calculate the osmotic coefficient
     *  units = dimensionless
     */
    virtual double osmoticCoefficient() const;

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
     */
    void getElectrochemPotentials(doublereal* mu) const {
      getChemPotentials(mu);
      double ve = Faraday * electricPotential();
      for (int k = 0; k < m_kk; k++) {
	mu[k] += ve*charge(k);
      }
    }

 
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
     * This method is used by the ChemEquil element-potential
     * based equilibrium solver.
     * It sets the state such that the chemical potentials of the 
     * species within the current phase satisfy
     * \f[ \frac{\mu_k}{\hat R T} = \sum_m A_{k,m}
     * \left(\frac{\lambda_m} {\hat R T}\right) \f] where 
     * \f$ \lambda_m \f$ is the element potential of element m. The
     * temperature is unchanged.  Any phase (ideal or not) that
     * implements this method can be equilibrated by ChemEquil.
     */
    virtual void setToEquilState(const doublereal* lambda_RT) {
      err("setToEquilState");
    }

    // called by function 'equilibrate' in ChemEquil.h to transfer
    // the element potentials to this object
    void setElementPotentials(const vector_fp& lambda) {
      m_lambda = lambda;
    }

    void getElementPotentials(doublereal* lambda) {
      copy(m_lambda.begin(), m_lambda.end(), lambda);
    }

    //@}


    /**
     * Set equation of state parameter values from XML
     * entries. This method is called by function importPhase in
     * file importCTML.cpp when processing a phase definition in
     * an input file. It should be overloaded in subclasses to set
     * any parameters that are specific to that particular phase
     * model. 
     *
     * The MolalityVPSSTP object defines a new method for setting
     * the concentrations of a phase. The new method is defined by a
     * block called "soluteMolalities". If this block
     * is found, the concentrations within that phase are
     * set to the "name":"molalities pairs found within that
     * XML block. The solvent concentration is then set
     * to everything else.
     *   
     * @param state An XML_Node object corresponding to
     * the "state" entry for this phase in the input file.
     *
     */
    virtual void setStateFromXML(const XML_Node& state);

    /// The following methods are used in the process of constructing
    /// the phase and setting its parameters from a specification in an 
    /// input file. They are not normally used in application programs.
    /// To see how they are used, see files importCTML.cpp and 
    /// ThermoFactory.cpp.

    /**
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
    void initThermoXML(XML_Node& phaseNode, string id);

    /** 
     * Set the temperature (K), pressure (Pa), and molalities
     * (gmol kg-1) of the solutes
     */
    void setState_TPM(doublereal t, doublereal p, 
		      const doublereal * const molalities);

    /** Set the temperature (K), pressure (Pa), and molalities.  */
    void setState_TPM(doublereal t, doublereal p, compositionMap& m);

    /** Set the temperature (K), pressure (Pa), and molalities.  */
    void setState_TPM(doublereal t, doublereal p, const string& m);

  private:
    void initLengths();
            
  protected:

    int        m_indexSolvent;
    doublereal m_weightSolvent;
    /*
     * In any molality implementation, it makes sense to have
     * a minimum solvent mole fraction requirement, since the 
     * implementation becomes singular in the xmolSolvent=0
     * limit. The default is to set it to 0.01.
     * We then modify the molality definition to ensure that
     * molal_solvent = 0 when xmol_solvent = 0. 
     */
    doublereal m_xmolSolventMIN;
    /*
     * This is the multiplication factor that goes inside 
     * log expressions involving the molalities of species.
     * Its equal to Wt_0 / 1000.
     *     where Wt_0 = weight of solvent (kg/kmol)
     */
    doublereal m_Mnaught;
    
    mutable vector_fp  m_molalities;
  private:
    doublereal err(string msg) const;

  };

}
        
#endif





