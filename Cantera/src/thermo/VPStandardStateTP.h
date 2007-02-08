/**
 *  @file VPStandardStateTP.h
 *
 * Header file for a derived class of ThermoPhase that handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties. These include most of the
 * methods for calculating liquid electrolyte thermodynamics.
 */
/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 *  $Author$
 *  $Date$
 *  $Revision$
 */

#ifndef CT_VPSTANDARDSTATETP_H
#define CT_VPSTANDARDSTATETP_H

#include "ThermoPhase.h"

namespace Cantera {

  class XML_Node;

  /**
   * @ingroup thermoprops
   *
   *  This is a filter class for ThermoPhase that implements some prepatory
   *  steps for efficiently handling
   *  a variable pressure standard state for species.
   *
   *  Several concepts are introduced. The first concept is there are temporary
   *  variables for holding the species standard values of Cp, H, S, and V at the
   *  last temperature and pressure called. These functions are not recalculated
   *  if a new call is made using the previous temperature and pressure.
   *
   *  There are also temporary
   *  variables for holding the species reference-state values of Cp, H, S, and V at the
   *  last temperature and reference pressure called. These functions are not recalculated
   *  if a new call is made using the previous temperature. 
   *   
   *  @nosubgrouping
   */
  class VPStandardStateTP : public ThermoPhase {

  public:

    /*!
     *   
     * @name Constructors and Duplicators for %VPStandardStateTP 
     *
     */   
    /// Constructor. 
    VPStandardStateTP();

    /// Copy Constructor.
    VPStandardStateTP(const VPStandardStateTP &);

    /// Assignment operator
    VPStandardStateTP& operator=(const VPStandardStateTP &);

    /// Destructor. 
    virtual ~VPStandardStateTP();

    /*
     * Duplication routine
     */
    virtual ThermoPhase *duplMyselfAsThermoPhase();

    //@}

    /**
     * @name  Utilities (VPStandardStateTP)
     */
    //@{
    /** 
     * Equation of state type flag. The base class returns
     * zero. Subclasses should define this to return a unique
     * non-zero value. Constants defined for this purpose are
     * listed in mix_defs.h.
     */
    virtual int eosType() const { return 0; }

    //@}
 

    /// @name  Partial Molar Properties of the Solution  (VPStandardStateTP)
    //@{

    /**
     * Get the array of non-dimensional species chemical potentials
     * These are partial molar Gibbs free energies.
     * \f$ \mu_k / \hat R T \f$.
     * Units: unitless
     *
     * We close the loop on this function, here, calling
     * getChemPotentials() and then dividing by RT.
     *
     * @param mu    Output vector of  non-dimensional species chemical potentials
     *              Length: m_kk.
     */
    virtual void getChemPotentials_RT(doublereal* mu) const;
  
    //@}

    /*!
     * @name  Properties of the Standard State of the Species in the Solution  (VPStandardStateTP)
     *
     *  Within VPStandardStateTP, these properties are calculated via a common routine, _updateStandardStateThermo(),
     *  which must be overloaded in inherited objects.
     *  The values are cached within this object, and are not recalculated unless
     *  the temperature or pressure changes.
     */
    //@{
    
    //!Get the array of chemical potentials at unit activity.
    /*!
     * These are the standard state chemical potentials \f$ \mu^0_k(T,P)
     * \f$. The values are evaluated at the current temperature and pressure.
     *
     * @param mu   Output vector of standard state chemical potentials.
     *             length = m_kk. units are J / kmol.
     */
    virtual void getStandardChemPotentials(doublereal* mu) const;

    /**
     * Get the nondimensional Enthalpy functions for the species
     * at their standard states at the current
     * <I>T</I> and <I>P</I> of the solution.
     *
     * @param hrt     Output vector of standard state enthalpies.
     *                length = m_kk. units are unitless.
     */
    virtual void getEnthalpy_RT(doublereal* hrt) const;

    /**
     * Get the array of nondimensional Enthalpy functions for the
     * standard state species
     * at the current <I>T</I> and <I>P</I> of the solution.
     *
     * @param sr     Output vector of nondimensional standard state
     *               entropies. length = m_kk.
     */
    virtual void getEntropy_R(doublereal* sr) const;

    /**
     * Get the nondimensional Gibbs functions for the species
     * at their standard states of solution at the current T and P
     * of the solution.
     *
     * @param grt    Output vector of nondimensional standard state
     *               Gibbs free energies. length = m_kk.
     */
    virtual void getGibbs_RT(doublereal* grt) const;

    /**
     * Get the nondimensional Gibbs functions for the standard
     * state of the species at the current T and P.
     *
     * @param gpure  Output vector of standard state
     *               Gibbs free energies. length = m_kk.
     *               units are J/kmol.
     */
    virtual void getPureGibbs(doublereal* gpure) const;

    /**
     *  Returns the vector of nondimensional
     *  internal Energies of the standard state at the current temperature
     *  and pressure of the solution for each species.
     * \f[
     *  u^{ss}_k(T,P) = h^{ss}_k(T)  - P * V^{ss}_k
     * \f]
     *
     * @param urt    Output vector of nondimensional standard state
     *               internal energies. length = m_kk.
     */
    virtual void getIntEnergy_RT(doublereal *urt) const;

    /**
     * Get the nondimensional Heat Capacities at constant
     * pressure for the standard state of the species 
     * at the current T and P. 
     *
     * This is redefined here to call the internal function,  _updateStandardStateThermo(),
     * which calculates all standard state properties at the same time.
     *
     * @param cpr    Output vector containing the 
     *               the nondimensional Heat Capacities at constant
     *               pressure for the standard state of the species.
     *               Length: m_kk. 
     */
    virtual void getCp_R(doublereal* cpr) const;

    /**
     * Get the molar volumes of each species in their standard
     * states at the current
     * <I>T</I> and <I>P</I> of the solution.
     * units = m^3 / kmol
     *
     * This is redefined here to call the internal function,  _updateStandardStateThermo(),
     * which calculates all standard state properties at the same time.
     *
     * @param vol Output vector of species volumes. length = m_kk.
     *            units =  m^3 / kmol
     */
    virtual void getStandardVolumes(doublereal *vol) const;

  protected:

    //! Updates the standard state thermodynamic functions at the current T and P of the solution.
    /*!
     * @internal
     *
     * This function gets called for every call to functions in this
     * class. It checks to see whether the temperature or pressure has changed and
     * thus the ss thermodynamics functions for all of the species
     * must be recalculated.
     *
     * This function is responsible for updating the following internal members:
     *
     *    m_hss_RT;
     *    m_cpss_R;
     *    m_gss_RT;
     *    m_sss_R;
     *    m_Vss
     *
     *  Note, this will throw an error. It must be reimplemented in derived classes.
     */                    
    virtual void _updateStandardStateThermo() const;

  public:
    //@}
    /// @name Thermodynamic Values for the Species Reference States (VPStandardStateTP)
    /*!
     *  There are also temporary
     *  variables for holding the species reference-state values of Cp, H, S, and V at the
     *  last temperature and reference pressure called. These functions are not recalculated
     *  if a new call is made using the previous temperature.
     *  All calculations are done within the routine  _updateRefStateThermo().
     */
    //@{

    /*!
     *  Returns the vector of nondimensional
     *  enthalpies of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *
     * @param hrt Output vector contains the nondimensional enthalpies
     *            of the reference state of the species
     *            length = m_kk, units = dimensionless.
     */
    virtual void getEnthalpy_RT_ref(doublereal *hrt) const;
     
    /*!
     *  Returns the vector of nondimensional
     *  Gibbs free energies of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *
     * @param grt Output vector contains the nondimensional Gibbs free energies
     *            of the reference state of the species
     *            length = m_kk, units = dimensionless.
     */
    virtual void getGibbs_RT_ref(doublereal *grt) const;
                   
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
    virtual void getGibbs_ref(doublereal *g) const;
      
    /*!
     *  Returns the vector of nondimensional
     *  entropies of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *
     * @param er  Output vector contain the nondimensional entropies
     *            of the species in their reference states
     *            length: m_kk, units: dimensionless.
     */
    virtual void getEntropy_R_ref(doublereal *er) const;
                 
    /**
     *  Returns the vector of nondimensional
     *  constant pressure heat capacities of the reference state
     *  at the current temperature of the solution
     *  and reference pressure for the species.
     *
     * @param cprt Output vector contains the nondimensional heat capacities
     *             of the species in their reference states
     *             length: m_kk, units: dimensionless.
     */
    virtual void getCp_R_ref(doublereal *cprt) const;

    //! Recalculate the Reference state thermo functions
    /*!
     * This function checks to see whether the temperature has changed and
     * thus the reference thermodynamics functions for all of the species
     * must be recalculated.
     * If the temperature has changed, the species thermo manager is called
     * to recalculate G, Cp, H, and S at the current temperature and at
     * the reference pressure.
     */

  protected:

    //! Recalculate the Reference state thermo functions
    /*!
     * This function checks to see whether the temperature has changed and
     * thus the reference thermodynamics functions for all of the species
     * must be recalculated.
     * If the temperature has changed, the species thermo manager is called
     * to recalculate G, Cp, H, and S at the current temperature and at
     * the reference pressure.
     */
    virtual void _updateRefStateThermo() const;
 
    //@}

	
  public:
 
    //! @name Initialization Methods - For Internal use (VPStandardState)
    /*!
     * The following methods are used in the process of constructing
     * the phase and setting its parameters from a specification in an 
     * input file. They are not normally used in application programs.
     * To see how they are used, see files importCTML.cpp and 
     * ThermoFactory.cpp.
     */
    //@{

    /**
     * Set equation of state parameter values from XML
     * entries. This method is called by function importPhase in
     * file importCTML.cpp when processing a phase definition in
     * an input file. It should be overloaded in subclasses to set
     * any parameters that are specific to that particular phase
     * model. 
     *   
     * @param eosdata An XML_Node object corresponding to
     *                the "thermo" entry for this phase in the input file.
     */
    virtual void setParametersFromXML(const XML_Node& eosdata) {}
  
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
     *
     * @see importCTML.cpp
     */
    virtual void initThermo();

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
    virtual void initThermoXML(XML_Node& phaseNode, std::string id);

  private:
    //!  @internal Initialize the internal lengths in this object.
    /*!
     * Note this is not a virtual function.
     */
    void initLengths();

   //@}

  protected:
    
    //! The last temperature at which the reference thermodynamic properties were calculated at.
    mutable doublereal    m_tlast;

    //! The last pressure at which the Standard State thermodynamic properties were calculated at.
    mutable doublereal    m_plast;

    /*!
     * Vector containing the species reference enthalpies at T = m_tlast
     * and P = p_ref.
     */
    mutable vector_fp      m_h0_RT;

    /**
     * Vector containing the species reference constant pressure
     * heat capacities at T = m_tlast    and P = p_ref.
     */
    mutable vector_fp      m_cp0_R;

    /**
     * Vector containing the species reference Gibbs functions
     * at T = m_tlast  and P = p_ref.
     */
    mutable vector_fp      m_g0_RT;

    /**
     * Vector containing the species reference entropies
     * at T = m_tlast and P = p_ref.
     */
    mutable vector_fp      m_s0_R;

    /**
     * Vector containing the species reference volumes
     * at T = m_tlast and P = p_ref
     */   
    mutable vector_fp      m_V0;

    /**
     * Vector containing the species Standard State enthalpies at T = m_tlast
     * and P = m_plast.
     */
    mutable vector_fp      m_hss_RT;

    /**
     * Vector containing the species Standard State constant pressure
     * heat capacities at T = m_tlast and P = m_plast.
     */
    mutable vector_fp      m_cpss_R;

    /**
     * Vector containing the species Standard State Gibbs functions
     * at T = m_tlast and P = m_plast.
     */
    mutable vector_fp      m_gss_RT;

    /**
     * Vector containing the species Standard State entropies
     * at T = m_tlast and P = m_plast.
     */
    mutable vector_fp      m_sss_R;

    /**
     * Vector containing the species standard state volumes
     * at T = m_tlast and P = m_plast
     */   
    mutable vector_fp      m_Vss;
  
  private:

    /*!
     * VPStandardStateTP has its own err routine
     */
    doublereal err(std::string msg) const;

  };
}
        
#endif





  
