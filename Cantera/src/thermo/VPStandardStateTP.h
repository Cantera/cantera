/**
 *  @file VPStandardStateTP.h
 * Header file for a derived class of ThermoPhase that handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties (see \ref thermoprops and
 * class \link Cantera::VPStandardStateTP VPStandardStateTP\endlink).
 *
 * These include most of the
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
   *  variables for holding the species standard state values 
   *  of Cp, H, S, G, and V at the
   *  last temperature and pressure called. These functions are not recalculated
   *  if a new call is made using the previous temperature and pressure.
   *
   *  There are also temporary
   *  variables for holding the species reference-state values of Cp, H, S, and G at the
   *  last temperature and reference pressure called. These functions are not recalculated
   *  if a new call is made using the previous temperature. 
   *
   *  To support the above functionality, pressure and temperature variables,
   *  m_plast and m_tlast, are kept which store the last pressure and temperature
   *  used in the evaluation of standard state properties. An optional utility is provided
   *  to store the results from the last temperature and pressure standard
   *  state calculation and use it on subsequent calculations, if the temperature
   *  and pressure are unchanged.
   *
   *  If #m_useTmpRefStateStorage is set to true, then the following internal
   *  arrays, containing information about the reference arrays,
   *  are calculated and kept up to date at every call.
   *
   *  - #m_h0_RT
   *  - #m_g0_RT
   *  - #m_s0_R
   *  - #m_cp0_R
   * 
   *  The virtual function #_updateRefStateThermo() is supplied to do this
   *  and may be reimplemented in child routines. A default implementation
   *  based on the speciesThermo class is supplied in this base class.
   *  #_updateStandardStateThermo() is called whenever a reference state property is needed.
   *
   *  When  #m_useTmpStandardStateStorage is true, then the following
   *  internal arrays, containing information on the standard state properties
   *  are calculated and kept up to date. 
   *
   *  -  #m_hss_RT;
   *  -  #m_cpss_R;
   *  -  #m_gss_RT;
   *  -  #m_sss_R;
   *  -  #m_Vss
   *
   *  The virtual function #_updateStandardStateThermo() is supplied to do this
   *  and must be reimplemented in child routines, when  #m_useTmpStandardStateStorage is true.
   *  It may be optionally reimplemented in child routines if
   *  #m_useTmpStandardStateStorage is false.
   *  #_updateStandardStateThermo() is called whenever a standard state property is needed.
   *
   *  This class is usually used for nearly incompressible phases. For those phases, it
   *  makes sense to change the equation of state independent variable from density to pressure.
   *
   * @todo
   *   Put some teeth into this level by overloading the setDensity() function. It should
   *   now throw an exception. Instead, setPressure routines should calculate the
   *   solution density and then call State:setDensity() directly.
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
     * @name  Properties of the Standard State of the Species in the Solution  (VPStandardStateTP)
     *
     *  Within VPStandardStateTP, these properties are calculated via a common routine, 
     *  _updateStandardStateThermo(),
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

   
    //! Get the nondimensional Gibbs functions for the standard
    //! state of the species at the current T and P.
    /*!
     *  (Note resolved at this level)
     *
     * @param gpure  Output vector of standard state
     *               Gibbs free energies. length = m_kk.
     *               units are J/kmol.
     *
     * @todo This could be eliminated. It doesn't fit into the current
     *       naming convention.
     */
    void getPureGibbs(doublereal* gpure) const;

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
     * If m_useTmpStandardStateStorage is true,
     * this function must be called for every call to functions in this
     * class. It checks to see whether the temperature or pressure has changed and
     * thus the ss thermodynamics functions for all of the species
     * must be recalculated.
     *
     * This function is responsible for updating the following internal members,
     * when  m_useTmpStandardStateStorage is true.
     *
     *  -  m_hss_RT;
     *  -  m_cpss_R;
     *  -  m_gss_RT;
     *  -  m_sss_R;
     *  -  m_Vss
     *
     *  If m_useTmpStandardStateStorage is not true, this function may be
     *  required to be called by child classes to update internal member data.
     *
     *  Note, this will throw an error. It must be reimplemented in derived classes.
     *
     * @param pres Pressure at which to carry out the calculation.
     *             The default is to use the current pressure, storred in m_Pcurrent.
     */                    
    virtual void _updateStandardStateThermo(doublereal pres = -1.0) const;

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
    virtual void getCp_R_ref(doublereal *cprt) const;

    //!  Get the molar volumes of the species reference states at the current
    //!  <I>T</I> and <I>P_ref</I> of the solution.
    /*!
     * units = m^3 / kmol
     *
     * @param vol     Output vector containing the standard state volumes.
     *                Length: m_kk.
     */
    virtual void getStandardVolumes_ref(doublereal *vol) const;

  protected:

    //! Recalculate the Reference state thermo functions
    /*!
     * This function checks to see whether the temperature has changed and
     * thus the reference thermodynamics functions for all of the species
     * must be recalculated.
     * It must be called for every reference state function evaluation,
     * if m_useTmpRefStateStorage is set to true.
     * If the temperature has changed, the species thermo manager is called
     * to recalculate the following internal arrays at the current temperature and at
     * the reference pressure:
     *
     *  - m_h0_RT
     *  - m_g0_RT
     *  - m_s0_R
     *  - m_cp0_R
     *
     * This function may be reimplemented in child objects. However, it doesn't
     * necessarily have to be, if the species thermo manager can carry
     * out the full calculation.
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

    //! The current pressure of the solution (Pa)
    /*!
     * It gets initialized to 1 atm.
     */
    mutable doublereal    m_Pcurrent;
    
    //! The last temperature at which the reference thermodynamic properties were calculated at.
    mutable doublereal    m_tlast;

    //! The last temperature at which the reference thermodynamic properties were calculated at.
    mutable doublereal    m_tlast_ref;

    //! The last pressure at which the Standard State thermodynamic properties were calculated at.
    mutable doublereal    m_plast;

    /*!
     * Reference pressure (Pa) must be the same for all species
     * - defaults to 1 atm.
     */
    doublereal m_p0;

    /*!
     * boolean indicating whether temporary reference state storage is used
     * -> default is false
     */
    bool m_useTmpRefStateStorage;

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

    /*!
     * boolean indicating whether temporary standard state storage is used
     * -> default is false
     */
    bool m_useTmpStandardStateStorage;

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

    //! VPStandardStateTP has its own err routine
    /*!
     * @param msg  Error message string
     */
    doublereal err(std::string msg) const;

  };
}
        
#endif





  
