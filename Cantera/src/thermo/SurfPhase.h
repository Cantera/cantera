/**
 *  @file SurfPhase.h
 *  Header for a simple thermoydnamics model of a surface phase derived from ThermoPhase, 
 *  assuming an ideal solution model
 *  (see \ref thermoprops and class \link Cantera::SurfPhase SurfPhase\endlink).
 */

/*  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2002 California Institute of Technology
 *
 */


#ifndef CT_SURFPHASE_H
#define CT_SURFPHASE_H

#include "mix_defs.h"
#include "ThermoPhase.h"


namespace Cantera {


    
  //! A simple thermoydnamics model for a surface phase, assuming an ideal solution model.
  /*!
   * The surface consists of a grid of equivalent sites. Surface species may be defined to
   * occupy one or more sites. The surface species are assumed to be
   * independent, and thus the species form an ideal solution.
   *
   * The density of surface sites is given by the variable \f$ n_0 \f$, which has MKS units
   * of kmol m-2. 
   *
   *
   * <b> Specification of Species Standard State Properties </b>
   * 
   *  It is assumed that the reference state thermodynamics may be
   *  obtained by a pointer to a populated species thermodynamic property
   *  manager class (see ThermoPhase::m_spthermo). How to relate pressure
   *  changes to the reference state thermodynamics is resolved at this level.
   *  
   *  Pressure is defined as an independent variable in this phase. However, it has
   *  no effect on any quantities, as the molar concentration is a constant.
   *
   * Therefore, The standard state internal energy for species  <I>k</I> is 
   * equal to the enthalpy for species <I>k</I>.
   *
   *       \f[
   *            u^o_k = h^o_k  
   *       \f]
   *  
   *   Also, the standard state chemical potentials, entropy, and heat capacities
   *   are independent of pressure. The standard state gibbs free energy is obtained
   * from the enthalpy and entropy functions.
   *   
   * <b> Specification of Solution Thermodynamic Properties </b>
   *
   * The activity of species defined in the phase is given by
   *       \f[
   *            a_k = \theta_k      
   *       \f]
   *
   * The chemical potential for species <I>k</I> is equal to 
   *       \f[
   *            \mu_k(T,P) = \mu^o_k(T) + R T \log(\theta_k)     
   *       \f]
   *
   * Pressure is defined as an independent variable in this phase. However, it has
   * no effect on any quantities, as the molar concentration is a constant.
   *
   * The internal energy for species k is equal to the enthalpy for species <I>k</I>
   *       \f[
   *            u_k = h_k  
   *       \f]
   *
   * The entropy for the phase is given by the following relation, which is
   * independent of the pressure:
   *
   *       \f[
   *            s_k(T,P) = s^o_k(T) - R \log(\theta_k)     
   *       \f]
   *
   * <b> Application within %Kinetics Managers </b>
   *
   * The activity concentration,\f$  C^a_k \f$, used by the kinetics manager, is equal to 
   * the actual concentration, \f$ C^s_k \f$, and is given by the following
   * expression.
   *       \f[
   *            C^a_k = C^s_k = \frac{\theta_k  n_0}{s_k}      
   *       \f]
   *
   * The standard concentration for species <I>k</I> is: 
   *        \f[
   *            C^0_k = \frac{n_0}{s_k}      
   *        \f]
   *
   * <b> Instantiation of the Class </b>
   *
   * The constructor for this phase is located in the default ThermoFactory
   * for Cantera. A new SurfPhase may be created by the following code snippet:
   *
   * @code
   *    XML_Node *xc = get_XML_File("diamond.xml"); 
   *    XML_Node * const xs = xc->findNameID("phase", "diamond_100");
   *    ThermoPhase *diamond100TP_tp = newPhase(*xs);
   *    SurfPhase *diamond100TP = dynamic_cast <SurfPhase *>(diamond100TP_tp);
   * @endcode
   *
   * or by the following constructor:
   *
   * @code
   *    XML_Node *xc = get_XML_File("diamond.xml"); 
   *    XML_Node * const xs = xc->findNameID("phase", "diamond_100");
   *    SurfPhase *diamond100TP = new SurfPhase(*xs);
   * @endcode
   *
   *   <b> XML Example </b>
   *
   *   An example of an XML Element named phase setting up a SurfPhase object named diamond_100
   *   is given below.
   *
   *  @verbatim
   * <phase dim="2" id="diamond_100">
   *    <elementArray datasrc="elements.xml">H C</elementArray>
   *    <speciesArray datasrc="#species_data">c6HH c6H* c6*H c6** c6HM c6HM* c6*M c6B </speciesArray>
   *    <reactionArray datasrc="#reaction_data"/>
   *    <state>
   *       <temperature units="K">1200.0</temperature>
   *       <coverages>c6H*:0.1, c6HH:0.9</coverages>
   *    </state>
   *    <thermo model="Surface">
   *       <site_density units="mol/cm2">3e-09</site_density>
   *    </thermo>
   *    <kinetics model="Interface"/>
   *    <transport model="None"/>
   *    <phaseArray> 
   *         gas_phase diamond_bulk 
   *    </phaseArray>
   * </phase>
   *
   *  @endverbatim
   *
   * The model attribute, "Surface", on the thermo element identifies the phase as being
   * a SurfPhase object.
   *
   * @ingroup thermoprops
   */
  class SurfPhase : public ThermoPhase  {
      
  public:

    //! Constructor.
    /*!
     *  @param n0 Site Density of the Surface Phase
     *            Units: kmol m-2.
     */
    SurfPhase(doublereal n0 = 0.0);

    //! Constructor.
    /*!
     *  @param xmlphase XML node pointing to a SurfPhase description
     */
    SurfPhase(XML_Node& xmlphase);
    

    //! Destructor.
    virtual ~SurfPhase();

    //----- reimplimented methods of class ThermoPhase ------

    //! Equation of state type flag.
    /*!
     *  Redefine this to return cSurf, listed in mix_defs.h.
     */
    virtual int eosType() const { return cSurf; }

    //! Return the Molar Enthalpy. Units: J/kmol.
    /*!
     * For an ideal solution,
     * \f[
     * \hat h(T,P) = \sum_k X_k \hat h^0_k(T),
     * \f]
     * and is a function only of temperature.
     * The standard-state pure-species Enthalpies 
     * \f$ \hat h^0_k(T) \f$ are computed by the species thermodynamic 
     * property manager.
     *
     * \see SpeciesThermo
     */
    virtual doublereal enthalpy_mole() const;

    //! Return the Molar Internal Energy. Units: J/kmol
    /**
     * For a surface phase, the pressure is not a relevant
     * thermodynamic variable, and so the Enthalpy is equal to the
     * Internal Energy.
     */
    virtual doublereal intEnergy_mole() const;

    //! Get the array of chemical potentials at unit activity for the
    //! standard state species at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * These are the standard state chemical potentials \f$ \mu^0_k(T,P)
     * \f$. The values are evaluated at the current
     * temperature and pressure of the solution
     *
     * @param mu0  Output vector of chemical potentials. 
     *             Length: m_kk.
     */
    virtual void getStandardChemPotentials(doublereal* mu0) const;

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

    //! Return a vector of activity concentrations for each species
    /*!
     *  For this phase the activity concentrations,\f$ C^a_k \f$,  are defined to be
     *  equal to the actual concentrations, \f$ C^s_k \f$.
     *  Activity concentrations are 
     *
     *      \f[
     *            C^a_k = C^s_k = \frac{\theta_k  n_0}{s_k}      
     *      \f]
     *
     *  where \f$ \theta_k \f$ is the surface site fraction for species k,
     *  \f$ n_0 \f$ is the surface site density for the phase, and
     *  \f$ s_k \f$ is the surface size of species k.
     *
     * \f$ C^a_k\f$ that are defined such that \f$ a_k = C^a_k /
     * C^0_k, \f$ where \f$ C^0_k \f$ is a standard concentration
     * defined below and \f$ a_k \f$ are activities used in
     * the thermodynamic functions.  These activity concentrations are used
     * by kinetics manager classes to compute the forward and
     * reverse rates of elementary reactions. Note that they may
     * or may not have units of concentration --- they might be
     * partial pressures, mole fractions, or surface coverages,
     *
     * @param c vector of activity concentration (kmol m-2).
     */
    virtual void getActivityConcentrations(doublereal* c) const;

    //! Return the standard concentration for the kth species
    /*!
     * The standard concentration \f$ C^0_k \f$ used to normalize
     * the activity (i.e., generalized) concentration.
     * For this phase, the standard concentration is species-
     * specific
     *
     *        \f[
     *            C^0_k = \frac{n_0}{s_k}      
     *        \f]
     *
     *  This definition implies that the activity is equal to \f$ \theta_k \f$.
     *
     * @param k Optional parameter indicating the species. The default
     *          is to assume this refers to species 0.
     * @return 
     *   Returns the standard Concentration in units of m3 kmol-1.
     */
    virtual doublereal standardConcentration(int k = 0) const;
  
    //! Return the log of the standard concentration for the kth species
    /*!
     * @param k species index (default 0)
     */
    virtual doublereal logStandardConc(int k=0) const;

    //! Set the equation of state parameters from the argument list
    /*!
     * @internal
     * Set equation of state parameters.
     *
     * @param n number of parameters. Must be one
     * @param c array of \a n coefficients
     *           c[0] = The site density (kmol m-2)
     */
    virtual void setParameters(int n, doublereal* c);

    //! Set the Equation-of-State parameters by reading an XML Node Input
    /*!
     *
     * The Equation-of-State data consists of one item, the site density.
     *
     * @param thermoData   Reference to an XML_Node named thermo
     *                     containing the equation-of-state data. The
     *                     XML_Node is within the phase XML_Node describing
     *                     the %SurfPhase object.
     *
     * An example of the contents of the thermoData XML_Node is provided
     * below. The units attribute is used to supply the units of the
     * site density in any convenient form. Internally it is changed
     * into MKS form.
     *
     * @verbatim
     *    <thermo model="Surface">
     *       <site_density units="mol/cm2"> 3e-09 </site_density>
     *    </thermo>
     * @endverbatim
     */
    virtual void setParametersFromXML(const XML_Node& thermoData);

   
    //! Initialize the SurfPhase object after all species have been set up
    /*!
     * @internal Initialize.
     *
     * This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species.  The base class implementation does nothing,
     * and subclasses that do not require initialization do not
     * need to overload this method.  When importing a CTML phase
     * description, this method is called from ThermoPhase::initThermoXML(),
     * which is called from importPhase(),
     * just prior to returning from function importPhase().
     *
     * @see importCTML.cpp
     */
    virtual void initThermo();


    //! Set the initial state of the Surface Phase from an XML_Node
    /*!
     *  State variables that can be set by this routine are
     *  the temperature and the surface site coverages.
     *
     * @param state  XML_Node containing the state information
     *
     * An example of the XML code block is given below.
     *
     * @verbatim
     *   <state>
     *      <temperature units="K">1200.0</temperature>
     *      <coverages>c6H*:0.1, c6HH:0.9</coverages>
     *   </state>
     * @endverbatim
     */
    virtual void setStateFromXML(const XML_Node& state);

    //! Returns the site density
    /*!
     * Site density kmol m-2
     */
    doublereal siteDensity(){ return m_n0; }

    //! Sets the potential energy of species k.
    /*!
     *
     * @param k    Species index
     * @param pe   Value of the potential energy (J kmol-1)
     */
    void setPotentialEnergy(int k, doublereal pe);

    //! Return the potential energy of species k.
    /*!
     * Returns the potential energy of species, k,
     * J kmol-1
     *
     * @param k  Species index
     */
    doublereal potentialEnergy(int k) {return m_pe[k];}

    //! Set the site density of the surface phase (kmol m-2)
    /*!
     *  @param n0 Site density of the surface phase (kmol m-2)
     */
    void setSiteDensity(doublereal n0);

    //! Get the nondimensional Enthalpy functions for the species standard states
    //! at their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param hrt      Output vector of  nondimensional standard state enthalpies.
     *                 Length: m_kk.
     */
    void getEnthalpy_RT(doublereal* hrt) const;

    //! Get the array of nondimensional Entropy functions for the
    //! species standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param sr   Output vector of  nondimensional standard state entropies.
     *             Length: m_kk.
     */
    void getEntropy_R(doublereal* sr) const;

    //! Return the thermodynamic pressure (Pa).
    /*!
     *  This method must be overloaded in derived classes. Since the
     *  mass density, temperature, and mass fractions are stored,
     *  this method should use these values to implement the
     *  mechanical equation of state \f$ P(T, \rho, Y_1, \dots,
     *  Y_K) \f$.
     */
    virtual doublereal pressure() const {
      return m_press;
    }

    //! Set the internally storred pressure (Pa) at constant
    //! temperature and composition
    /*!
     *   This method must be reimplemented in derived classes, where it
     *   may involve the solution of a nonlinear equation. Within %Cantera,
     *   the independent variable is the density. Therefore, this function
     *   solves for the density that will yield the desired input pressure.
     *   The temperature and composition iare held constant during this process.
     *
     *  This base class function will print an error, if not overwritten.
     *
     *  @param p input Pressure (Pa)
     */
    virtual void setPressure(doublereal p) {
      m_press = p;
    }


    //------- new methods defined in this class ----------
    
    //! Set the surface site fractions to a specified state. 
    /*!
     * This routine converts to concentrations
     * in kmol/m2, using m_n0, the surface site density,
     * and size(k), which is defined to be the number of
     * surface sites occupied by the kth molecule.
     * It then calls State::setConcentrations to set the
     * internal concentration in the object.
     *
     * @param theta    This is the surface site fraction
     *                 for the kth species in the surface phase.
     *                 This is a dimensionless quantity.
     * 
     * This routine normalizes the theta's to 1, before application
     */
    void setCoverages(const doublereal* theta);
    
    //! Set the surface site fractions to a specified state. 
    /*!
     * This routine converts to concentrations
     * in kmol/m2, using m_n0, the surface site density,
     * and size(k), which is defined to be the number of
     * surface sites occupied by the kth molecule.
     * It then calls State::setConcentrations to set the
     * internal concentration in the object.
     *
     * @param theta    This is the surface site fraction
     *                 for the kth species in the surface phase.
     *                 This is a dimensionless quantity.
     */
    void setCoveragesNoNorm(const doublereal* theta);

    
    //! Set the coverages from a string of colon-separated  name:value pairs.
    /*!
     *  @param cov  String containing colon-separated  name:value pairs
     */
    void setCoveragesByName(std::string cov);

    //! Return a vector of surface coverages
    /*!
     * Get the coverages.
     *
     * @param theta Array theta must be at least as long as
     *              the number of species.
     */
    void getCoverages(doublereal* theta) const;

  protected:

    //! Surface site density (kmol m-2)
    doublereal m_n0;

    //! log of the surface site density
    doublereal m_logn0;

    //! Minimum temperature for valid species standard state thermo props
    /*!
     * This is the minimum temperature at which all species have valid standard
     * state thermo props defined.
     */
    doublereal m_tmin;

    //! Maximum temperature for valid species standard state thermo props
    /*!
     * This is the maximum temperature at which all species have valid standard
     * state thermo props defined.
     */
    doublereal m_tmax;

    //! Current value of the pressure (Pa)
    doublereal m_press;

    //! Current value of the temperature (Kelvin)
    mutable doublereal    m_tlast;

    //! Temporary storage for the reference state enthalpies
    mutable array_fp      m_h0;

    //! Temporary storage for the reference state entropies
    mutable array_fp      m_s0;

    //! Temporary storage for the reference state heat capacities
    mutable array_fp      m_cp0;

    //! Temporary storage for the reference state gibbs energies
    mutable array_fp      m_mu0;

    //! Temporary work array
    mutable array_fp      m_work;

    //! Potential energy of each species in the surface phase
    /*!
     * @todo Fix potential energy
     * Note, the potential energy terms seem to be orphaned at the moment.
     * They are not connected to the Gibbs free energy calculation in
     * this object
     *
     * @deprecated
     */
    mutable array_fp      m_pe;

    //! vector storring the log of the size of each species.
    /*!
     * The size of each species is defined as the number of surface
     * sites each species occupies.
     */
    mutable array_fp      m_logsize;

  private:

    //! Update the species reference state thermodynamic functions
    /*!
     * The polynomials for the standard state functions are only
     * reevalulated if the temperature has changed.
     *
     * @param force  Boolean, which if true, forces a reevalulation
     *               of the thermo polynomials.
     *               default = false.
     */
    void _updateThermo(bool force=false) const;

  };
}
        
#endif





