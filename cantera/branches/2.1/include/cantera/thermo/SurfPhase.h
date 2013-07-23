/**
 *  @file SurfPhase.h
 *  Header for a simple thermodynamics model of a surface phase
 *  derived from ThermoPhase,
 *  assuming an ideal solution model
 *  (see \ref thermoprops and class \link Cantera::SurfPhase SurfPhase\endlink).
 */

//  Copyright 2002 California Institute of Technology

#ifndef CT_SURFPHASE_H
#define CT_SURFPHASE_H

#include "mix_defs.h"
#include "ThermoPhase.h"

namespace Cantera
{

//!  A simple thermodynamic model for a surface phase,
//!  assuming an ideal solution model.
/*!
 * The surface consists of a grid of equivalent sites.
 * Surface species may be defined to
 * occupy one or more sites. The surface species are assumed to be
 * independent, and thus the species form an ideal solution.
 *
 * The density of surface sites is given by the variable \f$ n_0 \f$,
 * which has SI units of kmol m-2.
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
 * <b> %Application within %Kinetics Managers </b>
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
 * for %Cantera. A new SurfPhase may be created by the following code snippet:
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
 * @code
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
 * @endcode
 *
 * The model attribute, "Surface", on the thermo element identifies the phase as being
 * a SurfPhase object.
 *
 * @ingroup thermoprops
 */
class SurfPhase : public ThermoPhase
{
public:
    //! Constructor.
    /*!
     *  @param n0 Site Density of the Surface Phase
     *            Units: kmol m-2.
     */
    SurfPhase(doublereal n0 = 0.0);

    //! Construct and initialize a SurfPhase ThermoPhase object
    //! directly from an ASCII input file
    /*!
     * @param infile name of the input file
     * @param id     name of the phase id in the file.
     *               If this is blank, the first phase in the file is used.
     */
    SurfPhase(const std::string& infile, std::string id);

    //! Construct and initialize a SurfPhase ThermoPhase object
    //! directly from an XML database
    /*!
     *  @param xmlphase XML node pointing to a SurfPhase description
     */
    SurfPhase(XML_Node& xmlphase);

    //! Copy Constructor
    /*!
     * Copy constructor for the object. Constructed
     * object will be a clone of this object, but will
     * also own all of its data.
     * This is a wrapper around the assignment operator
     *
     * @param right Object to be copied.
     */
    SurfPhase(const SurfPhase& right);

    //! Assignment operator
    /*!
     * Assignment operator for the object. Constructed
     * object will be a clone of this object, but will
     * also own all of its data.
     *
     * @param right Object to be copied.
     */
    SurfPhase& operator=(const SurfPhase& right);

    //! Duplicator from the %ThermoPhase parent class
    /*
     * Given a pointer to a %ThermoPhase object, this function will
     * duplicate the %ThermoPhase object and all underlying structures.
     * This is basically a wrapper around the copy constructor.
     *
     * @return returns a pointer to a %ThermoPhase
     */
    ThermoPhase* duplMyselfAsThermoPhase() const;

    //! Equation of state type flag.
    /*!
     *  Redefine this to return cSurf, listed in mix_defs.h.
     */
    virtual int eosType() const {
        return cSurf;
    }

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
    //! in the mixture. Units (J/kmol)
    /*!
     * @param hbar    Output vector of species partial molar enthalpies.
     *                Length: m_kk. units are J/kmol.
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;

    //! Returns an array of partial molar entropies of the species in the
    //! solution. Units: J/kmol/K.
    /*!
     * @param sbar    Output vector of species partial molar entropies.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;

    //! Return an array of partial molar heat capacities for the
    //! species in the mixture.  Units: J/kmol/K
    /*!
     * @param cpbar   Output vector of species partial molar heat
     *                capacities at constant pressure.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const;

    //! Return an array of partial molar volumes for the
    //! species in the mixture. Units: m^3/kmol.
    /*!
     *  @param vbar   Output vector of species partial molar volumes.
     *                Length = m_kk. units are m^3/kmol.
     */
    virtual void getPartialMolarVolumes(doublereal* vbar) const;

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
    virtual doublereal standardConcentration(size_t k = 0) const;

    //! Return the log of the standard concentration for the kth species
    /*!
     * @param k species index (default 0)
     */
    virtual doublereal logStandardConc(size_t k=0) const;

    //! Set the equation of state parameters from the argument list
    /*!
     * @internal
     * Set equation of state parameters.
     *
     * @param n number of parameters. Must be one
     * @param c array of \a n coefficients
     *           c[0] = The site density (kmol m-2)
     * @deprecated use setSiteDensity()
     */
    virtual void setParameters(int n, doublereal* const c);

    //! Set the Equation-of-State parameters by reading an XML Node Input
    /*!
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
     * @code
     *    <thermo model="Surface">
     *       <site_density units="mol/cm2"> 3e-09 </site_density>
     *    </thermo>
     * @endcode
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
     * @code
     *   <state>
     *      <temperature units="K">1200.0</temperature>
     *      <coverages>c6H*:0.1, c6HH:0.9</coverages>
     *   </state>
     * @endcode
     */
    virtual void setStateFromXML(const XML_Node& state);

    //! Returns the site density
    /*!
     * Site density kmol m-2
     */
    doublereal siteDensity() {
        return m_n0;
    }

    //! Set the site density of the surface phase (kmol m-2)
    /*!
     *  @param n0 Site density of the surface phase (kmol m-2)
     */
    void setSiteDensity(doublereal n0);

    //! Get the nondimensional Gibbs functions for the species
    //! in their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param grt  Output vector of nondimensional standard state gibbs free energies
     *             Length: m_kk.
     */
    virtual void getGibbs_RT(doublereal* grt) const;

    //! Get the nondimensional Enthalpy functions for the species standard states
    //! at their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param hrt      Output vector of  nondimensional standard state enthalpies.
     *                 Length: m_kk.
     */
    virtual void getEnthalpy_RT(doublereal* hrt) const;

    //! Get the array of nondimensional Entropy functions for the
    //! species standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param sr   Output vector of  nondimensional standard state entropies.
     *             Length: m_kk.
     */
    virtual void getEntropy_R(doublereal* sr) const;

    //! Get the nondimensional Heat Capacities at constant
    //! pressure for the species standard states
    //! at the current <I>T</I> and <I>P</I> of the solution
    /*!
     * @param cpr   Output vector of nondimensional standard state heat capacities
     *              Length: m_kk.
     */
    virtual void getCp_R(doublereal* cpr) const;

    //!  Get the molar volumes of the species standard states at the current
    //!  <I>T</I> and <I>P</I> of the solution.
    /*!
     * units = m^3 / kmol
     *
     * @param vol     Output vector containing the standard state volumes.
     *                Length: m_kk.
     */
    virtual void getStandardVolumes(doublereal* vol) const;

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

    //! Set the internally stored pressure (Pa) at constant
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

    //!  Returns the vector of nondimensional
    //!  Gibbs Free Energies of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    /*!
     * @param grt     Output vector containing the nondimensional reference state
     *                Gibbs Free energies.  Length: m_kk.
     */
    virtual void getGibbs_RT_ref(doublereal* grt) const;

    //!  Returns the vector of nondimensional
    //!  enthalpies of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    /*!
     * @param hrt      Output vector of  nondimensional standard state enthalpies.
     *                 Length: m_kk.
     */
    virtual void getEnthalpy_RT_ref(doublereal* hrt) const;

#ifdef H298MODIFY_CAPABILITY

    //! Modify the value of the 298 K Heat of Formation of one species in the phase (J kmol-1)
    /*!
     *   The 298K heat of formation is defined as the enthalpy change to create the standard state
     *   of the species from its constituent elements in their standard states at 298 K and 1 bar.
     *
     *   @param  k           Species k
     *   @param  Hf298New    Specify the new value of the Heat of Formation at 298K and 1 bar
     */
    virtual void modifyOneHf298SS(const size_t& k, const doublereal Hf298New) {
        m_spthermo->modifyOneHf298(k, Hf298New);
        m_tlast += 0.0001234;
    }
#endif

    //!  Returns the vector of nondimensional
    //!  entropies of the reference state at the current temperature
    //!  of the solution and the reference pressure for each species.
    /*!
     * @param er      Output vector containing the nondimensional reference state
     *                entropies.  Length: m_kk.
     */
    virtual void getEntropy_R_ref(doublereal* er) const;

    //! Returns the vector of nondimensional constant pressure heat capacities
    //! of the reference state at the current temperature of the solution and
    //! reference pressure for each species.
    /*!
     * @param cprt   Output vector of nondimensional reference state
     *               heat capacities at constant pressure for the species.
     *               Length: m_kk
     */
    virtual void getCp_R_ref(doublereal* cprt) const;

    //! Set the surface site fractions to a specified state.
    /*!
     * This routine converts to concentrations
     * in kmol/m2, using m_n0, the surface site density,
     * and size(k), which is defined to be the number of
     * surface sites occupied by the kth molecule.
     * It then calls Phase::setConcentrations to set the
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
     * It then calls Phase::setConcentrations to set the
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
    void setCoveragesByName(const std::string& cov);

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

    //! Current value of the pressure (Pa)
    doublereal m_press;

    //! Current value of the temperature (Kelvin)
    mutable doublereal    m_tlast;

    //! Temporary storage for the reference state enthalpies
    mutable vector_fp      m_h0;

    //! Temporary storage for the reference state entropies
    mutable vector_fp      m_s0;

    //! Temporary storage for the reference state heat capacities
    mutable vector_fp      m_cp0;

    //! Temporary storage for the reference state gibbs energies
    mutable vector_fp      m_mu0;

    //! Temporary work array
    mutable vector_fp      m_work;

    //! vector storing the log of the size of each species.
    /*!
     * The size of each species is defined as the number of surface
     * sites each species occupies.
     */
    mutable vector_fp      m_logsize;

private:
    //! Update the species reference state thermodynamic functions
    /*!
     * The polynomials for the standard state functions are only
     * reevaluated if the temperature has changed.
     *
     * @param force  Boolean, which if true, forces a reevaluation
     *               of the thermo polynomials.
     *               default = false.
     */
    void _updateThermo(bool force=false) const;
};
}

#endif
