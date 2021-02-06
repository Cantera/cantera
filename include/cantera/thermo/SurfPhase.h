/**
 *  @file SurfPhase.h
 *  Header for a simple thermodynamics model of a surface phase
 *  derived from ThermoPhase,
 *  assuming an ideal solution model
 *  (see \ref thermoprops and class \link Cantera::SurfPhase SurfPhase\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_SURFPHASE_H
#define CT_SURFPHASE_H

#include "ThermoPhase.h"

namespace Cantera
{

//! A simple thermodynamic model for a surface phase, assuming an ideal solution
//! model.
/*!
 * The surface consists of a grid of equivalent sites. Surface species may be
 * defined to occupy one or more sites. The surface species are assumed to be
 * independent, and thus the species form an ideal solution.
 *
 * The density of surface sites is given by the variable \f$ n_0 \f$,
 * which has SI units of kmol m-2.
 *
 * ## Specification of Species Standard State Properties
 *
 * It is assumed that the reference state thermodynamics may be obtained by a
 * pointer to a populated species thermodynamic property manager class (see
 * ThermoPhase::m_spthermo). How to relate pressure changes to the reference
 * state thermodynamics is resolved at this level.
 *
 * Pressure is defined as an independent variable in this phase. However, it has
 * no effect on any quantities, as the molar concentration is a constant.
 *
 * Therefore, The standard state internal energy for species *k* is equal to the
 * enthalpy for species *k*.
 *
 *       \f[
 *            u^o_k = h^o_k
 *       \f]
 *
 * Also, the standard state chemical potentials, entropy, and heat capacities
 * are independent of pressure. The standard state Gibbs free energy is obtained
 * from the enthalpy and entropy functions.
 *
 * ## Specification of Solution Thermodynamic Properties
 *
 * The activity of species defined in the phase is given by
 *       \f[
 *            a_k = \theta_k
 *       \f]
 *
 * The chemical potential for species *k* is equal to
 *       \f[
 *            \mu_k(T,P) = \mu^o_k(T) + R T \log(\theta_k)
 *       \f]
 *
 * Pressure is defined as an independent variable in this phase. However, it has
 * no effect on any quantities, as the molar concentration is a constant.
 *
 * The internal energy for species k is equal to the enthalpy for species *k*
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
 * ## %Application within Kinetics Managers
 *
 * The activity concentration,\f$  C^a_k \f$, used by the kinetics manager, is equal to
 * the actual concentration, \f$ C^s_k \f$, and is given by the following
 * expression.
 *       \f[
 *            C^a_k = C^s_k = \frac{\theta_k  n_0}{s_k}
 *       \f]
 *
 * The standard concentration for species *k* is:
 *        \f[
 *            C^0_k = \frac{n_0}{s_k}
 *        \f]
 *
 * ## Instantiation of the Class
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
 * ## XML Example
 *
 * An example of an XML Element named phase setting up a SurfPhase object named
 * diamond_100 is given below.
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
    SurfPhase(doublereal n0 = 1.0);

    //! Construct and initialize a SurfPhase ThermoPhase object directly from an
    //! ASCII input file
    /*!
     * @param infile name of the input file
     * @param id     name of the phase id in the file.
     *               If this is blank, the first phase in the file is used.
     */
    SurfPhase(const std::string& infile, const std::string& id);

    //! Construct and initialize a SurfPhase ThermoPhase object directly from an
    //! XML database
    /*!
     *  @param xmlphase XML node pointing to a SurfPhase description
     *
     * @deprecated The XML input format is deprecated and will be removed in
     *     Cantera 3.0.
     */
    SurfPhase(XML_Node& xmlphase);

    virtual std::string type() const {
        return "Surf";
    }

    //! Return the Molar Enthalpy. Units: J/kmol.
    /*!
     * For an ideal solution,
     * \f[
     * \hat h(T,P) = \sum_k X_k \hat h^0_k(T),
     * \f]
     * and is a function only of temperature. The standard-state pure-species
     * Enthalpies \f$ \hat h^0_k(T) \f$ are computed by the species
     * thermodynamic property manager.
     *
     * \see MultiSpeciesThermo
     */
    virtual doublereal enthalpy_mole() const;

    //! Return the Molar Internal Energy. Units: J/kmol
    /**
     * For a surface phase, the pressure is not a relevant thermodynamic
     * variable, and so the Enthalpy is equal to the Internal Energy.
     */
    virtual doublereal intEnergy_mole() const;

    //! Return the Molar Entropy. Units: J/kmol-K
    /**
     * \f[
     *  \hat s(T,P) = \sum_k X_k (\hat s^0_k(T) - R \log(\theta_k))
     * \f]
     */
    virtual doublereal entropy_mole() const;

    virtual doublereal cp_mole() const;
    virtual doublereal cv_mole() const;

    virtual void getChemPotentials(doublereal* mu) const;
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;
    virtual void getPartialMolarEntropies(doublereal* sbar) const;
    virtual void getPartialMolarCp(doublereal* cpbar) const;
    virtual void getPartialMolarVolumes(doublereal* vbar) const;
    virtual void getStandardChemPotentials(doublereal* mu0) const;

    //! Return a vector of activity concentrations for each species
    /*!
     * For this phase the activity concentrations,\f$ C^a_k \f$, are defined to
     * be equal to the actual concentrations, \f$ C^s_k \f$. Activity
     * concentrations are
     *
     *      \f[
     *            C^a_k = C^s_k = \frac{\theta_k  n_0}{s_k}
     *      \f]
     *
     * where \f$ \theta_k \f$ is the surface site fraction for species k,
     * \f$ n_0 \f$ is the surface site density for the phase, and
     * \f$ s_k \f$ is the surface size of species k.
     *
     * \f$ C^a_k\f$ that are defined such that \f$ a_k = C^a_k / C^0_k, \f$
     * where \f$ C^0_k \f$ is a standard concentration defined below and \f$ a_k
     * \f$ are activities used in the thermodynamic functions.  These activity
     * concentrations are used by kinetics manager classes to compute the
     * forward and reverse rates of elementary reactions. Note that they may or
     * may not have units of concentration --- they might be partial pressures,
     * mole fractions, or surface coverages,
     *
     * @param c vector of activity concentration (kmol m-2).
     */
    virtual void getActivityConcentrations(doublereal* c) const;

    //! Return the standard concentration for the kth species
    /*!
     * The standard concentration \f$ C^0_k \f$ used to normalize the activity
     * (i.e., generalized) concentration. For this phase, the standard
     * concentration is species- specific
     *
     *        \f[
     *            C^0_k = \frac{n_0}{s_k}
     *        \f]
     *
     * This definition implies that the activity is equal to \f$ \theta_k \f$.
     *
     * @param k Optional parameter indicating the species. The default
     *          is to assume this refers to species 0.
     * @return
     *   Returns the standard Concentration in units of m3 kmol-1.
     */
    virtual doublereal standardConcentration(size_t k = 0) const;
    virtual doublereal logStandardConc(size_t k=0) const;

    //! Set the equation of state parameters from the argument list
    /*!
     * @internal
     * Set equation of state parameters.
     *
     * @param n number of parameters. Must be one
     * @param c array of \a n coefficients
     *           c[0] = The site density (kmol m-2)
     */
    virtual void setParameters(int n, doublereal* const c);

    //! Set the Equation-of-State parameters by reading an XML Node Input
    /*!
     * The Equation-of-State data consists of one item, the site density.
     *
     * @param thermoData   Reference to an XML_Node named thermo containing the
     *                     equation-of-state data. The XML_Node is within the
     *                     phase XML_Node describing the SurfPhase object.
     *
     * An example of the contents of the thermoData XML_Node is provided below.
     * The units attribute is used to supply the units of the site density in
     * any convenient form. Internally it is changed into MKS form.
     *
     * @code
     *    <thermo model="Surface">
     *       <site_density units="mol/cm2"> 3e-09 </site_density>
     *    </thermo>
     * @endcode
     *
     * @deprecated The XML input format is deprecated and will be removed in
     *     Cantera 3.0.
     */
    virtual void setParametersFromXML(const XML_Node& thermoData);
    virtual void initThermo();

    virtual bool addSpecies(shared_ptr<Species> spec);

    //! Set the initial state of the Surface Phase from an XML_Node
    /*!
     * State variables that can be set by this routine are the temperature and
     * the surface site coverages.
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
     *
     * @deprecated The XML input format is deprecated and will be removed in
     *     Cantera 3.0.
     */
    virtual void setStateFromXML(const XML_Node& state);

    //! Returns the site density
    /*!
     * Site density kmol m-2
     */
    doublereal siteDensity() {
        return m_n0;
    }

    //! Returns the number of sites occupied by one molecule of species *k*.
    virtual double size(size_t k) const {
        return m_speciesSize[k];
    }

    //! Set the site density of the surface phase (kmol m-2)
    /*!
     *  @param n0 Site density of the surface phase (kmol m-2)
     */
    void setSiteDensity(doublereal n0);

    virtual void getGibbs_RT(doublereal* grt) const;
    virtual void getEnthalpy_RT(doublereal* hrt) const;
    virtual void getEntropy_R(doublereal* sr) const;
    virtual void getCp_R(doublereal* cpr) const;
    virtual void getStandardVolumes(doublereal* vol) const;

    //! Return the thermodynamic pressure (Pa).
    virtual doublereal pressure() const {
        return m_press;
    }

    //! Set the internally stored pressure (Pa) at constant temperature and
    //! composition
    /*!
     *  @param p input Pressure (Pa)
     */
    virtual void setPressure(doublereal p) {
        m_press = p;
    }

    virtual void getPureGibbs(doublereal* g) const;
    virtual void getGibbs_RT_ref(doublereal* grt) const;
    virtual void getEnthalpy_RT_ref(doublereal* hrt) const;
    virtual void getEntropy_R_ref(doublereal* er) const;
    virtual void getCp_R_ref(doublereal* cprt) const;

    //! Set the surface site fractions to a specified state.
    /*!
     * This routine converts to concentrations in kmol/m2, using m_n0, the
     * surface site density, and size(k), which is defined to be the number of
     * surface sites occupied by the kth molecule. It then calls
     * Phase::setConcentrations to set the internal concentration in the object.
     *
     * @param theta    This is the surface site fraction for the kth species in
     *                 the surface phase. This is a dimensionless quantity.
     *
     * This routine normalizes the theta's to 1, before application
     */
    void setCoverages(const doublereal* theta);

    //! Set the surface site fractions to a specified state.
    /*!
     * This routine converts to concentrations in kmol/m2, using m_n0, the
     * surface site density, and size(k), which is defined to be the number of
     * surface sites occupied by the kth molecule. It then calls
     * Phase::setConcentrations to set the internal concentration in the object.
     *
     * @param theta    This is the surface site fraction for the kth species in
     *                 the surface phase. This is a dimensionless quantity.
     */
    void setCoveragesNoNorm(const doublereal* theta);

    //! Set the coverages from a string of colon-separated name:value pairs.
    /*!
     *  @param cov  String containing colon-separated name:value pairs
     */
    void setCoveragesByName(const std::string& cov);

    //! Set the coverages from a map of name:value pairs
    void setCoveragesByName(const compositionMap& cov);

    //! Return a vector of surface coverages
    /*!
     * Get the coverages.
     *
     * @param theta Array theta must be at least as long as the number of
     *              species.
     */
    void getCoverages(doublereal* theta) const;

    //! @copydoc ThermoPhase::setState
    /*!
     * Additionally uses the key `coverages` to set the fractional coverages.
     */
    virtual void setState(const AnyMap& state);

protected:
    //! Surface site density (kmol m-2)
    doublereal m_n0;

    //! Vector of species sizes (number of sites occupied). length m_kk.
    vector_fp m_speciesSize;

    //! log of the surface site density
    doublereal m_logn0;

    //! Current value of the pressure (Pa)
    doublereal m_press;

    //! Temporary storage for the reference state enthalpies
    mutable vector_fp m_h0;

    //! Temporary storage for the reference state entropies
    mutable vector_fp m_s0;

    //! Temporary storage for the reference state heat capacities
    mutable vector_fp m_cp0;

    //! Temporary storage for the reference state Gibbs energies
    mutable vector_fp m_mu0;

    //! Temporary work array
    mutable vector_fp m_work;

    //! vector storing the log of the size of each species.
    /*!
     * The size of each species is defined as the number of surface sites each
     * species occupies.
     */
    mutable vector_fp m_logsize;

private:
    //! Update the species reference state thermodynamic functions
    /*!
     * The polynomials for the standard state functions are only reevaluated if
     * the temperature has changed.
     *
     * @param force  Boolean, which if true, forces a reevaluation of the thermo
     *               polynomials. default = false.
     */
    void _updateThermo(bool force=false) const;
};

}

#endif
