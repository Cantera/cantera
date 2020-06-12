/**
 *  @file LatticePhase.h Header for a simple thermodynamics model of a bulk
 *      phase derived from ThermoPhase, assuming a lattice of solid atoms (see
 *      \ref thermoprops and class \link Cantera::LatticePhase
 *      LatticePhase\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_LATTICE_H
#define CT_LATTICE_H

#include "ThermoPhase.h"

namespace Cantera
{

//! A simple thermodynamic model for a bulk phase, assuming a lattice of solid
//! atoms
/*!
 * The bulk consists of a matrix of equivalent sites whose molar density does
 * not vary with temperature or pressure. The thermodynamics obeys the ideal
 * solution laws. The phase and the pure species phases which comprise the
 * standard states of the species are assumed to have zero volume expansivity
 * and zero isothermal compressibility.
 *
 * The density of matrix sites is given by the variable \f$ C_o \f$, which has
 * SI units of kmol m-3.
 *
 * ## Specification of Species Standard State Properties
 *
 * It is assumed that the reference state thermodynamics may be obtained by a
 * pointer to a populated species thermodynamic property manager class (see
 * ThermoPhase::m_spthermo). However, how to relate pressure changes to the
 * reference state thermodynamics is within this class.
 *
 * Pressure is defined as an independent variable in this phase. However, it has
 * no effect on any quantities, as the molar concentration is a constant.
 *
 * The standard state enthalpy function is given by the following relation,
 * which has a weak dependence on the system pressure, \f$P\f$.
 *
 * \f[
 *     h^o_k(T,P) =
 *            h^{ref}_k(T) +  \left( \frac{P - P_{ref}}{C_o} \right)
 * \f]
 *
 * For an incompressible substance, the molar internal energy is independent of
 * pressure. Since the thermodynamic properties are specified by giving the
 * standard-state enthalpy, the term \f$ \frac{P_{ref}}{C_o} \f$ is subtracted
 * from the specified reference molar enthalpy to compute the standard state
 * molar internal energy:
 *
 * \f[
 *      u^o_k(T,P) = h^{ref}_k(T) - \frac{P_{ref}}{C_o}
 * \f]
 *
 * The standard state heat capacity, internal energy, and entropy are
 * independent of pressure. The standard state Gibbs free energy is obtained
 * from the enthalpy and entropy functions.
 *
 * The standard state molar volume is independent of temperature, pressure, and
 * species identity:
 *
 * \f[
 *      V^o_k(T,P) = \frac{1.0}{C_o}
 * \f]
 *
 * ## Specification of Solution Thermodynamic Properties
 *
 * The activity of species \f$ k \f$ defined in the phase, \f$ a_k \f$, is given
 * by the ideal solution law:
 *
 * \f[
 *      a_k = X_k ,
 * \f]
 *
 * where \f$ X_k \f$ is the mole fraction of species *k*. The chemical potential
 * for species *k* is equal to
 *
 * \f[
 *      \mu_k(T,P) = \mu^o_k(T, P) + R T \log(X_k)
 * \f]
 *
 * The partial molar entropy for species *k* is given by the following relation,
 *
 * \f[
 *      \tilde{s}_k(T,P) = s^o_k(T,P) - R \log(X_k) = s^{ref}_k(T) - R \log(X_k)
 * \f]
 *
 * The partial molar enthalpy for species *k* is
 *
 * \f[
 *      \tilde{h}_k(T,P) = h^o_k(T,P) = h^{ref}_k(T) + \left( \frac{P - P_{ref}}{C_o} \right)
 * \f]
 *
 * The partial molar Internal Energy for species *k* is
 *
 * \f[
 *      \tilde{u}_k(T,P) = u^o_k(T,P) = u^{ref}_k(T)
 * \f]
 *
 * The partial molar Heat Capacity for species *k* is
 *
 * \f[
 *      \tilde{Cp}_k(T,P) = Cp^o_k(T,P) = Cp^{ref}_k(T)
 * \f]
 *
 * The partial molar volume is independent of temperature, pressure, and species
 * identity:
 *
 * \f[
 *      \tilde{V}_k(T,P) =  V^o_k(T,P) = \frac{1.0}{C_o}
 * \f]
 *
 * It is assumed that the reference state thermodynamics may be obtained by a
 * pointer to a populated species thermodynamic property manager class (see
 * ThermoPhase::m_spthermo). How to relate pressure changes to the reference
 * state thermodynamics is resolved at this level.
 *
 * Pressure is defined as an independent variable in this phase. However, it
 * only has a weak dependence on the enthalpy, and doesn't effect the molar
 * concentration.
 *
 * ## %Application within Kinetics Managers
 *
 * \f$ C^a_k\f$ are defined such that \f$ C^a_k = a_k = X_k \f$. \f$ C^s_k \f$,
 * the standard concentration, is defined to be equal to one. \f$ a_k \f$ are
 * activities used in the thermodynamic functions.  These activity (or
 * generalized) concentrations are used by kinetics manager classes to compute
 * the forward and reverse rates of elementary reactions. The activity
 * concentration,\f$  C^a_k \f$, is given by the following expression.
 *
 * \f[
 *      C^a_k = C^s_k  X_k  =  X_k
 * \f]
 *
 * The standard concentration for species *k* is identically one
 *
 * \f[
 *     C^s_k =  C^s = 1.0
 * \f]
 *
 * For example, a bulk-phase binary gas reaction between species j and k,
 * producing a new species l would have the following equation for its rate of
 * progress variable, \f$ R^1 \f$, which has units of kmol m-3 s-1.
 *
 * \f[
 *    R^1 = k^1 C_j^a C_k^a =  k^1  X_j X_k
 * \f]
 *
 * The reverse rate constant can then be obtained from the law of microscopic
 * reversibility and the equilibrium expression for the system.
 *
 * \f[
 *     \frac{X_j X_k}{ X_l} = K_a^{o,1} = \exp(\frac{\mu^o_l - \mu^o_j - \mu^o_k}{R T} )
 * \f]
 *
 * \f$  K_a^{o,1} \f$ is the dimensionless form of the equilibrium constant,
 * associated with the pressure dependent standard states \f$ \mu^o_l(T,P) \f$
 * and their associated activities,
 * \f$ a_l \f$, repeated here:
 *
 * \f[
 *      \mu_l(T,P) = \mu^o_l(T, P) + R T \log(a_l)
 * \f]
 *
 * The concentration equilibrium constant, \f$ K_c \f$, may be obtained by
 * changing over to activity concentrations. When this is done:
 *
 * \f[
 *     \frac{C^a_j C^a_k}{ C^a_l} = C^o K_a^{o,1} = K_c^1 =
 *         \exp(\frac{\mu^{o}_l - \mu^{o}_j - \mu^{o}_k}{R T} )
 * \f]
 *
 * %Kinetics managers will calculate the concentration equilibrium constant, \f$
 * K_c \f$, using the second and third part of the above expression as a
 * definition for the concentration equilibrium constant.
 *
 * ## Instantiation of the Class
 *
 * The constructor for this phase is located in the default ThermoFactory for
 * %Cantera. A new LatticePhase object may be created by the following code
 * snippet:
 *
 * @code
 *    XML_Node *xc = get_XML_File("O_lattice_SiO2.xml");
 *    XML_Node * const xs = xc->findNameID("phase", "O_lattice_SiO2");
 *    ThermoPhase *tp = newPhase(*xs);
 *    LatticePhase *o_lattice = dynamic_cast <LatticPhase *>(tp);
 * @endcode
 *
 * or by the following constructor:
 *
 * @code
 *    XML_Node *xc = get_XML_File("O_lattice_SiO2.xml");
 *    XML_Node * const xs = xc->findNameID("phase", "O_lattice_SiO2");
 *    LatticePhase *o_lattice = new LatticePhase(*xs);
 * @endcode
 *
 * The XML file used in this example is listed in the next section
 *
 * ## XML Example
 *
 * An example of an XML Element named phase setting up a LatticePhase object
 * named "O_lattice_SiO2" is given below.
 *
 * @code
 * <!--     phase O_lattice_SiO2      -->
 *   <phase dim="3" id="O_lattice_SiO2">
 *     <elementArray datasrc="elements.xml"> Si  H  He </elementArray>
 *     <speciesArray datasrc="#species_data">
 *       O_O  Vac_O
 *     </speciesArray>
 *     <reactionArray datasrc="#reaction_data"/>
 *     <thermo model="Lattice">
 *       <site_density> 73.159 </site_density>
 *     </thermo>
 *     <kinetics model="BulkKinetics"/>
 *     <transport model="None"/>
 *  </phase>
 *  @endcode
 *
 * The model attribute "Lattice" of the thermo XML element identifies the phase
 * as being of the type handled by the LatticePhase object.
 *
 * @ingroup thermoprops
 */
class LatticePhase : public ThermoPhase
{
public:
    //! Base Empty constructor
    LatticePhase();

    //! Full constructor for a lattice phase
    /*!
     * @param inputFile String name of the input file
     * @param id        string id of the phase name
     */
    LatticePhase(const std::string& inputFile, const std::string& id = "");

    //! Full constructor for a water phase
    /*!
     * @param phaseRef  XML node referencing the lattice phase.
     * @param id        string id of the phase name
     *
     * @deprecated The XML input format is deprecated and will be removed in
     *     Cantera 3.0.
     */
    LatticePhase(XML_Node& phaseRef, const std::string& id = "");

    virtual std::string type() const {
        return "Lattice";
    }

    virtual bool isCompressible() const {
        return false;
    }

    std::map<std::string, size_t> nativeState() const {
        return { {"T", 0}, {"P", 1}, {"X", 2} };
    }

    //! @name Molar Thermodynamic Properties of the Solution
    //! @{

    //! Return the Molar Enthalpy. Units: J/kmol.
    /*!
     * For an ideal solution,
     *
     *   \f[
     *    \hat h(T,P) = \sum_k X_k \hat h^0_k(T,P),
     *   \f]
     *
     * The standard-state pure-species Enthalpies \f$ \hat h^0_k(T,P) \f$ are
     * computed first by the species reference state thermodynamic property
     * manager and then a small pressure dependent term is added in.
     *
     * \see MultiSpeciesThermo
     */
    virtual doublereal enthalpy_mole() const;

    //! Molar entropy of the solution. Units: J/kmol/K
    /*!
     * For an ideal, constant partial molar volume solution mixture with
     * pure species phases which exhibit zero volume expansivity:
     *   \f[
     *       \hat s(T, P, X_k) = \sum_k X_k \hat s^0_k(T)  - \hat R  \sum_k X_k log(X_k)
     *   \f]
     * The reference-state pure-species entropies \f$ \hat s^0_k(T,p_{ref}) \f$
     * are computed by the species thermodynamic property manager. The pure
     * species entropies are independent of pressure since the volume
     * expansivities are equal to zero.
     *
     * Units: J/kmol/K.
     *
     * @see MultiSpeciesThermo
     */
    virtual doublereal entropy_mole() const;

    //! Molar heat capacity at constant pressure of the solution.
    //! Units: J/kmol/K.
    /*!
     * For an ideal, constant partial molar volume solution mixture with
     * pure species phases which exhibit zero volume expansivity:
     *   \f[
     *    \hat c_p(T,P) = \sum_k X_k \hat c^0_{p,k}(T) .
     *   \f]
     * The heat capacity is independent of pressure. The reference-state pure-
     * species heat capacities \f$ \hat c^0_{p,k}(T) \f$ are computed by the
     * species thermodynamic property manager.
     *
     * @see MultiSpeciesThermo
     */
    virtual doublereal cp_mole() const;

    //! Molar heat capacity at constant volume of the solution.
    //! Units: J/kmol/K.
    /*!
     * For an ideal, constant partial molar volume solution mixture with
     * pure species phases which exhibit zero volume expansivity:
     *  \f[
     *     \hat c_v(T,P) = \hat c_p(T,P)
     *  \f]
     *
     * The two heat capacities are equal.
     */
    virtual doublereal cv_mole() const;

    //@}
    /// @name Mechanical Equation of State Properties
    /**
     * In this equation of state implementation, the density is a function only
     * of the mole fractions. Therefore, it can't be an independent variable.
     * Instead, the pressure is used as the independent variable. Functions
     * which try to set the thermodynamic state by calling setDensity() may
     * cause an exception to be thrown.
     */
    //@{

    //! Pressure. Units: Pa.
    /*!
     * For this incompressible system, we return the internally stored
     * independent value of the pressure.
     */
    virtual doublereal pressure() const {
        return m_Pcurrent;
    }

    //! Set the internally stored pressure (Pa) at constant temperature and
    //! composition
    /*!
     * This method sets the pressure within the object. The mass density is not
     * a function of pressure.
     *
     * @param p   Input Pressure (Pa)
     */
    virtual void setPressure(doublereal p);

    //! Calculate the density of the mixture using the partial molar volumes and
    //! mole fractions as input
    /*!
     * The formula for this is
     *
     * \f[
     *      \rho = \frac{\sum_k{X_k W_k}}{\sum_k{X_k V_k}}
     * \f]
     *
     * where \f$X_k\f$ are the mole fractions, \f$W_k\f$ are the molecular
     * weights, and \f$V_k\f$ are the pure species molar volumes.
     *
     * Note, the basis behind this formula is that in an ideal solution the
     * partial molar volumes are equal to the pure species molar volumes. We
     * have additionally specified in this class that the pure species molar
     * volumes are independent of temperature and pressure.
     */
    doublereal calcDensity();

    //@}
    /// @name Activities, Standard States, and Activity Concentrations
    /**
     * The activity \f$a_k\f$ of a species in solution is related to the
     * chemical potential by \f[ \mu_k = \mu_k^0(T) + \hat R T \log a_k. \f] The
     * quantity \f$\mu_k^0(T,P)\f$ is the chemical potential at unit activity,
     * which depends only on temperature and the pressure. Activity is assumed
     * to be molality-based here.
     */
    //@{

    virtual Units standardConcentrationUnits() const;
    virtual void getActivityConcentrations(doublereal* c) const;

    //! Return the standard concentration for the kth species
    /*!
     * The standard concentration \f$ C^0_k \f$ used to normalize
     * the activity (i.e., generalized) concentration for use
     *
     * For the time being, we will use the concentration of pure solvent for the
     * the standard concentration of all species. This has the effect of making
     * mass-action reaction rates based on the molality of species proportional
     * to the molality of the species.
     *
     * @param k Optional parameter indicating the species. The default is to
     *         assume this refers to species 0.
     * @returns the standard Concentration in units of m^3/kmol.
     *
     * @param k Species index
     */
    virtual doublereal standardConcentration(size_t k=0) const;
    virtual doublereal logStandardConc(size_t k=0) const;

    //! Get the array of non-dimensional activity coefficients at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     *  For this phase, the activity coefficients are all equal to one.
     *
     * @param ac Output vector of activity coefficients. Length: m_kk.
     */
    virtual void getActivityCoefficients(doublereal* ac) const;

    //@}
    /// @name  Partial Molar Properties of the Solution
    //@{

    //! Get the species chemical potentials. Units: J/kmol.
    /*!
     * This function returns a vector of chemical potentials of the species in
     * solid solution at the current temperature, pressure and mole fraction of
     * the solid solution.
     *
     * @param mu  Output vector of species chemical
     *            potentials. Length: m_kk. Units: J/kmol
     */
    virtual void getChemPotentials(doublereal* mu) const;

    //@}
    /// @name  Partial Molar Properties of the Solution
    //@{

    /**
     * Returns an array of partial molar enthalpies for the species in the
     * mixture. Units (J/kmol). For this phase, the partial molar enthalpies are
     * equal to the pure species enthalpies
     * \f[
     * \bar h_k(T,P) = \hat h^{ref}_k(T) + (P - P_{ref}) \hat V^0_k
     * \f]
     * The reference-state pure-species enthalpies, \f$ \hat h^{ref}_k(T) \f$,
     * at the reference pressure,\f$ P_{ref} \f$, are computed by the species
     * thermodynamic property manager. They are polynomial functions of
     * temperature.
     * @see MultiSpeciesThermo
     *
     * @param hbar Output vector containing partial molar enthalpies.
     *             Length: m_kk.
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;

    /**
     * Returns an array of partial molar entropies of the species in the
     * solution. Units: J/kmol/K. For this phase, the partial molar entropies
     * are equal to the pure species entropies plus the ideal solution
     * contribution.
     * \f[
     * \bar s_k(T,P) =  \hat s^0_k(T) - R log(X_k)
     * \f]
     * The reference-state pure-species entropies,\f$ \hat s^{ref}_k(T) \f$, at
     * the reference pressure, \f$ P_{ref} \f$, are computed by the species
     * thermodynamic property manager. They are polynomial functions of
     * temperature.
     * @see MultiSpeciesThermo
     *
     * @param sbar Output vector containing partial molar entropies.
     *             Length: m_kk.
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;

    /**
     * Returns an array of partial molar Heat Capacities at constant pressure of
     * the species in the solution. Units: J/kmol/K. For this phase, the partial
     * molar heat capacities are equal to the standard state heat capacities.
     *
     * @param cpbar  Output vector of partial heat capacities. Length: m_kk.
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const;

    virtual void getPartialMolarVolumes(doublereal* vbar) const;
    virtual void getStandardChemPotentials(doublereal* mu) const;
    virtual void getPureGibbs(doublereal* gpure) const;

    //@}
    /// @name  Properties of the Standard State of the Species in the Solution
    //@{

    //! Get the nondimensional Enthalpy functions for the species standard
    //! states at their standard states at the current *T* and *P* of the
    //! solution.
    /*!
     * A small pressure dependent term is added onto the reference state enthalpy
     * to get the pressure dependence of this term.
     *
     * \f[
     *    h^o_k(T,P) = h^{ref}_k(T) +  \left( \frac{P - P_{ref}}{C_o} \right)
     * \f]
     *
     * The reference state thermodynamics is obtained by a pointer to a
     * populated species thermodynamic property manager class (see
     * ThermoPhase::m_spthermo). How to relate pressure changes to the reference
     * state thermodynamics is resolved at this level.
     *
     * @param hrt      Output vector of nondimensional standard state enthalpies.
     *                 Length: m_kk.
     */
    virtual void getEnthalpy_RT(doublereal* hrt) const;

    //! Get the array of nondimensional Entropy functions for the species
    //! standard states at the current *T* and *P* of the solution.
    /*!
     * The entropy of the standard state is defined as independent of
     * pressure here.
     *
     * \f[
     *      s^o_k(T,P) = s^{ref}_k(T)
     * \f]
     *
     * The reference state thermodynamics is obtained by a pointer to a
     * populated species thermodynamic property manager class (see
     * ThermoPhase::m_spthermo). How to relate pressure changes to the reference
     * state thermodynamics is resolved at this level.
     *
     * @param sr   Output vector of nondimensional standard state entropies.
     *             Length: m_kk.
     */
    virtual void getEntropy_R(doublereal* sr) const;

    //! Get the nondimensional Gibbs functions for the species standard states
    //! at the current *T* and *P* of the solution.
    /*!
     * The standard Gibbs free energies are obtained from the enthalpy and
     * entropy formulation.
     *
     * \f[
     *      g^o_k(T,P) = h^{o}_k(T,P) - T s^{o}_k(T,P)
     * \f]
     *
     * @param grt  Output vector of nondimensional standard state Gibbs free
     *             energies. Length: m_kk.
     */
    virtual void getGibbs_RT(doublereal* grt) const;

    //! Get the nondimensional Heat Capacities at constant pressure for the
    //! species standard states at the current *T* and *P* of the solution
    /*!
     * The heat capacity of the standard state is independent of pressure
     *
     * \f[
     *      Cp^o_k(T,P) = Cp^{ref}_k(T)
     * \f]
     *
     * The reference state thermodynamics is obtained by a pointer to a
     * populated species thermodynamic property manager class (see
     * ThermoPhase::m_spthermo). How to relate pressure changes to the reference
     * state thermodynamics is resolved at this level.
     *
     * @param cpr   Output vector of nondimensional standard state heat
     *              capacities. Length: m_kk.
     */
    virtual void getCp_R(doublereal* cpr) const;

    //! Get the molar volumes of the species standard states at the current
    //! *T* and *P* of the solution.
    /*!
     * units = m^3 / kmol
     *
     * @param vol     Output vector containing the standard state volumes.
     *                Length: m_kk.
     */
    virtual void getStandardVolumes(doublereal* vol) const;

    //@}
    /// @name Thermodynamic Values for the Species Reference States
    //@{

    const vector_fp& enthalpy_RT_ref() const;

    //! Returns a reference to the dimensionless reference state Gibbs free
    //! energy vector.
    /*!
     * This function is part of the layer that checks/recalculates the reference
     * state thermo functions.
     */
    const vector_fp& gibbs_RT_ref() const;

    virtual void getGibbs_RT_ref(doublereal* grt) const;
    virtual void getGibbs_ref(doublereal* g) const;

    //! Returns a reference to the dimensionless reference state Entropy vector.
    /*!
     * This function is part of the layer that checks/recalculates the reference
     * state thermo functions.
     */
    const vector_fp& entropy_R_ref() const;

    //! Returns a reference to the dimensionless reference state Heat Capacity
    //! vector.
    /*!
     * This function is part of the layer that checks/recalculates the reference
     * state thermo functions.
     */
    const vector_fp& cp_R_ref() const;

    //@}
    /// @name  Utilities for Initialization of the Object
    //@{

    virtual bool addSpecies(shared_ptr<Species> spec);

    //! Set the density of lattice sites [kmol/m^3]
    void setSiteDensity(double sitedens);

    virtual void initThermo();

    //! Set equation of state parameter values from XML entries.
    /*!
     * This method is called by function importPhase() when processing a phase
     * definition in an input file. It should be overloaded in subclasses to set
     * any parameters that are specific to that particular phase
     * model. Note, this method is called before the phase is
     * initialized with elements and/or species.
     *
     * For this phase, the molar density of the phase is specified in this
     * block, and is a required parameter.
     *
     * @param eosdata An XML_Node object corresponding to
     *                the "thermo" entry for this phase in the input file.
     *
     * eosdata points to the thermo block, and looks like this:
     *
     * @code
     * <phase id="O_lattice_SiO2" >
     *   <thermo model="Lattice">
     *     <site_density units="kmol/m^3"> 73.159 </site_density>
     *   </thermo>
     * </phase>
     * @endcode
     *
     * @deprecated The XML input format is deprecated and will be removed in
     *     Cantera 3.0.
     */
    virtual void setParametersFromXML(const XML_Node& eosdata);
    //@}

protected:
    virtual void compositionChanged();

    //! Reference state pressure
    doublereal m_Pref;

    //! The current pressure
    /*!
     * Since the density isn't a function of pressure, but only of the mole
     * fractions, we need to independently specify the pressure. The density
     * variable which is inherited as part of the State class, m_dens, is always
     * kept current whenever T, P, or X[] change.
     */
    doublereal m_Pcurrent;

    //! Reference state enthalpies / RT
    mutable vector_fp m_h0_RT;

    //! Temporary storage for the reference state heat capacities
    mutable vector_fp m_cp0_R;

    //! Temporary storage for the reference state Gibbs energies
    mutable vector_fp m_g0_RT;

    //! Temporary storage for the reference state entropies at the current
    //! temperature
    mutable vector_fp m_s0_R;

    //! Vector of molar volumes for each species in the solution
    /**
     * Species molar volumes \f$ m^3 kmol^-1 \f$
     */
    vector_fp m_speciesMolarVolume;

    //! Site Density of the lattice solid
    /*!
     *  Currently, this is imposed as a function of T, P or composition
     *
     *  units are kmol m-3
     */
    doublereal m_site_density;

private:
    //! Update the species reference state thermodynamic functions
    /*!
     * The polynomials for the standard state functions are only reevaluated if
     * the temperature has changed.
     */
    void _updateThermo() const;
};
}

#endif
