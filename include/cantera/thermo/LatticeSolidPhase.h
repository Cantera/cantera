/**
 *  @file LatticeSolidPhase.h
 *  Header for a simple thermodynamics model of a bulk solid phase
 *  derived from ThermoPhase,
 *  assuming an ideal solution model based on a lattice of solid atoms
 *  (see \ref thermoprops and class \link Cantera::LatticeSolidPhase LatticeSolidPhase\endlink).
 */

//  Copyright 2005 California Institute of Technology

#ifndef CT_LATTICESOLID_H
#define CT_LATTICESOLID_H

#include "cantera/base/config.h"
#include "cantera/base/ct_defs.h"
#include "mix_defs.h"
#include "ThermoPhase.h"
#include "SpeciesThermo.h"
#include "LatticePhase.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

//! A phase that is comprised of a fixed additive combination of other lattice phases
/*!
 *  This is the main way %Cantera describes semiconductors and other solid phases.
 *  This %ThermoPhase object calculates its properties as a sum over other %LatticePhase objects. Each of the %LatticePhase
 *  objects is a %ThermoPhase object by itself.
 *
 *  The results from this LatticeSolidPhase model reduces to the LatticePhase model when there is one
 *  lattice phase and the molar densities of the sublattice and the molar density within the LatticeSolidPhase
 *  have the same values.
 *
 *  The mole fraction vector is redefined witin the the LatticeSolidPhase object. Each of the mole
 *  fractions sum to one on each of the sublattices.  The routine getMoleFraction() and setMoleFraction()
 *  have been redefined to use this convention.
 *
 * <HR>
 * <H2> Specification of Species Standard %State Properties </H2>
 * <HR>
 *
 *   The standard state properties are calculated in the normal way for each of the sublattices. The normal way
 *   here means that a thermodynamic polynomial in temperature is developed. Also, a constant volume approximation
 *   for the pressure dependence is assumed.  All of these properties are on a Joules per kmol of sublattice
 *   constituent basis.
 *
 * <HR>
 * <H2> Specification of Solution Thermodynamic Properties </H2>
 * <HR>

 *  The sum over the %LatticePhase objects is carried out by weighting each %LatticePhase object
 *  value with the molar density (kmol m-3) of its %LatticePhase. Then the resulting quantity is divided by
 *  the molar density of the total compound. The LatticeSolidPhase object therefore only contains a
 *  listing of the number of %LatticePhase object
 *  that comprises the solid, and it contains a value for the molar density of the entire mixture.
 *  This is the same thing as saying that
 *
 *  \f[
 *         L_i = L^{solid}   \theta_i
 *  \f]
 *
 *  \f$ L_i \f$ is the molar volume of the ith lattice. \f$ L^{solid} \f$ is the molar volume of the entire
 *  solid. \f$  \theta_i \f$ is a fixed weighting factor for the ith lattice representing the lattice
 *  stoichiometric coefficient. For this object the  \f$  \theta_i \f$ values are fixed.
 *
 *  Let's take FeS2 as an example, which may be thought of as a combination of two lattices: Fe and S lattice.
 *  The Fe sublattice has a molar density of 1 gmol cm-3. The S sublattice has a molar density of 2 gmol cm-3.
 *  We then define the LatticeSolidPhase object as having a nominal composition of FeS2, and having a
 *  molar density of 1 gmol cm-3.  All quantities pertaining to the FeS2 compound will be have weights
 *  associated with the sublattices. The Fe sublattice will have a weight of 1.0 associated with it. The
 *  S sublattice will have a weight of 2.0 associated with it.
 *
 * <HR>
 * <H3> Specification of Solution Density Properties </H3>
 * <HR>
 *
 *  Currently, molar density is not a constant within the object, even though the species molar volumes are a
 *  constant.  The basic idea is that a swelling of one of the sublattices will result in a swelling of
 *  of all of the lattices. Therefore, the molar volumes of the individual lattices are not independent of
 *  one another.
 *
 *   The molar volume of the Lattice solid is calculated from the following formula
 *
 *  \f[
 *         V = \sum_i{ \theta_i V_i^{lattice}}
 *  \f]
 *
 *  where \f$ V_i^{lattice} \f$ is the molar volume of the ith sublattice. This is calculated from the
 *  following standard formula.
 *
 *
 *  \f[
 *         V_i = \sum_k{ X_k V_k}
 *  \f]
 *
 *  where k is a species in the ith sublattice.
 *
 *  The mole fraction vector is redefined witin the the LatticeSolidPhase object. Each of the mole
 *  fractions sum to one on each of the sublattices.  The routine getMoleFraction() and setMoleFraction()
 *  have been redefined to use this convention.
 *
 *  (This object is still under construction)
 *
 */
class LatticeSolidPhase : public ThermoPhase
{
public:
    //! Base empty constructor
    LatticeSolidPhase();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    LatticeSolidPhase(const LatticeSolidPhase& right);

    //! Assignment operator
    /*!
     * @param right Object to be copied
     */
    LatticeSolidPhase& operator=(const LatticeSolidPhase& right);

    //! Destructor
    virtual ~LatticeSolidPhase();

    //! Duplication function
    /*!
     * This virtual function is used to create a duplicate of the
     * current phase. It's used to duplicate the phase when given
     * a ThermoPhase pointer to the phase.
     *
     * @return It returns a ThermoPhase pointer.
     */
    ThermoPhase* duplMyselfAsThermoPhase() const;

    //! Equation of state type flag.
    /*!
     *  Returns cLatticeSolid, listed in mix_defs.h.
     */
    virtual int eosType() const {
        return cLatticeSolid;
    }

    //! Minimum temperature for which the thermodynamic data for the species
    //! or phase are valid.
    /*!
     * If no argument is supplied, the
     * value returned will be the lowest temperature at which the
     * data for \e all species are valid. Otherwise, the value
     * will be only for species \a k. This function is a wrapper
     * that calls the species thermo minTemp function.
     *
     * @param k index of the species. Default is -1, which will return the max of the min value
     *          over all species.
     */
    virtual doublereal minTemp(size_t k = npos) const;

    //! Maximum temperature for which the thermodynamic data for the species
    //! are valid.
    /*!
     * If no argument is supplied, the
     * value returned will be the highest temperature at which the
     * data for \e all species are valid. Otherwise, the value
     * will be only for species \a k. This function is a wrapper
     * that calls the species thermo maxTemp function.
     *
     * @param k index of the species. Default is -1, which will return the min of the max value
     *          over all species.
     */
    virtual doublereal maxTemp(size_t k = npos) const;

    //! Returns the reference pressure in Pa. This function is a wrapper
    //! that calls the species thermo refPressure function.
    virtual doublereal refPressure() const ;

    //! This method returns the convention used in specification
    //! of the standard state, of which there are currently two,
    //! temperature based, and variable pressure based.
    /*!
     *  All of the thermo is determined by slave %ThermoPhase routines.
     */
    virtual int standardStateConvention() const {
        return cSS_CONVENTION_SLAVE;
    }

    //! Return the Molar Enthalpy. Units: J/kmol.
    /*!
     * The molar enthalpy is determined by the following formula, where \f$ \theta_n \f$ is the
     * lattice stoichiometric coefficient of the nth lattice
     *
     * \f[
     *   \tilde h(T,P) = {\sum_n \theta_n  \tilde h_n(T,P) }
     * \f]
     *
     * \f$ \tilde h_n(T,P) \f$ is the enthalpy of the n<SUP>th</SUP> lattice.
     *
     *  units J/kmol
     */
    virtual doublereal enthalpy_mole() const;

    //! Return the Molar Internal Energy. Units: J/kmol.
    /*!
     * The molar enthalpy is determined by the following formula, where \f$ \theta_n \f$ is the
     * lattice stoichiometric coefficient of the nth lattice
     *
     * \f[
     *   \tilde u(T,P) = {\sum_n \theta_n \tilde u_n(T,P) }
     * \f]
     *
     * \f$ \tilde u_n(T,P) \f$ is the internal energy of the n<SUP>th</SUP> lattice.
     *
     *  units J/kmol
     */
    virtual doublereal intEnergy_mole() const;

    //! Return the Molar Entropy. Units: J/kmol/K.
    /*!
     * The molar enthalpy is determined by the following formula, where \f$ \theta_n \f$ is the
     * lattice stoichiometric coefficient of the nth lattice
     *
     * \f[
     *   \tilde s(T,P) = \sum_n \theta_n \tilde s_n(T,P)
     * \f]
     *
     * \f$ \tilde s_n(T,P) \f$ is the molar entropy of the n<SUP>th</SUP> lattice.
     *
     *  units J/kmol/K
     */
    virtual doublereal entropy_mole() const;

    //! Return the Molar Gibbs energy. Units: J/kmol.
    /*!
     * The molar gibbs free energy is determined by the following formula, where \f$ \theta_n \f$ is the
     * lattice stoichiometric coefficient of the nth lattice
     *
     * \f[
     *   \tilde h(T,P) = {\sum_n \theta_n \tilde h_n(T,P) }
     * \f]
     *
     * \f$ \tilde h_n(T,P) \f$ is the enthalpy of the n<SUP>th</SUP> lattice.
     *
     *  units J/kmol
     */
    virtual doublereal gibbs_mole() const;

    //! Return the constant pressure heat capacity. Units: J/kmol/K
    /*!
     * The molar constant pressure heat capacity is determined by the following formula, where \f$ C_n \f$ is the
     * lattice molar density of the nth lattice, and \f$ C_T \f$ is the molar density
     * of the solid compound.
     *
     * \f[
     *   \tilde c_{p,n}(T,P) = \frac{\sum_n C_n \tilde c_{p,n}(T,P) }{C_T},
     * \f]
     *
     * \f$ \tilde c_{p,n}(T,P) \f$ is the heat capacity of the n<SUP>th</SUP> lattice.
     *
     *  units J/kmol/K
     */
    virtual doublereal cp_mole() const;

    //! Return the constant volume heat capacity. Units: J/kmol/K
    /*!
     * The molar constant volume heat capacity is determined by the following formula, where \f$ C_n \f$ is the
     * lattice molar density of the nth lattice, and \f$ C_T \f$ is the molar density
     * of the solid compound.
     *
     * \f[
     *   \tilde c_{v,n}(T,P) = \frac{\sum_n C_n \tilde c_{v,n}(T,P) }{C_T},
     * \f]
     *
     * \f$ \tilde c_{v,n}(T,P) \f$ is the heat capacity of the n<SUP>th</SUP> lattice.
     *
     *  units J/kmol/K
     */
    virtual doublereal cv_mole() const {
        return cp_mole();
    }

    //! Report the Pressure. Units: Pa.
    /*!
     *  This method simply returns the stored pressure value.
     */
    virtual doublereal pressure() const {
        return m_press;
    }

    //! Set the pressure at constant temperature. Units: Pa.
    /*!
     *
     * @param p Pressure (units - Pa)
     */
    virtual void setPressure(doublereal p);

    //! Calculate the density of the solid mixture
    /*!
     * The formula for this is
     *
     * \f[
     *      \rho = \sum_n{ \rho_n \theta_n }
     * \f]
     *
     * where \f$ \rho_n \f$  is the density of the nth sublattice
     */
    doublereal calcDensity();

    //! Set the mole fractions to the specified values, and then
    //! normalize them so that they sum to 1.0 for each of the subphases
    /*!
     *  On input, the mole fraction vector is assumed to sum to one for each of the sublattices. The sublattices
     *  are updated with this mole fraction vector. The mole fractions are also stored within this object, after
     *  they are normalized to one by dividing by the number of sublattices.
     *
     *    @param x  Input vector of mole fractions. There is no restriction
     *           on the sum of the mole fraction vector. Internally,
     *           this object will pass portions of this vector to the sublattices which assume that the portions
     *           individually sum to one.
     *           Length is m_kk.
     */
    virtual void setMoleFractions(const doublereal* const x);

    //! Get the species mole fraction vector.
    /*!
     * On output the mole fraction vector will sum to one for each of the subphases which make up this phase.
     *
     * @param x On return, x contains the mole fractions. Must have a
     *          length greater than or equal to the number of species.
     */
    virtual void getMoleFractions(doublereal* const x) const;

    //! The mole fraction of species k.
    /*!
     *   If k is outside the valid
     *   range, an exception will be thrown. Note that it is
     *   somewhat more efficient to call getMoleFractions if the
     *   mole fractions of all species are desired.
     *   @param k species index
     */
    doublereal moleFraction(const int k) const {
        return err("not implemented");
    }

    //! Get the species mass fractions.
    /*!
     * @param y On return, y contains the mass fractions. Array \a y must have a length
     *          greater than or equal to the number of species.
     */
    void getMassFractions(doublereal* const y) const {
        err("not implemented");
    }

    //! Mass fraction of species k.
    /*!
     *  If k is outside the valid range, an exception will be thrown. Note that it is
     *  somewhat more efficient to call getMassFractions if the mass fractions of all species are desired.
     *
     * @param k    species index
     */
    doublereal massFraction(const int k) const {
        return err("not implemented");
    }

    //! Set the mass fractions to the specified values, and then
    //! normalize them so that they sum to 1.0.
    /*!
     * @param y  Array of unnormalized mass fraction values (input).
     *           Must have a length greater than or equal to the number of species.
     *           Input vector of mass fractions. There is no restriction
     *           on the sum of the mass fraction vector. Internally,
     *           the State object will normalize this vector before
     *           storing its contents.
     *           Length is m_kk.
     */
    virtual void setMassFractions(const doublereal* const y) {
        err("not implemented");
    }

    //! Set the mass fractions to the specified values without normalizing.
    /*!
     * This is useful when the normalization
     * condition is being handled by some other means, for example
     * by a constraint equation as part of a larger set of equations.
     *
     * @param y  Input vector of mass fractions.
     *           Length is m_kk.
     */
    virtual void setMassFractions_NoNorm(const doublereal* const y) {
        err("not implemented");
    }

    void getConcentrations(doublereal* const c) const {
        err("not implemented");
    }

    doublereal concentration(int k) const {
        return err("not implemented");
    }

    virtual void setConcentrations(const doublereal* const conc) {
        err("not implemented");
    }

    //! This method returns an array of generalized activity concentrations
    /*!
     *  The generalized activity concentrations,
     * \f$ C^a_k \f$,  are defined such that \f$ a_k = C^a_k /
     * C^0_k, \f$ where \f$ C^0_k \f$ is a standard concentration
     * defined below and \f$ a_k \f$ are activities used in the
     * thermodynamic functions.  These activity (or generalized)
     * concentrations are used by kinetics manager classes to compute the forward and
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

    //! Get the array of non-dimensional molar-based activity coefficients at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     * @param ac Output vector of activity coefficients. Length: m_kk.
     */
    virtual void getActivityCoefficients(doublereal* ac) const;

    //! Get the species chemical potentials. Units: J/kmol.
    /*!
     * This function returns a vector of chemical potentials of the
     * species in solution at the current temperature, pressure
     * and mole fraction of the solution.
     *
     *  This returns the underlying lattice chemical potentials, as the units are kmol-1 of
     *  the sublattice species.
     *
     * @param mu  Output vector of species chemical
     *            potentials. Length: m_kk. Units: J/kmol
     */
    virtual void getChemPotentials(doublereal* mu) const;

    //! Returns an array of partial molar enthalpies for the species in the mixture.
    /*!
     * Units (J/kmol)
     * For this phase, the partial molar enthalpies are equal to the
     * pure species enthalpies
     *  \f[
     * \bar h_k(T,P) = \hat h^{ref}_k(T) + (P - P_{ref}) \hat V^0_k
     * \f]
     * The reference-state pure-species enthalpies, \f$ \hat h^{ref}_k(T) \f$,
     * at the reference pressure,\f$ P_{ref} \f$,
     * are computed by the species thermodynamic
     * property manager. They are polynomial functions of temperature.
     * @see SpeciesThermo
     *
     * @param hbar Output vector containing partial molar enthalpies.
     *             Length: m_kk.
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;

    /**
     * Returns an array of partial molar entropies of the species in the
     * solution. Units: J/kmol/K.
     * For this phase, the partial molar entropies are equal to the
     * pure species entropies plus the ideal solution contribution.
     *  \f[
     * \bar s_k(T,P) =  \hat s^0_k(T) - R log(X_k)
     * \f]
     * The reference-state pure-species entropies,\f$ \hat s^{ref}_k(T) \f$,
     * at the reference pressure, \f$ P_{ref} \f$,  are computed by the
     * species thermodynamic
     * property manager. They are polynomial functions of temperature.
     * @see SpeciesThermo
     *
     * @param sbar Output vector containing partial molar entropies.
     *             Length: m_kk.
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;

    /**
     * Returns an array of partial molar Heat Capacities at constant
     * pressure of the species in the
     * solution. Units: J/kmol/K.
     * For this phase, the partial molar heat capacities are equal
     * to the standard state heat capacities.
     *
     * @param cpbar  Output vector of partial heat capacities. Length: m_kk.
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const;

    /**
     * returns an array of partial molar volumes of the species
     * in the solution. Units: m^3 kmol-1.
     *
     * For this solution, thepartial molar volumes are equal to the
     * constant species molar volumes.
     *
     * @param vbar  Output vector of partial molar volumes. Length: m_kk.
     */
    virtual void getPartialMolarVolumes(doublereal* vbar) const;

    //! Get the array of standard state chemical potentials at unit activity for the species
    //! at their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * These are the standard state chemical potentials \f$ \mu^0_k(T,P)
     * \f$. The values are evaluated at the current
     * temperature and pressure of the solution.
     *
     *  This returns the underlying lattice standard chemical potentials, as the units are kmol-1 of
     *  the sublattice species.
     *
     * @param mu0    Output vector of chemical potentials.
     *                Length: m_kk. Units: J/kmol
     */
    virtual void getStandardChemPotentials(doublereal* mu0) const;

    //! Return the standard concentration for the kth species
    /*!
     * The standard concentration \f$ C^0_k \f$ used to normalize
     * the activity (i.e., generalized) concentration. In many cases, this quantity
     * will be the same for all species in a phase - for example,
     * for an ideal gas \f$ C^0_k = P/\hat R T \f$. For this
     * reason, this method returns a single value, instead of an
     * array.  However, for phases in which the standard
     * concentration is species-specific (e.g. surface species of
     * different sizes), this method may be called with an
     * optional parameter indicating the species.
     *
     * @param k Optional parameter indicating the species. The default
     *          is to assume this refers to species 0.
     * @return
     *   Returns the standard concentration. The units are by definition
     *   dependent on the ThermoPhase and kinetics manager representation.
     */
    virtual doublereal standardConcentration(size_t k=0) const;

    //! Natural logarithm of the standard concentration of the kth species.
    /*!
     * @param k    index of the species (defaults to zero)
     */
    virtual doublereal logStandardConc(size_t k=0) const;
    //@}
    /// @name Thermodynamic Values for the Species Reference States
    //@{

    //! Returns the vector of nondimensional enthalpies of the reference state at the current
    //!  temperature of the solution and the reference pressure for the species.
    /*!
     *  This function fills in its one entry in hrt[] by calling
     *  the underlying species thermo function for the
     *  dimensionless gibbs free energy, calculated from the
     *  dimensionless enthalpy and entropy.
     *
     *  @param grt  Vector of dimensionless Gibbs free energies of the reference state
     *              length = m_kk
     */
    virtual void getGibbs_RT_ref(doublereal* grt) const;

    //!  Returns the vector of the  gibbs function of the reference state at the current
    //! temperatureof the solution and the reference pressure for the species.
    /*!
     *  units = J/kmol
     *
     *  This function fills in its one entry in g[] by calling the underlying species thermo
     *  functions for the  gibbs free energy, calculated from enthalpy and the
     *  entropy, and the multiplying by RT.
     *
     *  @param g  Vector of Gibbs free energies of the reference state.
     *              length = m_kk
     */
    virtual void  getGibbs_ref(doublereal* g) const;

    //! Initialize the ThermoPhase object after all species have been set up
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

    //! Initialize vectors that depend on the number of species and sublattices
    void initLengths();

    //! Add in species from Slave phases
    /*!
     *  This hook is used for  cSS_CONVENTION_SLAVE phases
     *
     *  @param  phaseNode    XML_Node for the current phase
     */
    virtual void installSlavePhases(Cantera::XML_Node* phaseNode);

    //! Set equation of state parameter values from XML entries.
    /*!
     * This method is called by function importPhase() in
     * file importCTML.cpp when processing a phase definition in
     * an input file. It should be overloaded in subclasses to set
     * any parameters that are specific to that particular phase
     * model. Note, this method is called before the phase is
     * initialized with elements and/or species.
     *
     * @param eosdata An XML_Node object corresponding to
     *                the "thermo" entry for this phase in the input file.
     */
    virtual void setParametersFromXML(const XML_Node& eosdata);

    //! Set the Lattice mole fractions using a string
    /*!
     *  @param n     Integer value of the lattice whose mole fractions are being set
     *  @param x     string containing Name:value pairs that will specify the mole fractions
     *               of species on a particular lattice
     */
    void setLatticeMoleFractionsByName(int n, const std::string& x);

#ifdef H298MODIFY_CAPABILITY

    //! Modify the value of the 298 K Heat of Formation of one species in the phase (J kmol-1)
    /*!
     *   The 298K heat of formation is defined as the enthalpy change to create the standard state
     *   of the species from its constituent elements in their standard states at 298 K and 1 bar.
     *
     *   @param  k           Species k
     *   @param  Hf298New    Specify the new value of the Heat of Formation at 298K and 1 bar
     */
    virtual void modifyOneHf298SS(const size_t& k, const doublereal Hf298New);
#endif

private:
    //! error routine
    /*!
     *  @param msg Message
     *
     *  @return nothing
     */
    doublereal err(const std::string& msg) const;

protected:
    //! Last temperature at which the reference thermo was calculated
    mutable doublereal  m_tlast;

    //! Current value of the pressure
    doublereal m_press;

    //! Current value of the molar density
    doublereal m_molar_density;

    //! Number of sublattice phases
    size_t m_nlattice;

    //! Vector of sublattic ThermoPhase objects
    std::vector<LatticePhase*> m_lattice;

    //! Vector of mole fractions
    /*!
     *  Note these mole fractions sum to one when summed over all phases.
     *  However, this is not what's passed down to the lower m_lattice objects.
     */
    mutable vector_fp m_x;

    //! Lattice stoichiometric coefficients
    std::vector<doublereal> theta_;

    //! Temporary vector
    mutable vector_fp tmpV_;

    std::vector<size_t> lkstart_;

private:
    //! Update the reference thermodynamic functions
    void _updateThermo() const;
};
}

#endif
