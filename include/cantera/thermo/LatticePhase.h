/**
 *  @file LatticePhase.h Header for a simple thermodynamics model of a bulk
 *      phase derived from ThermoPhase, assuming a lattice of solid atoms (see
 *      @ref thermoprops and class @link Cantera::LatticePhase
 *      LatticePhase@endlink).
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
 * The density of matrix sites is given by the variable @f$ C_o @f$, which has
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
 * which has a weak dependence on the system pressure, @f$ P @f$.
 *
 * @f[
 *     h^o_k(T,P) =
 *            h^{ref}_k(T) +  \left( \frac{P - P_{ref}}{C_o} \right)
 * @f]
 *
 * For an incompressible substance, the molar internal energy is independent of
 * pressure. Since the thermodynamic properties are specified by giving the
 * standard-state enthalpy, the term @f$ \frac{P_{ref}}{C_o} @f$ is subtracted
 * from the specified reference molar enthalpy to compute the standard state
 * molar internal energy:
 *
 * @f[
 *      u^o_k(T,P) = h^{ref}_k(T) - \frac{P_{ref}}{C_o}
 * @f]
 *
 * The standard state heat capacity, internal energy, and entropy are
 * independent of pressure. The standard state Gibbs free energy is obtained
 * from the enthalpy and entropy functions.
 *
 * The standard state molar volume is independent of temperature, pressure, and
 * species identity:
 *
 * @f[
 *      V^o_k(T,P) = \frac{1.0}{C_o}
 * @f]
 *
 * ## Specification of Solution Thermodynamic Properties
 *
 * The activity of species @f$ k @f$ defined in the phase, @f$ a_k @f$, is given
 * by the ideal solution law:
 *
 * @f[
 *      a_k = X_k ,
 * @f]
 *
 * where @f$ X_k @f$ is the mole fraction of species *k*. The chemical potential
 * for species *k* is equal to
 *
 * @f[
 *      \mu_k(T,P) = \mu^o_k(T, P) + R T \ln X_k
 * @f]
 *
 * The partial molar entropy for species *k* is given by the following relation,
 *
 * @f[
 *      \tilde{s}_k(T,P) = s^o_k(T,P) - R \ln X_k = s^{ref}_k(T) - R \ln X_k
 * @f]
 *
 * The partial molar enthalpy for species *k* is
 *
 * @f[
 *      \tilde{h}_k(T,P) = h^o_k(T,P) = h^{ref}_k(T) + \left( \frac{P - P_{ref}}{C_o} \right)
 * @f]
 *
 * The partial molar Internal Energy for species *k* is
 *
 * @f[
 *      \tilde{u}_k(T,P) = u^o_k(T,P) = u^{ref}_k(T)
 * @f]
 *
 * The partial molar Heat Capacity for species *k* is
 *
 * @f[
 *      \tilde{Cp}_k(T,P) = Cp^o_k(T,P) = Cp^{ref}_k(T)
 * @f]
 *
 * The partial molar volume is independent of temperature, pressure, and species
 * identity:
 *
 * @f[
 *      \tilde{V}_k(T,P) =  V^o_k(T,P) = \frac{1.0}{C_o}
 * @f]
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
 * ## Application within Kinetics Managers
 *
 * @f$ C^a_k @f$ are defined such that @f$ C^a_k = a_k = X_k @f$. @f$ C^s_k @f$,
 * the standard concentration, is defined to be equal to one. @f$ a_k @f$ are
 * activities used in the thermodynamic functions.  These activity (or
 * generalized) concentrations are used by kinetics manager classes to compute
 * the forward and reverse rates of elementary reactions. The activity
 * concentration,@f$  C^a_k @f$, is given by the following expression.
 *
 * @f[
 *      C^a_k = C^s_k  X_k  =  X_k
 * @f]
 *
 * The standard concentration for species *k* is identically one
 *
 * @f[
 *     C^s_k =  C^s = 1.0
 * @f]
 *
 * For example, a bulk-phase binary gas reaction between species j and k,
 * producing a new species l would have the following equation for its rate of
 * progress variable, @f$ R^1 @f$, which has units of kmol m-3 s-1.
 *
 * @f[
 *    R^1 = k^1 C_j^a C_k^a =  k^1  X_j X_k
 * @f]
 *
 * The reverse rate constant can then be obtained from the law of microscopic
 * reversibility and the equilibrium expression for the system.
 *
 * @f[
 *     \frac{X_j X_k}{ X_l} = K_a^{o,1} = \exp(\frac{\mu^o_l - \mu^o_j - \mu^o_k}{R T} )
 * @f]
 *
 * @f$  K_a^{o,1} @f$ is the dimensionless form of the equilibrium constant,
 * associated with the pressure dependent standard states @f$ \mu^o_l(T,P) @f$
 * and their associated activities,
 * @f$ a_l @f$, repeated here:
 *
 * @f[
 *      \mu_l(T,P) = \mu^o_l(T, P) + R T \ln a_l
 * @f]
 *
 * The concentration equilibrium constant, @f$ K_c @f$, may be obtained by
 * changing over to activity concentrations. When this is done:
 *
 * @f[
 *     \frac{C^a_j C^a_k}{ C^a_l} = C^o K_a^{o,1} = K_c^1 =
 *         \exp(\frac{\mu^{o}_l - \mu^{o}_j - \mu^{o}_k}{R T} )
 * @f]
 *
 * %Kinetics managers will calculate the concentration equilibrium constant, @f$
 * K_c @f$, using the second and third part of the above expression as a
 * definition for the concentration equilibrium constant.
 *
 * @ingroup thermoprops
 */
class LatticePhase : public ThermoPhase
{
public:
    //! Full constructor for a lattice phase
    /*!
     * @param inputFile String name of the input file. If blank,
     *                  an empty phase will be created.
     * @param id        string id of the phase name
     */
    explicit LatticePhase(const string& inputFile="", const string& id="");

    string type() const override {
        return "lattice";
    }

    bool isCompressible() const override {
        return false;
    }

    map<string, size_t> nativeState() const override {
        return { {"T", 0}, {"P", 1}, {"X", 2} };
    }

    //! @name Molar Thermodynamic Properties of the Solution
    //! @{

    //! Return the Molar Enthalpy. Units: J/kmol.
    /*!
     * For an ideal solution,
     *
     *   @f[
     *    \hat h(T,P) = \sum_k X_k \hat h^0_k(T,P),
     *   @f]
     *
     * The standard-state pure-species Enthalpies @f$ \hat h^0_k(T,P) @f$ are
     * computed first by the species reference state thermodynamic property
     * manager and then a small pressure dependent term is added in.
     *
     * \see MultiSpeciesThermo
     */
    double enthalpy_mole() const override;

    //! Molar entropy of the solution. Units: J/kmol/K
    /*!
     * For an ideal, constant partial molar volume solution mixture with
     * pure species phases which exhibit zero volume expansivity:
     *   @f[
     *      \hat s(T, P, X_k) = \sum_k X_k \hat s^0_k(T)  - \hat R  \sum_k X_k \ln(X_k)
     *   @f]
     * The reference-state pure-species entropies @f$ \hat s^0_k(T,p_{ref}) @f$
     * are computed by the species thermodynamic property manager. The pure
     * species entropies are independent of pressure since the volume
     * expansivities are equal to zero.
     *
     * Units: J/kmol/K.
     *
     * @see MultiSpeciesThermo
     */
    double entropy_mole() const override;

    //! Molar heat capacity at constant pressure of the solution.
    //! Units: J/kmol/K.
    /*!
     * For an ideal, constant partial molar volume solution mixture with
     * pure species phases which exhibit zero volume expansivity:
     *   @f[
     *    \hat c_p(T,P) = \sum_k X_k \hat c^0_{p,k}(T) .
     *   @f]
     * The heat capacity is independent of pressure. The reference-state pure-
     * species heat capacities @f$ \hat c^0_{p,k}(T) @f$ are computed by the
     * species thermodynamic property manager.
     *
     * @see MultiSpeciesThermo
     */
    double cp_mole() const override;

    //! Molar heat capacity at constant volume of the solution.
    //! Units: J/kmol/K.
    /*!
     * For an ideal, constant partial molar volume solution mixture with
     * pure species phases which exhibit zero volume expansivity:
     *  @f[
     *     \hat c_v(T,P) = \hat c_p(T,P)
     *  @f]
     *
     * The two heat capacities are equal.
     */
    double cv_mole() const override;

    //! @}
    //! @name Mechanical Equation of State Properties
    //!
    //! In this equation of state implementation, the density is a function only
    //! of the mole fractions. Therefore, it can't be an independent variable.
    //! Instead, the pressure is used as the independent variable. Functions
    //! which try to set the thermodynamic state by calling setDensity() may
    //! cause an exception to be thrown.
    //! @{

    //! Pressure. Units: Pa.
    /*!
     * For this incompressible system, we return the internally stored
     * independent value of the pressure.
     */
    double pressure() const override {
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
    void setPressure(double p) override;

    //! Calculate the density of the mixture using the partial molar volumes and
    //! mole fractions as input
    /*!
     * The formula for this is
     *
     * @f[
     *      \rho = \frac{\sum_k{X_k W_k}}{\sum_k{X_k V_k}}
     * @f]
     *
     * where @f$ X_k @f$ are the mole fractions, @f$ W_k @f$ are the molecular
     * weights, and @f$ V_k @f$ are the pure species molar volumes.
     *
     * Note, the basis behind this formula is that in an ideal solution the
     * partial molar volumes are equal to the pure species molar volumes. We
     * have additionally specified in this class that the pure species molar
     * volumes are independent of temperature and pressure.
     */
    double calcDensity();

    //! @}
    //! @name Activities, Standard States, and Activity Concentrations
    //!
    //! The activity @f$ a_k @f$ of a species in solution is related to the
    //! chemical potential by @f[ \mu_k = \mu_k^0(T) + \hat R T \ln a_k. @f] The
    //! quantity @f$ \mu_k^0(T,P) @f$ is the chemical potential at unit activity,
    //! which depends only on temperature and the pressure. Activity is assumed
    //! to be molality-based here.
    //! @{

    Units standardConcentrationUnits() const override;
    void getActivityConcentrations(double* c) const override;

    //! Return the standard concentration for the kth species
    /*!
     * The standard concentration @f$ C^0_k @f$ used to normalize
     * the activity (that is, generalized) concentration for use
     *
     * For the time being, we will use the concentration of pure solvent for the
     * the standard concentration of all species. This has the effect of making
     * mass-action reaction rates based on the molality of species proportional
     * to the molality of the species.
     *
     * @param k Optional parameter indicating the species. The default is to
     *         assume this refers to species 0.
     * @returns the standard Concentration in units of m^3/kmol.
     */
    double standardConcentration(size_t k=0) const override;
    double logStandardConc(size_t k=0) const override;

    //! Get the array of non-dimensional activity coefficients at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     *  For this phase, the activity coefficients are all equal to one.
     *
     * @param ac Output vector of activity coefficients. Length: m_kk.
     */
    void getActivityCoefficients(double* ac) const override;

    //! @}
    //! @name  Partial Molar Properties of the Solution
    //! @{

    //! Get the species chemical potentials. Units: J/kmol.
    /*!
     * This function returns a vector of chemical potentials of the species in
     * solid solution at the current temperature, pressure and mole fraction of
     * the solid solution.
     *
     * @param mu  Output vector of species chemical
     *            potentials. Length: m_kk. Units: J/kmol
     */
    void getChemPotentials(double* mu) const override;

    /**
     * Returns an array of partial molar enthalpies for the species in the
     * mixture. Units (J/kmol). For this phase, the partial molar enthalpies are
     * equal to the pure species enthalpies
     * @f[
     * \bar h_k(T,P) = \hat h^{ref}_k(T) + (P - P_{ref}) \hat V^0_k
     * @f]
     * The reference-state pure-species enthalpies, @f$ \hat h^{ref}_k(T) @f$,
     * at the reference pressure,@f$ P_{ref} @f$, are computed by the species
     * thermodynamic property manager. They are polynomial functions of
     * temperature.
     * @see MultiSpeciesThermo
     *
     * @param hbar Output vector containing partial molar enthalpies.
     *             Length: m_kk.
     */
    void getPartialMolarEnthalpies(double* hbar) const override;

    /**
     * Returns an array of partial molar entropies of the species in the
     * solution. Units: J/kmol/K. For this phase, the partial molar entropies
     * are equal to the pure species entropies plus the ideal solution
     * contribution.
     * @f[
     * \bar s_k(T,P) =  \hat s^0_k(T) - R \ln(X_k)
     * @f]
     * The reference-state pure-species entropies,@f$ \hat s^{ref}_k(T) @f$, at
     * the reference pressure, @f$ P_{ref} @f$, are computed by the species
     * thermodynamic property manager. They are polynomial functions of
     * temperature.
     * @see MultiSpeciesThermo
     *
     * @param sbar Output vector containing partial molar entropies.
     *             Length: m_kk.
     */
    void getPartialMolarEntropies(double* sbar) const override;

    /**
     * Returns an array of partial molar Heat Capacities at constant pressure of
     * the species in the solution. Units: J/kmol/K. For this phase, the partial
     * molar heat capacities are equal to the standard state heat capacities.
     *
     * @param cpbar  Output vector of partial heat capacities. Length: m_kk.
     */
    void getPartialMolarCp(double* cpbar) const override;

    void getPartialMolarVolumes(double* vbar) const override;
    void getStandardChemPotentials(double* mu) const override;
    void getPureGibbs(double* gpure) const override;

    //! @}
    //! @name  Properties of the Standard State of the Species in the Solution
    //! @{

    //! Get the nondimensional Enthalpy functions for the species standard
    //! states at their standard states at the current *T* and *P* of the
    //! solution.
    /*!
     * A small pressure dependent term is added onto the reference state enthalpy
     * to get the pressure dependence of this term.
     *
     * @f[
     *    h^o_k(T,P) = h^{ref}_k(T) +  \left( \frac{P - P_{ref}}{C_o} \right)
     * @f]
     *
     * The reference state thermodynamics is obtained by a pointer to a
     * populated species thermodynamic property manager class (see
     * ThermoPhase::m_spthermo). How to relate pressure changes to the reference
     * state thermodynamics is resolved at this level.
     *
     * @param hrt      Output vector of nondimensional standard state enthalpies.
     *                 Length: m_kk.
     */
    void getEnthalpy_RT(double* hrt) const override;

    //! Get the array of nondimensional Entropy functions for the species
    //! standard states at the current *T* and *P* of the solution.
    /*!
     * The entropy of the standard state is defined as independent of
     * pressure here.
     *
     * @f[
     *      s^o_k(T,P) = s^{ref}_k(T)
     * @f]
     *
     * The reference state thermodynamics is obtained by a pointer to a
     * populated species thermodynamic property manager class (see
     * ThermoPhase::m_spthermo). How to relate pressure changes to the reference
     * state thermodynamics is resolved at this level.
     *
     * @param sr   Output vector of nondimensional standard state entropies.
     *             Length: m_kk.
     */
    void getEntropy_R(double* sr) const override;

    //! Get the nondimensional Gibbs functions for the species standard states
    //! at the current *T* and *P* of the solution.
    /*!
     * The standard Gibbs free energies are obtained from the enthalpy and
     * entropy formulation.
     *
     * @f[
     *      g^o_k(T,P) = h^{o}_k(T,P) - T s^{o}_k(T,P)
     * @f]
     *
     * @param grt  Output vector of nondimensional standard state Gibbs free
     *             energies. Length: m_kk.
     */
    void getGibbs_RT(double* grt) const override;

    //! Get the nondimensional Heat Capacities at constant pressure for the
    //! species standard states at the current *T* and *P* of the solution
    /*!
     * The heat capacity of the standard state is independent of pressure
     *
     * @f[
     *      Cp^o_k(T,P) = Cp^{ref}_k(T)
     * @f]
     *
     * The reference state thermodynamics is obtained by a pointer to a
     * populated species thermodynamic property manager class (see
     * ThermoPhase::m_spthermo). How to relate pressure changes to the reference
     * state thermodynamics is resolved at this level.
     *
     * @param cpr   Output vector of nondimensional standard state heat
     *              capacities. Length: m_kk.
     */
    void getCp_R(double* cpr) const override;

    void getStandardVolumes(double* vol) const override;

    //! @}
    //! @name Thermodynamic Values for the Species Reference States
    //! @{

    const vector<double>& enthalpy_RT_ref() const;

    //! Returns a reference to the dimensionless reference state Gibbs free
    //! energy vector.
    /*!
     * This function is part of the layer that checks/recalculates the reference
     * state thermo functions.
     */
    const vector<double>& gibbs_RT_ref() const;

    void getGibbs_RT_ref(double* grt) const override;
    void getGibbs_ref(double* g) const override;

    //! Returns a reference to the dimensionless reference state Entropy vector.
    /*!
     * This function is part of the layer that checks/recalculates the reference
     * state thermo functions.
     */
    const vector<double>& entropy_R_ref() const;

    //! Returns a reference to the dimensionless reference state Heat Capacity
    //! vector.
    /*!
     * This function is part of the layer that checks/recalculates the reference
     * state thermo functions.
     */
    const vector<double>& cp_R_ref() const;

    //! @}
    //! @name  Utilities for Initialization of the Object
    //! @{

    bool addSpecies(shared_ptr<Species> spec) override;

    //! Set the density of lattice sites [kmol/m^3]
    void setSiteDensity(double sitedens);

    void initThermo() override;
    void getParameters(AnyMap& phaseNode) const override;
    void getSpeciesParameters(const string& name, AnyMap& speciesNode) const override;

    //! @}

protected:
    void compositionChanged() override;

    //! Reference state pressure
    double m_Pref = OneAtm;

    //! The current pressure
    /*!
     * Since the density isn't a function of pressure, but only of the mole
     * fractions, we need to independently specify the pressure. The density
     * variable which is inherited as part of the State class, m_dens, is always
     * kept current whenever T, P, or X[] change.
     */
    double m_Pcurrent = OneAtm;

    //! Reference state enthalpies / RT
    mutable vector<double> m_h0_RT;

    //! Temporary storage for the reference state heat capacities
    mutable vector<double> m_cp0_R;

    //! Temporary storage for the reference state Gibbs energies
    mutable vector<double> m_g0_RT;

    //! Temporary storage for the reference state entropies at the current
    //! temperature
    mutable vector<double> m_s0_R;

    //! Vector of molar volumes for each species in the solution
    /**
     * Species molar volumes @f$ m^3 kmol^-1 @f$
     */
    vector<double> m_speciesMolarVolume;

    //! Site Density of the lattice solid
    /*!
     *  Currently, this is imposed as a function of T, P or composition
     *
     *  units are kmol m-3
     */
    double m_site_density = 0.0;

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
