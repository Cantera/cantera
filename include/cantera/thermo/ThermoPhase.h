/**
 *  @file ThermoPhase.h
 * Header file for class ThermoPhase, the base class for phases with
 * thermodynamic properties, and the text for the Module thermoprops
 * (see \ref thermoprops and class \link Cantera::ThermoPhase ThermoPhase\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_THERMOPHASE_H
#define CT_THERMOPHASE_H

#include "Phase.h"
#include "MultiSpeciesThermo.h"
#include "cantera/base/global.h"

namespace Cantera
{

/*!
 * @name CONSTANTS - Specification of the Molality convention
 */
//@{
//! Standard state uses the molar convention
const int cAC_CONVENTION_MOLAR = 0;
//! Standard state uses the molality convention
const int cAC_CONVENTION_MOLALITY = 1;
//@}

/*!
 * @name CONSTANTS - Specification of the SS convention
 */
//@{
//! Standard state uses the molar convention
const int cSS_CONVENTION_TEMPERATURE = 0;
//! Standard state uses the molality convention
const int cSS_CONVENTION_VPSS = 1;
//! Standard state thermodynamics is obtained from slave ThermoPhase objects
const int cSS_CONVENTION_SLAVE = 2;
//@}

//! Base class for a phase with thermodynamic properties.
/*!
 * Class ThermoPhase is the base class for the family of classes that represent
 * phases of matter of any type. It defines a common public interface, and
 * implements a few methods. Most of the methods, however, are declared virtual
 * and are meant to be overloaded in derived classes.  The standard way used
 * throughout Cantera to compute properties of phases of matter is through
 * pointers of type ThermoPhase* that point to objects of subclasses of
 * ThermoPhase.
 *
 * Class ThermoPhase extends class Phase by adding methods to compute
 * thermodynamic properties in addition to the ones that are used to define the
 * state of a substance (temperature, density/pressure and composition). The
 * distinction is that the methods declared in ThermoPhase require knowing the
 * particular equation of state of the phase of interest, while those of class
 * Phase do not, since they only involve data values stored within the object.
 *
 * Instances of subclasses of ThermoPhase should be created using the factory
 * class ThermoFactory, not by calling the constructor directly. This allows new
 * classes to be used with the various Cantera language interfaces.
 *
 * To implement a new equation of state, derive a class from ThermoPhase and
 * overload the virtual methods in ThermoPhase. Methods that are not needed can
 * be left unimplemented, which will cause an exception to be thrown if it is
 * called.
 *
 * Relationship with the kinetics operator:
 *
 * Describe activity coefficients.
 *
 * Describe K_a, K_p, and K_c, These are three different equilibrium
 * constants.
 *
 *   K_a is the calculation of the equilibrium constant from the standard state
 *   Gibbs free energy values. It is by definition dimensionless.
 *
 *   K_p is the calculation of the equilibrium constant from the reference state
 *   Gibbs free energy values. It is by definition dimensionless. The pressure
 *   dependence is handled entirely on the RHS of the equilibrium expression.
 *
 *   K_c is the equilibrium constant calculated from the activity
 *   concentrations. The dimensions depend on the number of products and
 *   reactants.
 *
 * The kinetics manager requires the calculation of K_c for the calculation of
 * the reverse rate constant
 *
 * @ingroup thermoprops
 * @ingroup phases
 */
class ThermoPhase : public Phase
{
public:
    //! Constructor. Note that ThermoPhase is meant to be used as a base class,
    //! so this constructor should not be called explicitly.
    ThermoPhase();

    virtual ~ThermoPhase();

    //! @name  Information Methods
    //! @{

    virtual std::string type() const {
        return "ThermoPhase";
    }

    //! String indicating the mechanical phase of the matter in this Phase.
    /*!
     * Options for the string are:
     *   * `unspecified`
     *   * `supercritical`
     *   * `gas`
     *   * `liquid`
     *   * `solid`
     *   * `solid-liquid-mix`
     *   * `solid-gas-mix`
     *   * `liquid-gas-mix`
     *   * `solid-liquid-gas-mix`
     *
     * `unspecified` is the default and should be used when the Phase does not
     * distinguish between mechanical phases or does not have enough information to
     * determine which mechanical phase(s) are present.
     *
     * @todo Needs to be implemented for all phase types. Currently only implemented for
     * PureFluidPhase.
     */
    virtual std::string phaseOfMatter() const {
        return "unspecified";
    }

    /**
     * Returns the reference pressure in Pa. This function is a wrapper
     * that calls the species thermo refPressure function.
     */
    virtual doublereal refPressure() const {
        return m_spthermo.refPressure();
    }

    //! Minimum temperature for which the thermodynamic data for the species
    //! or phase are valid.
    /*!
     * If no argument is supplied, the value returned will be the lowest
     * temperature at which the data for \e all species are valid. Otherwise,
     * the value will be only for species \a k. This function is a wrapper that
     * calls the species thermo minTemp function.
     *
     * @param k index of the species. Default is -1, which will return the max
     *          of the min value over all species.
     */
    virtual doublereal minTemp(size_t k = npos) const {
        return m_spthermo.minTemp(k);
    }

    //! Report the 298 K Heat of Formation of the standard state of one species
    //! (J kmol-1)
    /*!
     *  The 298K Heat of Formation is defined as the enthalpy change to create
     *  the standard state of the species from its constituent elements in their
     *  standard states at 298 K and 1 bar.
     *
     *   @param k    species index
     *   @returns    the current value of the Heat of Formation at 298K
     *       and 1 bar
     */
    doublereal Hf298SS(const size_t k) const {
        return m_spthermo.reportOneHf298(k);
    }

    //! Modify the value of the 298 K Heat of Formation of one species in the
    //! phase (J kmol-1)
    /*!
     *  The 298K heat of formation is defined as the enthalpy change to create
     *  the standard state of the species from its constituent elements in their
     *  standard states at 298 K and 1 bar.
     *
     *   @param  k           Species k
     *   @param  Hf298New    Specify the new value of the Heat of Formation at
     *       298K and 1 bar
     */
    virtual void modifyOneHf298SS(const size_t k, const doublereal Hf298New) {
        m_spthermo.modifyOneHf298(k, Hf298New);
        invalidateCache();
    }

    //! Restore the original heat of formation of one or more species
    /*!
     *  Resets changes made by modifyOneHf298SS(). If the species index is not
     *  specified, the heats of formation for all species are restored.
     */
    virtual void resetHf298(const size_t k=npos);

    //! Maximum temperature for which the thermodynamic data for the species
    //! are valid.
    /*!
     * If no argument is supplied, the value returned will be the highest
     * temperature at which the data for \e all species are valid. Otherwise,
     * the value will be only for species \a k. This function is a wrapper that
     * calls the species thermo maxTemp function.
     *
     * @param k index of the species. Default is -1, which will return the min
     *          of the max value over all species.
     */
    virtual doublereal maxTemp(size_t k = npos) const {
        return m_spthermo.maxTemp(k);
    }

    //! Returns the chargeNeutralityNecessity boolean
    /*!
     * Some phases must have zero net charge in order for their thermodynamics
     * functions to be valid. If this is so, then the value returned from this
     * function is true. If this is not the case, then this is false. Now, ideal
     * gases have this parameter set to false, while solution with molality-
     * based activity coefficients have this parameter set to true.
     */
    bool chargeNeutralityNecessary() const {
        return m_chargeNeutralityNecessary;
    }

    //! @}
    //! @name  Molar Thermodynamic Properties of the Solution
    //! @{

    /// Molar enthalpy. Units: J/kmol.
    virtual doublereal enthalpy_mole() const {
        throw NotImplementedError("ThermoPhase::enthalpy_mole");
    }

    /// Molar internal energy. Units: J/kmol.
    virtual doublereal intEnergy_mole() const {
        return enthalpy_mole() - pressure()* molarVolume();
    }

    /// Molar entropy. Units: J/kmol/K.
    virtual doublereal entropy_mole() const {
        throw NotImplementedError("ThermoPhase::entropy_mole");
    }

    /// Molar Gibbs function. Units: J/kmol.
    virtual doublereal gibbs_mole() const {
        return enthalpy_mole() - temperature()*entropy_mole();
    }

    /// Molar heat capacity at constant pressure. Units: J/kmol/K.
    virtual doublereal cp_mole() const {
        throw NotImplementedError("ThermoPhase::cp_mole");
    }

    /// Molar heat capacity at constant volume. Units: J/kmol/K.
    virtual doublereal cv_mole() const {
        throw NotImplementedError("ThermoPhase::cv_mole");
    }

    //! @}
    //! @name Mechanical Properties
    //! @{

    //! Returns the isothermal compressibility. Units: 1/Pa.
    /*!
     * The isothermal compressibility is defined as
     * \f[
     * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
     * \f]
     *  or
     * \f[
     * \kappa_T = \frac{1}{\rho}\left(\frac{\partial \rho}{\partial P}\right)_T
     * \f]
     */
    virtual doublereal isothermalCompressibility() const {
        throw NotImplementedError("ThermoPhase::isothermalCompressibility");
    }

    //! Return the volumetric thermal expansion coefficient. Units: 1/K.
    /*!
     * The thermal expansion coefficient is defined as
     * \f[
     * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
     * \f]
     */
    virtual doublereal thermalExpansionCoeff() const {
        throw NotImplementedError("ThermoPhase::thermalExpansionCoeff");
    }

    /**
     * @}
     * @name Electric Potential
     *
     * The phase may be at some non-zero electrical potential. These methods
     * set or get the value of the electric potential.
     */
    //@{

    //! Set the electric potential of this phase (V).
    /*!
     * This is used by classes InterfaceKinetics and EdgeKinetics to
     * compute the rates of charge-transfer reactions, and in computing
     * the electrochemical potentials of the species.
     *
     * Each phase may have its own electric potential.
     *
     * @param v Input value of the electric potential in Volts
     */
    void setElectricPotential(doublereal v) {
        m_phi = v;
    }

    //! Returns the electric potential of this phase (V).
    /*!
     *  Units are Volts (which are Joules/coulomb)
     */
    doublereal electricPotential() const {
        return m_phi;
    }

    /**
     * @}
     * @name Activities, Standard States, and Activity Concentrations
     *
     * The activity \f$a_k\f$ of a species in solution is related to the
     * chemical potential by \f[ \mu_k = \mu_k^0(T,P) + \hat R T \log a_k. \f]
     * The quantity \f$\mu_k^0(T,P)\f$ is the standard chemical potential at
     * unit activity, which depends on temperature and pressure, but not on
     * composition. The activity is dimensionless.
     * @{
     */

    //! This method returns the convention used in specification of the
    //! activities, of which there are currently two, molar- and molality-based
    //! conventions.
    /*!
     * Currently, there are two activity conventions:
     *  - Molar-based activities
     *       %Unit activity of species at either a hypothetical pure
     *       solution of the species or at a hypothetical
     *       pure ideal solution at infinite dilution
     *   cAC_CONVENTION_MOLAR 0
     *      - default
     *
     *  - Molality-based activities
     *       (unit activity of solutes at a hypothetical 1 molal
     *        solution referenced to infinite dilution at all
     *        pressures and temperatures).
     *   cAC_CONVENTION_MOLALITY 1
     */
    virtual int activityConvention() const;

    //! This method returns the convention used in specification of the standard
    //! state, of which there are currently two, temperature based, and variable
    //! pressure based.
    /*!
     * Currently, there are two standard state conventions:
     *  - Temperature-based activities
     *   cSS_CONVENTION_TEMPERATURE 0
     *      - default
     *
     *  -  Variable Pressure and Temperature -based activities
     *   cSS_CONVENTION_VPSS 1
     *
     *  -  Thermodynamics is set via slave ThermoPhase objects with
     *     nothing being carried out at this ThermoPhase object level
     *   cSS_CONVENTION_SLAVE 2
     */
    virtual int standardStateConvention() const;

    //! Returns the units of the "standard concentration" for this phase
    /*!
     * These are the units of the values returned by the functions
     * getActivityConcentrations() and standardConcentration(), which can
     * vary between different ThermoPhase-derived classes, or change within
     * a single class depending on input options. See the documentation for
     * standardConcentration() for the derived class for specific details.
     */
    virtual Units standardConcentrationUnits() const;

    //! This method returns an array of generalized concentrations
    /*!
     * \f$ C^a_k\f$ are defined such that \f$ a_k = C^a_k / C^0_k, \f$ where
     * \f$ C^0_k \f$ is a standard concentration defined below and \f$ a_k \f$
     * are activities used in the thermodynamic functions. These activity (or
     * generalized) concentrations are used by kinetics manager classes to
     * compute the forward and reverse rates of elementary reactions. Note that
     * they may or may not have units of concentration --- they might be partial
     * pressures, mole fractions, or surface coverages, for example.
     *
     * @param c Output array of generalized concentrations. The units depend
     *           upon the implementation of the reaction rate expressions within
     *           the phase.
     */
    virtual void getActivityConcentrations(doublereal* c) const {
        throw NotImplementedError("ThermoPhase::getActivityConcentrations");
    }

    //! Return the standard concentration for the kth species
    /*!
     * The standard concentration \f$ C^0_k \f$ used to normalize the activity
     * (i.e., generalized) concentration. In many cases, this quantity will be
     * the same for all species in a phase - for example, for an ideal gas \f$
     * C^0_k = P/\hat R T \f$. For this reason, this method returns a single
     * value, instead of an array.  However, for phases in which the standard
     * concentration is species-specific (e.g. surface species of different
     * sizes), this method may be called with an optional parameter indicating
     * the species.
     *
     * @param k Optional parameter indicating the species. The default
     *          is to assume this refers to species 0.
     * @return
     *   Returns the standard concentration. The units are by definition
     *   dependent on the ThermoPhase and kinetics manager representation.
     */
    virtual doublereal standardConcentration(size_t k=0) const {
        throw NotImplementedError("ThermoPhase::standardConcentration");
    }

    //! Natural logarithm of the standard concentration of the kth species.
    /*!
     * @param k    index of the species (defaults to zero)
     */
    virtual doublereal logStandardConc(size_t k=0) const;

    //! Get the array of non-dimensional activities at the current solution
    //! temperature, pressure, and solution concentration.
    /*!
     * Note, for molality based formulations, this returns the molality based
     * activities.
     *
     * We resolve this function at this level by calling on the
     * activityConcentration function. However, derived classes may want to
     * override this default implementation.
     *
     * @param a   Output vector of activities. Length: m_kk.
     */
    virtual void getActivities(doublereal* a) const;

    //! Get the array of non-dimensional molar-based activity coefficients at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     * @param ac Output vector of activity coefficients. Length: m_kk.
     */
    virtual void getActivityCoefficients(doublereal* ac) const {
        if (m_kk == 1) {
            ac[0] = 1.0;
        } else {
            throw NotImplementedError("ThermoPhase::getActivityCoefficients");
        }
    }

    //! Get the array of non-dimensional molar-based ln activity coefficients at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     * @param lnac Output vector of ln activity coefficients. Length: m_kk.
     */
    virtual void getLnActivityCoefficients(doublereal* lnac) const;

    //@}
    /// @name  Partial Molar Properties of the Solution
    //@{

    /**
     * Get the array of non-dimensional species chemical potentials
     * These are partial molar Gibbs free energies.
     * \f$ \mu_k / \hat R T \f$.
     * Units: unitless
     *
     * @param mu  Output vector of dimensionless chemical potentials.
     *            Length: m_kk.
     */
    virtual void getChemPotentials_RT(doublereal* mu) const {
        throw NotImplementedError("ThermoPhase::getChemPotentials_RT");
    }

    //! Get the species chemical potentials. Units: J/kmol.
    /*!
     * This function returns a vector of chemical potentials of the species in
     * solution at the current temperature, pressure and mole fraction of the
     * solution.
     *
     * @param mu  Output vector of species chemical
     *            potentials. Length: m_kk. Units: J/kmol
     */
    virtual void getChemPotentials(doublereal* mu) const {
        throw NotImplementedError("ThermoPhase::getChemPotentials");
    }

    //!  Get the species electrochemical potentials.
    /*!
     *  These are partial molar quantities.  This method adds a term \f$ F z_k
     *  \phi_p \f$ to each chemical potential. The electrochemical potential of
     *  species k in a phase p, \f$ \zeta_k \f$, is related to the chemical
     *  potential via the following equation,
     *
     *       \f[
     *            \zeta_{k}(T,P) = \mu_{k}(T,P) + F z_k \phi_p
     *       \f]
     *
     * @param mu  Output vector of species electrochemical
     *            potentials. Length: m_kk. Units: J/kmol
     */
    void getElectrochemPotentials(doublereal* mu) const;

    //! Returns an array of partial molar enthalpies for the species
    //! in the mixture. Units (J/kmol)
    /*!
     * @param hbar    Output vector of species partial molar enthalpies.
     *                Length: m_kk. units are J/kmol.
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const {
        throw NotImplementedError("ThermoPhase::getPartialMolarEnthalpies");
    }

    //! Returns an array of partial molar entropies of the species in the
    //! solution. Units: J/kmol/K.
    /*!
     * @param sbar    Output vector of species partial molar entropies.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const {
        throw NotImplementedError("ThermoPhase::getPartialMolarEntropies");
    }

    //! Return an array of partial molar internal energies for the
    //! species in the mixture.  Units: J/kmol.
    /*!
     * @param ubar    Output vector of species partial molar internal energies.
     *                Length = m_kk. units are J/kmol.
     */
    virtual void getPartialMolarIntEnergies(doublereal* ubar) const {
        throw NotImplementedError("ThermoPhase::getPartialMolarIntEnergies");
    }

    //! Return an array of partial molar heat capacities for the
    //! species in the mixture.  Units: J/kmol/K
    /*!
     * @param cpbar   Output vector of species partial molar heat
     *                capacities at constant pressure.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const {
        throw NotImplementedError("ThermoPhase::getPartialMolarCp");
    }

    //! Return an array of partial molar volumes for the
    //! species in the mixture. Units: m^3/kmol.
    /*!
     *  @param vbar   Output vector of species partial molar volumes.
     *                Length = m_kk. units are m^3/kmol.
     */
    virtual void getPartialMolarVolumes(doublereal* vbar) const {
        throw NotImplementedError("ThermoPhase::getPartialMolarVolumes");
    }

    //@}
    /// @name Properties of the Standard State of the Species in the Solution
    //@{

    //! Get the array of chemical potentials at unit activity for the species at
    //! their standard states at the current *T* and *P* of the solution.
    /*!
     * These are the standard state chemical potentials \f$ \mu^0_k(T,P)
     * \f$. The values are evaluated at the current temperature and pressure of
     * the solution
     *
     * @param mu      Output vector of chemical potentials.
     *                Length: m_kk.
     */
    virtual void getStandardChemPotentials(doublereal* mu) const {
        throw NotImplementedError("ThermoPhase::getStandardChemPotentials");
    }

    //! Get the nondimensional Enthalpy functions for the species at their
    //! standard states at the current *T* and *P* of the solution.
    /*!
     * @param hrt      Output vector of nondimensional standard state enthalpies.
     *                 Length: m_kk.
     */
    virtual void getEnthalpy_RT(doublereal* hrt) const {
        throw NotImplementedError("ThermoPhase::getEnthalpy_RT");
    }

    //! Get the array of nondimensional Entropy functions for the standard state
    //! species at the current *T* and *P* of the solution.
    /*!
     * @param sr   Output vector of nondimensional standard state entropies.
     *             Length: m_kk.
     */
    virtual void getEntropy_R(doublereal* sr) const {
        throw NotImplementedError("ThermoPhase::getEntropy_R");
    }

    //! Get the nondimensional Gibbs functions for the species in their standard
    //! states at the current *T* and *P* of the solution.
    /*!
     * @param grt  Output vector of nondimensional standard state Gibbs free
     *             energies. Length: m_kk.
     */
    virtual void getGibbs_RT(doublereal* grt) const {
        throw NotImplementedError("ThermoPhase::getGibbs_RT");
    }

    //! Get the Gibbs functions for the standard state of the species at the
    //! current *T* and *P* of the solution
    /*!
     * Units are Joules/kmol
     * @param gpure  Output vector of standard state Gibbs free energies.
     *               Length: m_kk.
     */
    virtual void getPureGibbs(doublereal* gpure) const {
        throw NotImplementedError("ThermoPhase::getPureGibbs");
    }

    //! Returns the vector of nondimensional Internal Energies of the standard
    //! state species at the current *T* and *P* of the solution
    /*!
     * @param urt  output vector of nondimensional standard state internal energies
     *             of the species. Length: m_kk.
     */
    virtual void getIntEnergy_RT(doublereal* urt) const {
        throw NotImplementedError("ThermoPhase::getIntEnergy_RT");
    }

    //! Get the nondimensional Heat Capacities at constant pressure for the
    //! species standard states at the current *T* and *P* of the
    //! solution
    /*!
     * @param cpr   Output vector of nondimensional standard state heat
     *              capacities. Length: m_kk.
     */
    virtual void getCp_R(doublereal* cpr) const {
        throw NotImplementedError("ThermoPhase::getCp_R");
    }

    //!  Get the molar volumes of the species standard states at the current
    //!  *T* and *P* of the solution.
    /*!
     * units = m^3 / kmol
     *
     * @param vol     Output vector containing the standard state volumes.
     *                Length: m_kk.
     */
    virtual void getStandardVolumes(doublereal* vol) const {
        throw NotImplementedError("ThermoPhase::getStandardVolumes");
    }

    //@}
    /// @name Thermodynamic Values for the Species Reference States
    //@{

    //! Returns the vector of nondimensional enthalpies of the reference state
    //! at the current temperature of the solution and the reference pressure
    //! for the species.
    /*!
     * @param hrt     Output vector containing the nondimensional reference
     *                state enthalpies. Length: m_kk.
     */
    virtual void getEnthalpy_RT_ref(doublereal* hrt) const {
        throw NotImplementedError("ThermoPhase::getEnthalpy_RT_ref");
    }

    //! Returns the vector of nondimensional Gibbs Free Energies of the
    //! reference state at the current temperature of the solution and the
    //! reference pressure for the species.
    /*!
     * @param grt     Output vector containing the nondimensional reference state
     *                Gibbs Free energies.  Length: m_kk.
     */
    virtual void getGibbs_RT_ref(doublereal* grt) const {
        throw NotImplementedError("ThermoPhase::getGibbs_RT_ref");
    }

    //! Returns the vector of the Gibbs function of the reference state at the
    //! current temperature of the solution and the reference pressure for the
    //! species.
    /*!
     * @param g       Output vector containing the reference state
     *                Gibbs Free energies. Length: m_kk. Units: J/kmol.
     */
    virtual void getGibbs_ref(doublereal* g) const {
        throw NotImplementedError("ThermoPhase::getGibbs_ref");
    }

    //! Returns the vector of nondimensional entropies of the reference state at
    //! the current temperature of the solution and the reference pressure for
    //! each species.
    /*!
     * @param er      Output vector containing the nondimensional reference
     *                state entropies. Length: m_kk.
     */
    virtual void getEntropy_R_ref(doublereal* er) const {
        throw NotImplementedError("ThermoPhase::getEntropy_R_ref");
    }

    //! Returns the vector of nondimensional internal Energies of the reference
    //! state at the current temperature of the solution and the reference
    //! pressure for each species.
    /*!
     * @param urt    Output vector of nondimensional reference state internal
     *               energies of the species. Length: m_kk
     */
    virtual void getIntEnergy_RT_ref(doublereal* urt) const {
        throw NotImplementedError("ThermoPhase::getIntEnergy_RT_ref");
    }

    //! Returns the vector of nondimensional constant pressure heat capacities
    //! of the reference state at the current temperature of the solution and
    //! reference pressure for each species.
    /*!
     * @param cprt   Output vector of nondimensional reference state
     *               heat capacities at constant pressure for the species.
     *               Length: m_kk
     */
    virtual void getCp_R_ref(doublereal* cprt) const {
        throw NotImplementedError("ThermoPhase::getCp_R_ref");
    }

    //! Get the molar volumes of the species reference states at the current
    //! *T* and *P_ref* of the solution.
    /*!
     * units = m^3 / kmol
     *
     * @param vol     Output vector containing the standard state volumes.
     *                Length: m_kk.
     */
    virtual void getStandardVolumes_ref(doublereal* vol) const {
        throw NotImplementedError("ThermoPhase::getStandardVolumes_ref");
    }

    // The methods below are not virtual, and should not be overloaded.

    //@}
    //! @name Specific Properties
    //@{

    //! Specific enthalpy. Units: J/kg.
    doublereal enthalpy_mass() const {
        return enthalpy_mole()/meanMolecularWeight();
    }

    //! Specific internal energy. Units: J/kg.
    doublereal intEnergy_mass() const {
        return intEnergy_mole()/meanMolecularWeight();
    }

    //! Specific entropy. Units: J/kg/K.
    doublereal entropy_mass() const {
        return entropy_mole()/meanMolecularWeight();
    }

    //! Specific Gibbs function. Units: J/kg.
    doublereal gibbs_mass() const {
        return gibbs_mole()/meanMolecularWeight();
    }

    //! Specific heat at constant pressure. Units: J/kg/K.
    doublereal cp_mass() const {
        return cp_mole()/meanMolecularWeight();
    }

    //! Specific heat at constant volume. Units: J/kg/K.
    doublereal cv_mass() const {
        return cv_mole()/meanMolecularWeight();
    }
    //@}

    //! Return the Gas Constant multiplied by the current temperature
    /*!
     *  The units are Joules kmol-1
     */
    doublereal RT() const {
        return temperature() * GasConstant;
    }

    /**
     * @name Setting the State
     *
     * These methods set all or part of the thermodynamic state.
     * @{
     */

    //! Set the temperature (K), pressure (Pa), and mole fractions.
    /*!
     * Note, the mole fractions are set first before the pressure is set.
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param t    Temperature (K)
     * @param p    Pressure (Pa)
     * @param x    Vector of mole fractions.
     *             Length is equal to m_kk.
     */
    virtual void setState_TPX(doublereal t, doublereal p, const doublereal* x);

    //! Set the temperature (K), pressure (Pa), and mole fractions.
    /*!
     * Note, the mole fractions are set first before the pressure is set.
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param t    Temperature (K)
     * @param p    Pressure (Pa)
     * @param x    Composition map of mole fractions. Species not in
     *             the composition map are assumed to have zero mole fraction
     */
    virtual void setState_TPX(doublereal t, doublereal p, const compositionMap& x);

    //! Set the temperature (K), pressure (Pa), and mole fractions.
    /*!
     * Note, the mole fractions are set first before the pressure is set.
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param t    Temperature (K)
     * @param p    Pressure (Pa)
     * @param x    String containing a composition map of the mole fractions.
     *             Species not in the composition map are assumed to have zero
     *             mole fraction
     */
    virtual void setState_TPX(doublereal t, doublereal p, const std::string& x);

    //! Set the internally stored temperature (K), pressure (Pa), and mass
    //! fractions of the phase.
    /*!
     * Note, the mass fractions are set first before the pressure is set.
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param t    Temperature (K)
     * @param p    Pressure (Pa)
     * @param y    Vector of mass fractions.
     *             Length is equal to m_kk.
     */
    virtual void setState_TPY(doublereal t, doublereal p, const doublereal* y);

    //! Set the internally stored temperature (K), pressure (Pa), and mass
    //! fractions of the phase
    /*!
     * Note, the mass fractions are set first before the pressure is set.
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param t    Temperature (K)
     * @param p    Pressure (Pa)
     * @param y    Composition map of mass fractions. Species not in
     *             the composition map are assumed to have zero mass fraction
     */
    virtual void setState_TPY(doublereal t, doublereal p, const compositionMap& y);

    //! Set the internally stored temperature (K), pressure (Pa), and mass
    //! fractions of the phase
    /*!
     * Note, the mass fractions are set first before the pressure is set.
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param t    Temperature (K)
     * @param p    Pressure (Pa)
     * @param y    String containing a composition map of the mass fractions.
     *             Species not in the composition map are assumed to have zero
     *             mass fraction
     */
    virtual void setState_TPY(doublereal t, doublereal p, const std::string& y);

    //! Set the temperature (K) and pressure (Pa)
    /*!
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param t    Temperature (K)
     * @param p    Pressure (Pa)
     */
    virtual void setState_TP(doublereal t, doublereal p);

    //! Set the pressure (Pa) and mole fractions.
    /*!
     * Note, the mole fractions are set first before the pressure is set.
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param p    Pressure (Pa)
     * @param x    Vector of mole fractions.
     *             Length is equal to m_kk.
     */
    virtual void setState_PX(doublereal p, doublereal* x);

    //! Set the internally stored pressure (Pa) and mass fractions.
    /*!
     * Note, the temperature is held constant during this operation. Note, the
     * mass fractions are set first before the pressure is set. Setting the
     * pressure may involve the solution of a nonlinear equation.
     *
     * @param p    Pressure (Pa)
     * @param y    Vector of mass fractions.
     *             Length is equal to m_kk.
     */
    virtual void setState_PY(doublereal p, doublereal* y);

    //! Set the internally stored specific enthalpy (J/kg) and pressure (Pa) of
    //! the phase.
    /*!
     * @param h     Specific enthalpy (J/kg)
     * @param p     Pressure (Pa)
     * @param tol   Optional parameter setting the tolerance of the calculation.
     *              Important for some applications where numerical Jacobians
     *              are being calculated.
     */
    virtual void setState_HP(double h, double p, double tol=1e-9);

    //! Set the specific internal energy (J/kg) and specific volume (m^3/kg).
    /*!
     * This function fixes the internal state of the phase so that the specific
     * internal energy and specific volume have the value of the input
     * parameters.
     *
     * @param u    specific internal energy (J/kg)
     * @param v    specific volume (m^3/kg).
     * @param tol  Optional parameter setting the tolerance of the calculation.
     *             Important for some applications where numerical Jacobians
     *             are being calculated.
     */
    virtual void setState_UV(double u, double v, double tol=1e-9);

    //! Set the specific entropy (J/kg/K) and pressure (Pa).
    /*!
     * This function fixes the internal state of the phase so that the specific
     * entropy and the pressure have the value of the input parameters.
     *
     * @param s    specific entropy (J/kg/K)
     * @param p    specific pressure (Pa).
     * @param tol   Optional parameter setting the tolerance of the calculation.
     *              Important for some applications where numerical Jacobians
     *              are being calculated.
     */
    virtual void setState_SP(double s, double p, double tol=1e-9);

    //! Set the specific entropy (J/kg/K) and specific volume (m^3/kg).
    /*!
     * This function fixes the internal state of the phase so that the specific
     * entropy and specific volume have the value of the input parameters.
     *
     * @param s    specific entropy (J/kg/K)
     * @param v    specific volume (m^3/kg).
     * @param tol   Optional parameter setting the tolerance of the calculation.
     *              Important for some applications where numerical Jacobians
     *              are being calculated.
     */
    virtual void setState_SV(double s, double v, double tol=1e-9);

    //! Set the specific entropy (J/kg/K) and temperature (K).
    /*!
     * This function fixes the internal state of the phase so that the specific
     * entropy and temperature have the value of the input parameters.
     * This base class function will throw an exception if not overridden.
     *
     * @param s    specific entropy (J/kg/K)
     * @param t    temperature (K)
     * @param tol   Optional parameter setting the tolerance of the calculation.
     *              Important for some applications where numerical Jacobians
     *              are being calculated.
     */
    virtual void setState_ST(double s, double t, double tol=1e-9) {
        throw NotImplementedError("ThermoPhase::setState_ST");
    }

    //! Set the temperature (K) and specific volume (m^3/kg).
    /*!
     * This function fixes the internal state of the phase so that the
     * temperature and specific volume have the value of the input parameters.
     * This base class function will throw an exception if not overridden.
     *
     * @param t    temperature (K)
     * @param v    specific volume (m^3/kg)
     * @param tol   Optional parameter setting the tolerance of the calculation.
     *              Important for some applications where numerical Jacobians
     *              are being calculated.
     */
    virtual void setState_TV(double t, double v, double tol=1e-9) {
        throw NotImplementedError("ThermoPhase::setState_TV");
    }

    //! Set the pressure (Pa) and specific volume (m^3/kg).
    /*!
     * This function fixes the internal state of the phase so that the
     * pressure and specific volume have the value of the input parameters.
     * This base class function will throw an exception if not overridden.
     *
     * @param p    pressure (Pa)
     * @param v    specific volume (m^3/kg)
     * @param tol   Optional parameter setting the tolerance of the calculation.
     *              Important for some applications where numerical Jacobians
     *              are being calculated.
     */
    virtual void setState_PV(double p, double v, double tol=1e-9) {
        throw NotImplementedError("ThermoPhase::setState_PV");
    }

    //! Set the specific internal energy (J/kg) and pressure (Pa).
    /*!
     * This function fixes the internal state of the phase so that the specific
     * internal energy and pressure have the value of the input parameters.
     * This base class function will throw an exception if not overridden.
     *
     * @param u    specific internal energy (J/kg)
     * @param p    pressure (Pa)
     * @param tol   Optional parameter setting the tolerance of the calculation.
     *              Important for some applications where numerical Jacobians
     *              are being calculated.
     */
    virtual void setState_UP(double u, double p, double tol=1e-9) {
        throw NotImplementedError("ThermoPhase::setState_UP");
    }

    //! Set the specific volume (m^3/kg) and the specific enthalpy (J/kg)
    /*!
     * This function fixes the internal state of the phase so that the specific
     * volume and the specific enthalpy have the value of the input parameters.
     * This base class function will throw an exception if not overridden.
     *
     * @param v    specific volume (m^3/kg)
     * @param h    specific enthalpy (J/kg)
     * @param tol   Optional parameter setting the tolerance of the calculation.
     *              Important for some applications where numerical Jacobians
     *              are being calculated.
     */
    virtual void setState_VH(double v, double h, double tol=1e-9) {
        throw NotImplementedError("ThermoPhase::setState_VH");
    }

    //! Set the temperature (K) and the specific enthalpy (J/kg)
    /*!
     * This function fixes the internal state of the phase so that the
     * temperature and specific enthalpy have the value of the input parameters.
     * This base class function will throw an exception if not overridden.
     *
     * @param t    temperature (K)
     * @param h    specific enthalpy (J/kg)
     * @param tol   Optional parameter setting the tolerance of the calculation.
     *              Important for some applications where numerical Jacobians
     *              are being calculated.
     */
    virtual void setState_TH(double t, double h, double tol=1e-9) {
        throw NotImplementedError("ThermoPhase::setState_TH");
    }

    //! Set the specific entropy (J/kg/K) and the specific enthalpy (J/kg)
    /*!
     * This function fixes the internal state of the phase so that the
     * temperature and pressure have the value of the input parameters.
     * This base class function will throw an exception if not overridden.
     *
     * @param s    specific entropy (J/kg/K)
     * @param h    specific enthalpy (J/kg)
     * @param tol   Optional parameter setting the tolerance of the calculation.
     *              Important for some applications where numerical Jacobians
     *              are being calculated.
     */
    virtual void setState_SH(double s, double h, double tol=1e-9) {
        throw NotImplementedError("ThermoPhase::setState_SH");
    }

    //! Set the density (kg/m**3) and pressure (Pa) at constant composition
    /*!
     * This method must be reimplemented in derived classes, where it may
     * involve the solution of a nonlinear equation. Within %Cantera, the
     * independent variable is the density. Therefore, this function solves for
     * the temperature that will yield the desired input pressure and density.
     * The composition is held constant during this process.
     *
     * This base class function will print an error, if not overridden.
     *
     * @param rho Density (kg/m^3)
     * @param p   Pressure (Pa)
     */
    virtual void setState_RP(doublereal rho, doublereal p) {
        throw NotImplementedError("ThermoPhase::setState_RP");
    }

    //! Set the density (kg/m**3), pressure (Pa) and mole fractions
    /*!
     * Note, the mole fractions are set first before the density and pressure
     * are set. Setting the pressure may involve the solution of a nonlinear
     * equation.
     *
     * @param rho  Density (kg/m^3)
     * @param p    Pressure (Pa)
     * @param x    Vector of mole fractions.
     *             Length is equal to m_kk.
     */
    virtual void setState_RPX(doublereal rho, doublereal p, const doublereal* x);

    //! Set the density (kg/m**3), pressure (Pa) and mole fractions
    /*!
     * Note, the mole fractions are set first before the density and pressure
     * are set. Setting the pressure may involve the solution of a nonlinear
     * equation.
     *
     * @param rho  Density (kg/m^3)
     * @param p    Pressure (Pa)
     * @param x    Composition map of mole fractions. Species not in
     *              the composition map are assumed to have zero mole fraction
     */
    virtual void setState_RPX(doublereal rho, doublereal p, const compositionMap& x);

    //! Set the density (kg/m**3), pressure (Pa) and mole fractions
    /*!
     * Note, the mole fractions are set first before the density and pressure
     * are set. Setting the pressure may involve the solution of a nonlinear
     * equation.
     *
     * @param rho  Density (kg/m^3)
     * @param p    Pressure (Pa)
     * @param x    String containing a composition map of the mole fractions.
     *              Species not in the composition map are assumed to have zero
     *              mole fraction
     */
    virtual void setState_RPX(doublereal rho, doublereal p, const std::string& x);

    //! Set the density (kg/m**3), pressure (Pa) and mass fractions
    /*!
     * Note, the mass fractions are set first before the density and pressure
     * are set. Setting the pressure may involve the solution of a nonlinear
     * equation.
     *
     * @param rho  Density (kg/m^3)
     * @param p    Pressure (Pa)
     * @param y    Vector of mole fractions.
     *              Length is equal to m_kk.
     */
    virtual void setState_RPY(doublereal rho, doublereal p, const doublereal* y);

    //! Set the density (kg/m**3), pressure (Pa) and mass fractions
    /*!
     * Note, the mass fractions are set first before the density and pressure
     * are set. Setting the pressure may involve the solution of a nonlinear
     * equation.
     *
     * @param rho Density (kg/m^3)
     * @param p   Pressure (Pa)
     * @param y   Composition map of mole fractions. Species not in
     *             the composition map are assumed to have zero mole fraction
     */
    virtual void setState_RPY(doublereal rho, doublereal p, const compositionMap& y);

    //! Set the density (kg/m**3), pressure (Pa) and mass fractions
    /*!
     * Note, the mass fractions are set first before the density and pressure
     * are set. Setting the pressure may involve the solution of a nonlinear
     * equation.
     *
     * @param rho  Density (kg/m^3)
     * @param p    Pressure (Pa)
     * @param y    String containing a composition map of the mole fractions.
     *              Species not in the composition map are assumed to have zero
     *              mole fraction
     */
    virtual void setState_RPY(doublereal rho, doublereal p, const std::string& y);

    //! Set the state using an AnyMap containing any combination of properties
    //! supported by the thermodynamic model
    /*!
     * Accepted keys are:
     * * `X` (mole fractions)
     * * `Y` (mass fractions)
     * * `T` or `temperature`
     * * `P` or `pressure` [Pa]
     * * `H` or `enthalpy` [J/kg]
     * * `U` or `internal-energy` [J/kg]
     * * `S` or `entropy` [J/kg/K]
     * * `V` or `specific-volume` [m^3/kg]
     * * `D` or `density` [kg/m^3]
     *
     * Composition can be specified as either an AnyMap of species names to
     * values or as a composition string. All other values can be given as
     * floating point values in Cantera's default units, or as strings with the
     * units specified, which will be converted using the Units class.
     *
     * If no thermodynamic property pair is given, or only one of temperature or
     * pressure is given, then 298.15 K and 101325 Pa will be used as necessary
     * to fully set the state.
     */
    virtual void setState(const AnyMap& state);

    //@}

private:
    //! Carry out work in HP and UV calculations.
    /*!
     * @param h     Specific enthalpy or internal energy (J/kg)
     * @param p     Pressure (Pa) or specific volume (m^3/kg)
     * @param tol   Optional parameter setting the tolerance of the calculation.
     *              Important for some applications where numerical Jacobians
     *              are being calculated.
     * @param doUV  True if solving for UV, false for HP.
     */
    void setState_HPorUV(doublereal h, doublereal p,
                         doublereal tol=1e-9, bool doUV = false);

    //! Carry out work in SP and SV calculations.
    /*!
     * @param s     Specific entropy (J/kg)
     * @param p     Pressure (Pa) or specific volume (m^3/kg)
     * @param tol   Optional parameter setting the tolerance of the calculation.
     *              Important for some applications where numerical Jacobians
     *              are being calculated.
     * @param doSV  True if solving for SV, false for SP.
     */
    void setState_SPorSV(double s, double p, double tol=1e-9, bool doSV = false);

    //! Helper function used by setState_HPorUV and setState_SPorSV.
    //! Sets the temperature and (if set_p is true) the pressure.
    void setState_conditional_TP(doublereal t, doublereal p, bool set_p);

public:
    /**
     * @name Chemical Equilibrium
     * Chemical equilibrium.
     * @{
     */

    //! Equilibrate a ThermoPhase object
    /*!
     *  Set this phase to chemical equilibrium by calling one of several
     *  equilibrium solvers. The XY parameter indicates what two thermodynamic
     *  quantities are to be held constant during the equilibration process.
     *
     *  @param XY      String representation of what two properties are being
     *                 held constant
     *  @param solver  Name of the solver to be used to equilibrate the phase.
     *      If solver = 'element_potential', the ChemEquil element potential
     *      solver will be used. If solver = 'vcs', the VCS solver will be used.
     *      If solver = 'gibbs', the MultiPhaseEquil solver will be used. If
     *      solver = 'auto', the solvers will be tried in order if the initial
     *      solver(s) fail.
     *  @param rtol      Relative tolerance
     *  @param max_steps Maximum number of steps to take to find the solution
     *  @param max_iter  For the 'gibbs' and 'vcs' solvers, this is the maximum
     *      number of outer temperature or pressure iterations to take when T
     *      and/or P is not held fixed.
     *  @param estimate_equil For MultiPhaseEquil solver, an integer indicating
     *      whether the solver should estimate its own initial condition. If 0,
     *      the initial mole fraction vector in the ThermoPhase object is used
     *      as the initial condition. If 1, the initial mole fraction vector is
     *      used if the element abundances are satisfied. If -1, the initial
     *      mole fraction vector is thrown out, and an estimate is formulated.
     *  @param log_level  loglevel Controls amount of diagnostic output.
     *      log_level=0 suppresses diagnostics, and increasingly-verbose
     *      messages are written as loglevel increases.
     *
     * @ingroup equilfunctions
     */
    void equilibrate(const std::string& XY, const std::string& solver="auto",
                     double rtol=1e-9, int max_steps=50000, int max_iter=100,
                     int estimate_equil=0, int log_level=0);

    //!This method is used by the ChemEquil equilibrium solver.
    /*!
     * It sets the state such that the chemical potentials satisfy
     * \f[ \frac{\mu_k}{\hat R T} = \sum_m A_{k,m}
     * \left(\frac{\lambda_m} {\hat R T}\right) \f] where
     * \f$ \lambda_m \f$ is the element potential of element m. The
     * temperature is unchanged.  Any phase (ideal or not) that
     * implements this method can be equilibrated by ChemEquil.
     *
     * @param lambda_RT Input vector of dimensionless element potentials
     *                  The length is equal to nElements().
     */
    virtual void setToEquilState(const doublereal* lambda_RT) {
        throw NotImplementedError("ThermoPhase::setToEquilState");
    }

    //! Indicates whether this phase type can be used with class MultiPhase for
    //! equilibrium calculations. Returns `false` for special phase types which
    //! already represent multi-phase mixtures, namely PureFluidPhase.
    virtual bool compatibleWithMultiPhase() const {
        return true;
    }

    //@}
    /// @name Critical State Properties.
    /// These methods are only implemented by subclasses that implement
    /// liquid-vapor equations of state.
    //@{

    /// Critical temperature (K).
    virtual doublereal critTemperature() const {
        throw NotImplementedError("ThermoPhase::critTemperature");
    }

    /// Critical pressure (Pa).
    virtual doublereal critPressure() const {
        throw NotImplementedError("ThermoPhase::critPressure");
    }

    /// Critical volume (m3/kmol).
    virtual doublereal critVolume() const {
        throw NotImplementedError("ThermoPhase::critVolume");
    }

    /// Critical compressibility (unitless).
    virtual doublereal critCompressibility() const {
        throw NotImplementedError("ThermoPhase::critCompressibility");
    }

    /// Critical density (kg/m3).
    virtual doublereal critDensity() const {
        throw NotImplementedError("ThermoPhase::critDensity");
    }

    //@}

    /** @name Saturation Properties.
     *
     *  These methods are only implemented by subclasses that implement full
     *  liquid-vapor equations of state.
     */
    //@{

    //! Return the saturation temperature given the pressure
    /*!
     * @param p Pressure (Pa)
     */
    virtual doublereal satTemperature(doublereal p) const {
        throw NotImplementedError("ThermoPhase::satTemperature");
    }

    //! Return the saturation pressure given the temperature
    /*!
     * @param t Temperature (Kelvin)
     */
    virtual doublereal satPressure(doublereal t) {
        throw NotImplementedError("ThermoPhase::satPressure");
    }

    //! Return the fraction of vapor at the current conditions
    virtual doublereal vaporFraction() const {
        throw NotImplementedError("ThermoPhase::vaporFraction");
    }

    //! Set the state to a saturated system at a particular temperature
    /*!
     * @param t  Temperature (kelvin)
     * @param x  Fraction of vapor
     */
    virtual void setState_Tsat(doublereal t, doublereal x) {
        throw NotImplementedError("ThermoPhase::setState_Tsat");
    }

    //! Set the state to a saturated system at a particular pressure
    /*!
     * @param p  Pressure (Pa)
     * @param x  Fraction of vapor
     */
    virtual void setState_Psat(doublereal p, doublereal x) {
        throw NotImplementedError("ThermoPhase::setState_Psat");
    }

    //! Set the temperature, pressure, and vapor fraction (quality).
    /*!
     * An exception is thrown if the thermodynamic state is not consistent.
     *
     * For temperatures below the critical temperature, if the vapor fraction is
     * not 0 or 1, the pressure and temperature must fall on the saturation
     * line.
     *
     * Above the critical temperature, the vapor fraction must be 1 if the
     * pressure is less than the critical pressure. Above the critical pressure,
     * the vapor fraction is not defined, and its value is ignored.
     *
     * @param T    Temperature (K)
     * @param P    Pressure (Pa)
     * @param Q    vapor fraction
     */
    void setState_TPQ(double T, double P, double Q);

    //@}

    //! @name Initialization Methods - For Internal Use (ThermoPhase)
    /*!
     * The following methods are used in the process of constructing
     * the phase and setting its parameters from a specification in an
     * input file. They are not normally used in application programs.
     * To see how they are used, see importPhase().
     */
    //@{

    virtual bool addSpecies(shared_ptr<Species> spec);

    virtual void modifySpecies(size_t k, shared_ptr<Species> spec);

    //! Store a reference pointer to the XML tree containing the species data
    //! for this phase.
    /*!
     *   This is used to access data needed to construct transport manager later.
     *   @internal
     *
     * @param k      Species index
     * @param data   Pointer to the XML_Node data containing
     *               information about the species in the phase.
     *
     * @deprecated The XML input format is deprecated and will be removed in
     *     Cantera 3.0.
     */
    void saveSpeciesData(const size_t k, const XML_Node* const data);

    //!  Return a pointer to the vector of XML nodes containing the species
    //!  data for this phase.
    //!
    //! @deprecated The XML input format is deprecated and will be removed in
    //!     Cantera 3.0.
    const std::vector<const XML_Node*> & speciesData() const;

    //! Return a changeable reference to the calculation manager for species
    //! reference-state thermodynamic properties
    /*!
     * @param k   Species id. The default is -1, meaning return the default
     *
     * @internal
     */
    virtual MultiSpeciesThermo& speciesThermo(int k = -1);
    
    virtual const MultiSpeciesThermo& speciesThermo(int k = -1) const;

    /**
     * @internal
     * Initialize a ThermoPhase object using an input file.
     *
     * Used to implement constructors for derived classes which take a
     * file name and phase name as arguments.
     *
     * @param inputFile Input file containing the description of the phase
     * @param id  Optional parameter identifying the name of the phase. If
     *            blank, the first phase definition encountered will be used.
     */
    virtual void initThermoFile(const std::string& inputFile,
                                const std::string& id);

    //!Import and initialize a ThermoPhase object using an XML tree.
    /*!
     * @internal
     *
     * Here we read extra information about the XML description of a phase.
     * Regular information about elements and species and their reference state
     * thermodynamic information have already been read at this point. For
     * example, we do not need to call this function for ideal gas equations of
     * state. This function is called from importPhase() after the elements and
     * the species are initialized with default ideal solution level data.
     *
     * The default implementation in ThermoPhase calls the virtual function
     * initThermo() and then sets the "state" of the phase by looking for an XML
     * element named "state", and then interpreting its contents by calling the
     * virtual function setStateFromXML().
     *
     * @param phaseNode This object must be the phase node of a complete XML
     *     tree description of the phase, including all of the species data. In
     *     other words while "phase" must point to an XML phase object, it must
     *     have sibling nodes "speciesData" that describe the species in the
     *     phase.
     * @param id   ID of the phase. If nonnull, a check is done to see if
     *             phaseNode is pointing to the phase with the correct id.
     *
     * @deprecated The XML input format is deprecated and will be removed in
     *     Cantera 3.0.
     */
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

    //! Initialize the ThermoPhase object after all species have been set up
    /*!
     * @internal Initialize.
     *
     * This method is provided to allow subclasses to perform any initialization
     * required after all species have been added. For example, it might be used
     * to resize internal work arrays that must have an entry for each species.
     * The base class implementation does nothing, and subclasses that do not
     * require initialization do not need to overload this method. Derived
     * classes which do override this function should call their parent class's
     * implementation of this function as their last action.
     *
     * When importing a CTML phase description, this method is called from
     * initThermoXML(), which is called from importPhase(), just prior to
     * returning from function importPhase().
     *
     * When importing from an AnyMap phase description (or from a YAML file),
     * this method is responsible for setting model parameters from the data
     * stored in #m_input.
     */
    virtual void initThermo();

    //! Set the equation of state parameters
    /*!
     * @internal The number and meaning of these depends on the subclass.
     *
     * @param n number of parameters
     * @param c array of \a n coefficients
     */
    virtual void setParameters(int n, doublereal* const c) {
    }

    //! Get the equation of state parameters in a vector
    /*!
     * @internal The number and meaning of these depends on the subclass.
     *
     * @param n number of parameters
     * @param c array of \a n coefficients
     */
    virtual void getParameters(int& n, doublereal* const c) const {
    }

    //! Set equation of state parameters from an AnyMap phase description.
    //! Phases that need additional parameters from the root node should
    //! override this method.
    virtual void setParameters(const AnyMap& phaseNode,
                               const AnyMap& rootNode=AnyMap());

    //! Access input data associated with the phase description
    const AnyMap& input() const;
    AnyMap& input();

    //! Set equation of state parameter values from XML entries.
    /*!
     * This method is called by function importPhase() when processing a phase
     * definition in an input file. It should be overloaded in subclasses to set
     * any parameters that are specific to that particular phase model. Note,
     * this method is called before the phase is initialized with elements
     * and/or species.
     *
     * @param eosdata An XML_Node object corresponding to
     *                the "thermo" entry for this phase in the input file.
     *
     * @deprecated The XML input format is deprecated and will be removed in
     *     Cantera 3.0.
     */
    virtual void setParametersFromXML(const XML_Node& eosdata) {}

    //! Set the initial state of the phase to the conditions specified in the
    //! state XML element.
    /*!
     * This method sets the temperature, pressure, and mole fraction vector to a
     * set default value.
     *
     * @param state AN XML_Node object corresponding to the "state" entry for
     *              this phase in the input file.
     *
     * @deprecated The XML input format is deprecated and will be removed in
     *     Cantera 3.0.
     */
    virtual void setStateFromXML(const XML_Node& state);

    virtual void invalidateCache();

    //! @}
    //! @name  Derivatives of Thermodynamic Variables needed for Applications
    //! @{

    //! Get the change in activity coefficients wrt changes in state (temp, mole
    //! fraction, etc) along a line in parameter space or along a line in
    //! physical space
    /*!
     * @param dTds           Input of temperature change along the path
     * @param dXds           Input vector of changes in mole fraction along the
     *                       path. length = m_kk Along the path length it must
     *                       be the case that the mole fractions sum to one.
     * @param dlnActCoeffds  Output vector of the directional derivatives of the
     *                       log Activity Coefficients along the path. length =
     *                       m_kk units are 1/units(s). if s is a physical
     *                       coordinate then the units are 1/m.
     */
    virtual void getdlnActCoeffds(const doublereal dTds, const doublereal* const dXds,
                                  doublereal* dlnActCoeffds) const {
        throw NotImplementedError("ThermoPhase::getdlnActCoeffds");
    }

    //! Get the array of ln mole fraction derivatives of the log activity
    //! coefficients - diagonal component only
    /*!
     * For ideal mixtures (unity activity coefficients), this can return zero.
     * Implementations should take the derivative of the logarithm of the
     * activity coefficient with respect to the logarithm of the mole fraction
     * variable that represents the standard state. This quantity is to be used
     * in conjunction with derivatives of that mole fraction variable when the
     * derivative of the chemical potential is taken.
     *
     * units = dimensionless
     *
     * @param dlnActCoeffdlnX_diag    Output vector of derivatives of the log
     *     Activity Coefficients wrt the mole fractions. length = m_kk
     */
    virtual void getdlnActCoeffdlnX_diag(doublereal* dlnActCoeffdlnX_diag) const {
        throw NotImplementedError("ThermoPhase::getdlnActCoeffdlnX_diag");
    }

    //! Get the array of log species mole number derivatives of the log activity
    //! coefficients
    /*!
     * For ideal mixtures  (unity activity coefficients), this can return zero.
     * Implementations should take the derivative of the logarithm of the
     * activity coefficient with respect to the logarithm of the concentration-
     * like variable (i.e. moles) that represents the standard state. This
     * quantity is to be used in conjunction with derivatives of that species
     * mole number variable when the derivative of the chemical potential is
     * taken.
     *
     * units = dimensionless
     *
     * @param dlnActCoeffdlnN_diag    Output vector of derivatives of the
     *                                log Activity Coefficients. length = m_kk
     */
    virtual void getdlnActCoeffdlnN_diag(doublereal* dlnActCoeffdlnN_diag) const {
        throw NotImplementedError("ThermoPhase::getdlnActCoeffdlnN_diag");
    }

    //! Get the array of derivatives of the log activity coefficients with
    //! respect to the log of the species mole numbers
    /*!
     * Implementations should take the derivative of the logarithm of the
     * activity coefficient with respect to a species log mole number (with all
     * other species mole numbers held constant). The default treatment in the
     * ThermoPhase object is to set this vector to zero.
     *
     *  units = 1 / kmol
     *
     * dlnActCoeffdlnN[ ld * k  + m] will contain the derivative of log
     * act_coeff for the *m*-th species with respect to the number of moles of
     * the *k*-th species.
     *
     * \f[
     *     \frac{d \ln(\gamma_m) }{d \ln( n_k ) }\Bigg|_{n_i}
     * \f]
     *
     * @param ld                 Number of rows in the matrix
     * @param dlnActCoeffdlnN    Output vector of derivatives of the
     *                           log Activity Coefficients. length = m_kk * m_kk
     */
    virtual void getdlnActCoeffdlnN(const size_t ld, doublereal* const dlnActCoeffdlnN);

    virtual void getdlnActCoeffdlnN_numderiv(const size_t ld, doublereal* const dlnActCoeffdlnN);

    //! @}
    //! @name Printing
    //! @{

    //! returns a summary of the state of the phase as a string
    /*!
     * @param show_thermo If true, extra information is printed out
     *                    about the thermodynamic state of the system.
     * @param threshold   Show information about species with mole fractions
     *                    greater than *threshold*.
     */
    virtual std::string report(bool show_thermo=true,
                               doublereal threshold=-1e-14) const;

    //! returns a summary of the state of the phase to a comma separated file.
    /*!
     * To customize the data included in the report, derived classes should
     * override the getCsvReportData method.
     *
     * @param csvFile  ofstream file to print comma separated data for the phase
     */
    virtual void reportCSV(std::ofstream& csvFile) const;

    //@}

protected:
    //! Fills `names` and `data` with the column names and species thermo
    //! properties to be included in the output of the reportCSV method.
    virtual void getCsvReportData(std::vector<std::string>& names,
                                  std::vector<vector_fp>& data) const;

    //! Pointer to the calculation manager for species reference-state
    //! thermodynamic properties
    /*!
     * This class is called when the reference-state thermodynamic properties
     * of all the species in the phase needs to be evaluated.
     */
    MultiSpeciesThermo m_spthermo;

    //! Data supplied via setParameters. When first set, this may include
    //! parameters used by different phase models when initThermo() is called.
    AnyMap m_input;

    //! Vector of pointers to the species databases.
    /*!
     * This is used to access data needed to construct the transport manager and
     * other properties later in the initialization process. We create a copy of
     * the XML_Node data read in here. Therefore, we own this data.
     *
     * @deprecated The XML input format is deprecated and will be removed in
     *     Cantera 3.0.
     */
    std::vector<const XML_Node*> m_speciesData;

    //! Stored value of the electric potential for this phase. Units are Volts.
    doublereal m_phi;

    //! Boolean indicating whether a charge neutrality condition is a necessity
    /*!
     * Note, the charge neutrality condition is not a necessity for ideal gas
     * phases. There may be a net charge in those phases, because the NASA
     * polynomials for ionized species in Ideal gases take this condition into
     * account. However, liquid phases usually require charge neutrality in
     * order for their derived thermodynamics to be valid.
     */
    bool m_chargeNeutralityNecessary;

    //! Contains the standard state convention
    int m_ssConvention;

    //! last value of the temperature processed by reference state
    mutable doublereal m_tlast;
};

//! typedef for the ThermoPhase class
typedef ThermoPhase thermo_t;

}

#endif
