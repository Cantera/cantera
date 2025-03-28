/**
 *  @file ThermoPhase.h
 * Header file for class ThermoPhase, the base class for phases with
 * thermodynamic properties, and the text for the Module thermoprops
 * (see @ref thermoprops and class @link Cantera::ThermoPhase ThermoPhase@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_THERMOPHASE_H
#define CT_THERMOPHASE_H

#include "Phase.h"
#include "MultiSpeciesThermo.h"
#include "cantera/base/Units.h"
#include "cantera/base/AnyMap.h"
#include "cantera/base/Solution.h"

namespace Cantera
{

/**
 * @defgroup thermoprops Thermodynamic Properties
 *
 * These classes are used to compute the thermodynamic properties of phases of matter.
 * The main base class for describing thermodynamic properties of phases within %Cantera
 * is called ThermoPhase. %ThermoPhase is a large class that describes the interface
 * within %Cantera to thermodynamic functions for a phase.
 *
 * ## Categorizing the Different ThermoPhase Objects
 *
 * ThermoPhase objects may be cataloged into four general bins.
 *
 * The first type are those whose underlying species have a reference state associated
 * with them. The reference state describes the thermodynamic functions for a species at
 * a single reference pressure, @f$ p_0 @f$. The thermodynamic functions are specified
 * via derived objects of the SpeciesThermoInterpType object class, and usually consist
 * of polynomials in temperature such as the NASA polynomial or the Shomate polynomial.
 * Calculators for these reference states, which manage the calculation for all of the
 * species in a phase, are all derived from the virtual base class
 * SpeciesThermoInterpType. Calculators are needed because the actual calculation of the
 * reference state thermodynamics has been shown to be relatively expensive. A great
 * deal of work has gone into devising efficient schemes for calculating the
 * thermodynamic polynomials of a set of species in a phase, in particular gas species
 * in ideal gas phases whose reference state thermodynamics is specified by NASA
 * polynomials.
 *
 * The reference state thermodynamics combined with the mixing rules and an assumption
 * about the pressure dependence yields the thermodynamic functions for the phase.
 * Expressions involving the specification of the fugacities of species would fall into
 * this category of %ThermoPhase objects. Note, however, that at this time, we do not
 * have any nontrivial examples of these types of phases. In general, the independent
 * variables that completely describe the state of the system  for this class are
 * temperature, the phase density, and @f$ N - 1 @f$ species mole or mass fractions.
 * Additionally, if the phase involves charged species, the phase electric potential is
 * an added independent variable. Examples of this first class of %ThermoPhase models,
 * which includes the IdealGasPhase object, the most commonly used object with %Cantera,
 * include:
 *
 * - IdealGasPhase
 * - StoichSubstance
 * - SurfPhase
 * - EdgePhase
 * - LatticePhase
 * - LatticeSolidPhase
 * - PureFluidPhase
 * - IdealSolidSolnPhase
 * - VPStandardStateTP
 *
 * The second class of objects are all derivatives of the VPStandardStateTP class listed
 * above. These classes assume that there exists a standard state for each species in
 * the phase, where the thermodynamic functions are specified as a function of
 * temperature and pressure. Standard state objects for each species are all derived
 * from the PDSS virtual base class. In turn, these standard states may employ reference
 * state calculation to aid in their calculations. However, there are some PDSS objects
 * which do not employ reference state calculations. An example of this is real equation
 * of state for liquid water used within the calculation of brine thermodynamics. In
 * general, the independent variables that completely describe the state of the system
 * for this class are temperature, the phase pressure, and @f$ N - 1 @f$ species mole or
 * mass fractions or molalities. The standard state thermodynamics combined with the
 * mixing rules yields the thermodynamic functions for the phase. Mixing rules are given
 * in terms of specifying the molar-base activity coefficients or activities. Lists of
 * phases which belong to this group are given below
 *
 * - IdealSolnGasVPSS
 * - MolalityVPSSTP
 *
 * Note, the ideal gas and ideal solution approximations are lumped together in the
 * class IdealSolnGasVPSS, because at this level they look alike having the same mixing
 * rules with respect to the specification of the excess thermodynamic properties.
 *
 * The third class of objects are all derivatives of the MolalityVPSSTP object. They
 * assume that the standard states are temperature and pressure dependent but they also
 * assume that the standard states are molality-based. In other words, they assume that
 * the standard state of the solute species are in a pseudo state of 1 molality but at
 * infinite dilution. A solvent must be specified in these calculations, defined as the
 * first species in the phase, and its standard state is the pure solvent state. Phases
 * which belong to this group include:
 *
 * - DebyeHuckel
 * - IdealMolalSoln
 * - HMWSoln
 *
 * The fourth class of %ThermoPhase objects are stoichiometric phases. Stoichiometric
 * phases are phases which consist of one and only one species. The class
 * SingleSpeciesTP is the base class for these substances. Within the class, the general
 * %ThermoPhase interface is dumbed down so that phases consisting of one species may be
 * succinctly described. These phases may have PDSS classes or SpeciesThermoInterpType
 * calculators associated with them. In general, the independent variables that
 * completely describe the state of the system for this class are temperature and either
 * the phase density or the phase pressure. Classes in this group include:
 *
 * - StoichSubstance
 * - WaterSSTP
 *
 * ## Creating ThermoPhase objects
 *
 * Instances of subclasses of ThermoPhase should be created using the factory methods
 * newThermo(const string&, const string&), newThermo(const AnyMap&, const AnyMap&), or
 * newThermoModel(). This allows new classes to be used with the various %Cantera
 * language interfaces.
 *
 * ## Defining new thermodynamic models
 *
 * To implement a new equation of state, derive a class from ThermoPhase or a relevant
 * existing derived class and overload the virtual methods in ThermoPhase. Methods that
 * are not needed can be left unimplemented, which will cause an exception to be thrown
 * if they are called.
 */

//! @name CONSTANTS - Specification of the Molality convention
//! @{

//! Standard state uses the molar convention
const int cAC_CONVENTION_MOLAR = 0;
//! Standard state uses the molality convention
const int cAC_CONVENTION_MOLALITY = 1;

//! @}
//! @name CONSTANTS - Specification of the SS convention
//! @{

//! Standard state uses the molar convention
const int cSS_CONVENTION_TEMPERATURE = 0;
//! Standard state uses the molality convention
const int cSS_CONVENTION_VPSS = 1;
//! Standard state thermodynamics is obtained from slave ThermoPhase objects
const int cSS_CONVENTION_SLAVE = 2;
//! @}

//! Differentiate between mole fractions and mass fractions for input mixture
//! composition
enum class ThermoBasis
{
    mass,
    molar
};

//! Base class for a phase with thermodynamic properties.
/*!
 * Class ThermoPhase is the base class for the family of classes that represent
 * phases of matter of any type. It defines a common public interface, and
 * implements a few methods. Most of the methods, however, are declared virtual
 * and are meant to be overloaded in derived classes.  The standard way used
 * throughout %Cantera to compute properties of phases of matter is through
 * pointers of type `ThermoPhase*` that point to objects of subclasses of
 * ThermoPhase.
 *
 * Class ThermoPhase extends class Phase by adding methods to compute
 * thermodynamic properties in addition to the ones that are used to define the
 * state of a substance (temperature, density/pressure and composition). The
 * distinction is that the methods declared in ThermoPhase require knowing the
 * particular equation of state of the phase of interest, while those of class
 * Phase do not, since they only involve data values stored within the object.
 * These methods are then implemented by the classes derived from ThermoPhase to
 * represent a phase with a specific equation of state.
 *
 * ## Calculating and accessing thermodynamic properties
 *
 * The calculation of thermodynamic functions within %ThermoPhase is broken down roughly
 * into two or more steps. First, the standard state properties of all of the species
 * are calculated at the current temperature and at either the current pressure or at a
 * reference pressure. If the calculation is carried out at a reference pressure instead
 * of at the current pressure the calculation is called a "reference state properties"
 * calculation, just to make the distinction (even though it may be considered to be a
 * fixed-pressure standard-state calculation). The next step is to adjust the reference
 * state calculation to the current pressure. The thermodynamic functions then are
 * considered to be at the standard state of each species. Lastly the mixing
 * contributions are added to arrive at the thermodynamic functions for the solution.
 *
 * The %ThermoPhase class provides interfaces to thermodynamic properties calculated for
 * the reference state of each species, the standard state values for each species, the
 * thermodynamic functions for solution values, both on a per mole of solution basis
 * (such as ThermoPhase::enthalpy_mole()), on a per kg of solution basis, and on a
 * partial molar basis for each species (such as
 * ThermoPhase::getPartialMolarEnthalpies). At each level, functions for the enthalpy,
 * entropy, Gibbs free energy, internal energy, and volume are provided. So, 5 levels
 * (reference state, standard state, partial molar, per mole of solution, and per mass
 * of solution) and 5 functions multiplied together makes 25 possible functions. That's
 * why %ThermoPhase is such a large class.
 *
 * ## Setting the State of the phase
 *
 * Typically, the way the ThermoPhase object works is that there are a set of functions
 * that set the state of the phase via setting the internal independent variables. Then,
 * there are another set of functions that query the thermodynamic functions evaluated
 * at the current %State of the phase. Internally, most of the intermediate work
 * generally occurs at the point where the internal state of the system is set and not
 * at the time when individual thermodynamic functions are queried (though the actual
 * breakdown in work is dependent on the individual derived ThermoPhase object).
 * Therefore, for efficiency, the user should lump together queries of thermodynamic
 * functions after setting the state. Moreover, in setting the state, if the density is
 * the independent variable, the following order should be used:
 *
 * - Set the temperature
 * - Set the mole or mass fractions or set the molalities
 * - set the pressure.
 *
 * For classes which inherit from VPStandardStateTP, the above order may be used, or the
 * following order may be used. It's not important.
 *
 * - Set the temperature
 * - Set the pressure
 * - Set the mole or mass fractions or set the molalities
 *
 * See the @ref sec-thermophase-set-state "list of methods" that can be used to set
 * the complete state of ThermoPhase objects.
 *
 * ## Treatment of the phase potential and the electrochemical potential of a species
 *
 * The electrochemical potential of species k in a phase p, @f$ \zeta_k @f$, is related
 * to the chemical potential as:
 *
 * @f[
 *     \zeta_{k}(T,P) = \mu_{k}(T,P) + z_k \phi_p
 * @f]
 *
 * where @f$ \nu_k @f$ is the charge of species k, and @f$ \phi_p @f$ is the electric
 * potential of phase p.
 *
 * The potential @f$ \phi_p @f$ is tracked and internally stored within the base
 * ThermoPhase object. It constitutes a specification of the internal state of the
 * phase; it's the third state variable, the first two being temperature and density
 * (or, pressure, for incompressible equations of state). It may be set with the
 * function, setElectricPotential(), and may be queried with the function
 * electricPotential().
 *
 * Note, the overall electrochemical potential of a phase may not be changed by the
 * potential because many phases enforce charge neutrality:
 *
 * @f[
 *     0 = \sum_k z_k X_k
 * @f]
 *
 * Whether charge neutrality is necessary for a phase is also specified within the
 * ThermoPhase object, by the function call chargeNeutralityNecessary(). Note, that it
 * is not necessary for the ideal gas phase, currently. However, it is necessary for
 * liquid phases such as DebyeHuckel and HMWSoln for the proper specification of the
 * chemical potentials.
 *
 * This equation, when applied to the @f$ \zeta_k @f$ equation described above, results
 * in a zero net change in the effective Gibbs free energy of the phase. However,
 * specific charged species in the phase may increase or decrease their electrochemical
 * potentials, which will have an effect on interfacial reactions involving charged
 * species, when there is a potential drop between phases. This effect is used within
 * the InterfaceKinetics and EdgeKinetics classes.
 *
 * ## Specification of Activities and Activity Conventions
 *
 * The activity @f$ a_k @f$ and activity coefficient @f$ \gamma_k @f$ of a species in
 * solution is related to the chemical potential by
 *
 * @f[
 *   \mu_k = \mu_k^0(T,P) + \hat R T \ln a_k = \mu_k^0(T,P) + \hat R T \ln x_k \gamma_k
 * @f]
 *
 * The quantity @f$ \mu_k^0(T,P) @f$ is the standard chemical potential at unit
 * activity, which depends on the temperature and pressure, but not on the composition.
 * The activity is dimensionless. Within liquid electrolytes it's common to use a
 * molality convention, where solute species employ the molality-based activity
 * coefficients:
 *
 * @f[
 *   \mu_k =  \mu_k^\triangle(T,P) + R T \ln a_k^{\triangle} =
 *            \mu_k^\triangle(T,P) + R T \ln \frac{\gamma_k^{\triangle} m_k}{m^\triangle}
 * @f]
 *
 * And the solvent employs the convention
 * @f[
 *   \mu_o = \mu^o_o(T,P) + RT \ln a_o
 * @f]
 *
 * where @f$ a_o @f$ is often redefined in terms of the osmotic coefficient @f$ \phi
 * @f$:
 *
 * @f[
 *   \phi = \frac{- \ln a_o}{\tilde{M}_o \sum_{i \ne o} m_i}
 * @f]
 *
 * ThermoPhase classes which employ the molality based convention are all derived from
 * the MolalityVPSSTP class. See the class description for further information on its
 * capabilities.
 *
 * The activity convention used by a ThermoPhase object may be queried via the
 * activityConvention() function. A zero means molar based, while a one
 * means molality based.
 *
 * The function getActivities() returns a vector of activities. Whether these are
 * molar-based or molality-based depends on the value of activityConvention().
 *
 * The function getActivityCoefficients() always returns molar-based activity
 * coefficients regardless of the activity convention used. The function
 * MolalityVPSSTP::getMolalityActivityCoefficients() returns molality
 * based activity coefficients for those ThermoPhase objects derived
 * from the MolalityVPSSTP class. The function MolalityVPSSTP::osmoticCoefficient()
 * returns the osmotic coefficient.

 * ## Activity Concentrations: Relationship of ThermoPhase to Kinetics Expressions
 *
 * %Cantera can handle both thermodynamics and kinetics mechanisms. Reversible kinetics
 * mechanisms within %Cantera must be compatible with thermodynamics in the sense that
 * at equilibrium, or at infinite times, the concentrations of species must conform to
 * thermodynamics. This means that for every valid reversible kinetics reaction in a
 * mechanism, it must be reducible to an expression involving the ratio of the product
 * activity to the reactant activities being equal to the exponential of the
 * dimensionless standard state gibbs free energies of reaction. Irreversible kinetics
 * reactions do not have this requirement; however, their usage can yield unexpected and
 * inconsistent results in many situations.
 *
 * The actual units used in a kinetics expression depend on the context or the relative
 * field of study. For example, in gas phase kinetics, species in kinetics expressions
 * are expressed in terms of concentrations, for example, gmol cm-3. In solid phase
 * studies, however, kinetics is usually expressed in terms of unitless activities,
 * which most often equate to solid phase mole fractions. In order to accommodate
 * variability here, %Cantera has come up with the idea of activity concentrations,
 * @f$ C^a_k @f$. Activity concentrations are the expressions used directly in kinetics
 * expressions. These activity (or generalized) concentrations are used by kinetics
 * manager classes to compute the forward and reverse rates of elementary reactions.
 * Note that they may or may not have units of concentration --- they might be partial
 * pressures, mole fractions, or surface coverages, The activity concentrations for
 * species *k*, @f$ C^a_k @f$, are related to the activity for species *k*, @f$ a_k @f$,
 * via the expression:
 *
 * @f[
 *   a_k = C^a_k / C^0_k
 * @f]
 *
 * @f$ C^0_k @f$ are called standard concentrations. They serve as multiplicative
 * factors between the activities and the generalized concentrations. Standard
 * concentrations may be different for each species. They may depend on both the
 * temperature and the pressure. However, they may not depend on the composition of the
 * phase. For example, for the IdealGasPhase object the standard concentration is
 * defined as
 *
 * @f[
 *   C^0_k = \frac{P}{RT}
 * @f]
 *
 * while in many solid phase kinetics problems,
 *
 * @f[
 *   C^0_k = 1.0
 * @f]
 *
 * is employed making the units for activity concentrations in solids unitless.
 *
 * ThermoPhase member functions dealing with this concept include
 * getActivityConcentrations(), which provides a vector of the current activity
 * concentrations. The function standardConcentration() returns the standard
 * concentration of the kth species. The function logStandardConc(), returns the natural
 * log of the kth standard concentration. The function standardConcentrationUnits()
 * returns the units of the standard concentration.
 *
 * ### Equilibrium constants
 *
 * - @f$ K_a @f$ is the equilibrium constant defined in terms of the standard state
 *   Gibbs free energy values. It is by definition dimensionless.
 *
 * - @f$ K_p @f$ is the equilibrium constant defined in terms of the reference state
 *   Gibbs free energy values. It is by definition dimensionless. The pressure
 *   dependence is handled entirely on the RHS of the equilibrium expression.
 *
 * - @f$ K_c @f$ is the equilibrium constant defined in terms of the activity
 *   concentrations. The dimensions depend on the number of products and reactants.
 *
 * The kinetics manager requires the calculation of @f$ K_c @f$ for the calculation of
 * the reverse rate constant.
 *
 * @ingroup thermoprops
 */
class ThermoPhase : public Phase
{
public:
    //! Constructor. Note that ThermoPhase is meant to be used as a base class,
    //! so this constructor should not be called explicitly.
    ThermoPhase() = default;

    //! @name  Information Methods
    //! @{

    string type() const override {
        return "none";
    }

    //! Boolean indicating whether phase is ideal
    virtual bool isIdeal() const {
        return false;
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
    virtual string phaseOfMatter() const {
        return "unspecified";
    }

    /**
     * Returns the reference pressure in Pa. This function is a wrapper
     * that calls the species thermo refPressure function.
     */
    virtual double refPressure() const {
        return m_spthermo.refPressure();
    }

    //! Minimum temperature for which the thermodynamic data for the species
    //! or phase are valid.
    /*!
     * If no argument is supplied, the value returned will be the lowest
     * temperature at which the data for @e all species are valid. Otherwise,
     * the value will be only for species @e k. This function is a wrapper that
     * calls the species thermo minTemp function.
     *
     * @param k index of the species. Default is -1, which will return the max
     *          of the min value over all species.
     */
    virtual double minTemp(size_t k = npos) const {
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
    double Hf298SS(const size_t k) const {
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
    virtual void modifyOneHf298SS(const size_t k, const double Hf298New) {
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
     * temperature at which the data for @e all species are valid. Otherwise,
     * the value will be only for species @e k. This function is a wrapper that
     * calls the species thermo maxTemp function.
     *
     * @param k index of the species. Default is -1, which will return the min
     *          of the max value over all species.
     */
    virtual double maxTemp(size_t k = npos) const {
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

    //! Molar enthalpy. Units: J/kmol.
    /**
     * Returns the amount of enthalpy per mole,
     * @f[
     * \hat{h} = \sum_k X_k \hat{h}_k
     * @f]
     * @see getPartialMolarEnthalpies()
     */
    virtual double enthalpy_mole() const {
        getPartialMolarEnthalpies(m_workS.data());
        return mean_X(m_workS);
    }

    //! Molar internal energy. Units: J/kmol.
    virtual double intEnergy_mole() const {
        return enthalpy_mole() - pressure()* molarVolume();
    }

    //! Molar entropy. Units: J/kmol/K.
    /**
     * Returns the amount of entropy per mole,
     * @f[
     * \hat{s} = \sum_k X_k \hat{s}_k
     * @f]
     * @see getPartialMolarEnthalpies()
     */
    virtual double entropy_mole() const {
        getPartialMolarEntropies(m_workS.data());
        return mean_X(m_workS);
    }

    //! Molar Gibbs function. Units: J/kmol.
    /*!
     * Returns the Gibbs free energy per mole,
     * @f[
     * \hat{g} = \sum_k X_k \mu_k
     * @f]
     * @see getChemPotentials()
     */
    virtual double gibbs_mole() const {
        getChemPotentials(m_workS.data());
        return mean_X(m_workS);
    }

    //! Molar heat capacity at constant pressure. Units: J/kmol/K.
    /*!
     * @f[
     * \hat{c}_p = \sum_k X_k \hat{c}_{p,k}
     * @f]
     * @see getPartialMolarCp()
     */
    virtual double cp_mole() const {
        getPartialMolarCp(m_workS.data());
        return mean_X(m_workS);
    }

    //! Molar heat capacity at constant volume. Units: J/kmol/K.
    virtual double cv_mole() const {
        throw NotImplementedError("ThermoPhase::cv_mole",
                                  "Not implemented for phase type '{}'", type());
    }

    //! @}
    //! @name Mechanical Properties
    //! @{

    //! Returns the isothermal compressibility. Units: 1/Pa.
    /*!
     * The isothermal compressibility is defined as
     * @f[
     * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
     * @f]
     *  or
     * @f[
     * \kappa_T = \frac{1}{\rho}\left(\frac{\partial \rho}{\partial P}\right)_T
     * @f]
     */
    virtual double isothermalCompressibility() const {
        throw NotImplementedError("ThermoPhase::isothermalCompressibility",
                                  "Not implemented for phase type '{}'", type());
    }

    //! Return the volumetric thermal expansion coefficient. Units: 1/K.
    /*!
     * The thermal expansion coefficient is defined as
     * @f[
     * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
     * @f]
     */
    virtual double thermalExpansionCoeff() const {
        throw NotImplementedError("ThermoPhase::thermalExpansionCoeff",
                                  "Not implemented for phase type '{}'", type());
    }

    //! Return the speed of sound. Units: m/s.
    /*!
     * The speed of sound is defined as
     * @f[
     * c = \sqrt{\left(\frac{\partial P}{\partial\rho}\right)_s}
     * @f]
     */
    virtual double soundSpeed() const {
        throw NotImplementedError("ThermoPhase::soundSpeed",
                                  "Not implemented for phase type '{}'", type());
    }

    //! @}
    //! @name Electric Potential
    //!
    //! The phase may be at some non-zero electrical potential. These methods
    //! set or get the value of the electric potential.
    //! @{

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
    void setElectricPotential(double v) {
        m_phi = v;
        invalidateCache();
    }

    //! Returns the electric potential of this phase (V).
    /*!
     *  Units are Volts (which are Joules/coulomb)
     */
    double electricPotential() const {
        return m_phi;
    }

    //! @}
    //! @name Activities, Standard States, and Activity Concentrations
    //!
    //! The activity @f$ a_k @f$ of a species in solution is related to the
    //! chemical potential by @f[ \mu_k = \mu_k^0(T,P) + \hat R T \ln a_k. @f]
    //! The quantity @f$ \mu_k^0(T,P) @f$ is the standard chemical potential at
    //! unit activity, which depends on temperature and pressure, but not on
    //! composition. The activity is dimensionless.
    //! @{

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
     * @f$ C^a_k @f$ are defined such that @f$ a_k = C^a_k / C^0_k, @f$ where
     * @f$ C^0_k @f$ is a standard concentration defined below and @f$ a_k @f$
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
    virtual void getActivityConcentrations(double* c) const {
        throw NotImplementedError("ThermoPhase::getActivityConcentrations",
                                  "Not implemented for phase type '{}'", type());
    }

    //! Return the standard concentration for the kth species
    /*!
     * The standard concentration @f$ C^0_k @f$ used to normalize the activity
     * (that is, generalized) concentration. In many cases, this quantity will be
     * the same for all species in a phase - for example, for an ideal gas @f$
     * C^0_k = P/\hat R T @f$. For this reason, this method returns a single
     * value, instead of an array.  However, for phases in which the standard
     * concentration is species-specific (such as surface species of different
     * sizes), this method may be called with an optional parameter indicating
     * the species.
     *
     * @param k Optional parameter indicating the species. The default
     *          is to assume this refers to species 0.
     * @return
     *   Returns the standard concentration. The units are by definition
     *   dependent on the ThermoPhase and kinetics manager representation.
     */
    virtual double standardConcentration(size_t k=0) const {
        throw NotImplementedError("ThermoPhase::standardConcentration",
                                  "Not implemented for phase type '{}'", type());
    }

    //! Natural logarithm of the standard concentration of the kth species.
    /*!
     * @param k    index of the species (defaults to zero)
     */
    virtual double logStandardConc(size_t k=0) const;

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
    virtual void getActivities(double* a) const;

    //! Get the array of non-dimensional molar-based activity coefficients at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     * @param ac Output vector of activity coefficients. Length: m_kk.
     */
    virtual void getActivityCoefficients(double* ac) const {
        if (m_kk == 1) {
            ac[0] = 1.0;
        } else {
            throw NotImplementedError("ThermoPhase::getActivityCoefficients",
                                      "Not implemented for phase type '{}'", type());
        }
    }

    //! Get the array of non-dimensional molar-based ln activity coefficients at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     * @param lnac Output vector of ln activity coefficients. Length: m_kk.
     */
    virtual void getLnActivityCoefficients(double* lnac) const;

    //! @}
    //! @name  Partial Molar Properties of the Solution
    //! @{

    //! Get the species chemical potentials. Units: J/kmol.
    /*!
     * This function returns a vector of chemical potentials of the species in
     * solution at the current temperature, pressure and mole fraction of the
     * solution.
     *
     * @param mu  Output vector of species chemical
     *            potentials. Length: m_kk. Units: J/kmol
     */
    virtual void getChemPotentials(double* mu) const {
        throw NotImplementedError("ThermoPhase::getChemPotentials",
                                  "Not implemented for phase type '{}'", type());
    }

    //!  Get the species electrochemical potentials.
    /*!
     *  These are partial molar quantities.  This method adds a term @f$ F z_k
     *  \phi_p @f$ to each chemical potential. The electrochemical potential of
     *  species k in a phase p, @f$ \zeta_k @f$, is related to the chemical
     *  potential via the following equation,
     *
     *  @f[
     *            \zeta_{k}(T,P) = \mu_{k}(T,P) + F z_k \phi_p
     *  @f]
     *
     * @param mu  Output vector of species electrochemical
     *            potentials. Length: m_kk. Units: J/kmol
     */
    void getElectrochemPotentials(double* mu) const;

    //! Returns an array of partial molar enthalpies for the species
    //! in the mixture. Units (J/kmol)
    /*!
     * @param hbar    Output vector of species partial molar enthalpies.
     *                Length: m_kk. units are J/kmol.
     */
    virtual void getPartialMolarEnthalpies(double* hbar) const {
        throw NotImplementedError("ThermoPhase::getPartialMolarEnthalpies",
                                  "Not implemented for phase type '{}'", type());
    }

    //! Returns an array of partial molar entropies of the species in the
    //! solution. Units: J/kmol/K.
    /*!
     * @param sbar    Output vector of species partial molar entropies.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarEntropies(double* sbar) const {
        throw NotImplementedError("ThermoPhase::getPartialMolarEntropies",
                                  "Not implemented for phase type '{}'", type());
    }

    //! Return an array of partial molar internal energies for the
    //! species in the mixture.  Units: J/kmol.
    /*!
     * @param ubar    Output vector of species partial molar internal energies.
     *                Length = m_kk. units are J/kmol.
     */
    virtual void getPartialMolarIntEnergies(double* ubar) const {
        throw NotImplementedError("ThermoPhase::getPartialMolarIntEnergies",
                                  "Not implemented for phase type '{}'", type());
    }

    //! Return an array of partial molar heat capacities for the
    //! species in the mixture.  Units: J/kmol/K
    /*!
     * @param cpbar   Output vector of species partial molar heat
     *                capacities at constant pressure.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarCp(double* cpbar) const {
        throw NotImplementedError("ThermoPhase::getPartialMolarCp",
                                  "Not implemented for phase type '{}'", type());
    }

    //! Return an array of partial molar volumes for the
    //! species in the mixture. Units: m^3/kmol.
    /*!
     *  @param vbar   Output vector of species partial molar volumes.
     *                Length = m_kk. units are m^3/kmol.
     */
    virtual void getPartialMolarVolumes(double* vbar) const {
        throw NotImplementedError("ThermoPhase::getPartialMolarVolumes",
                                  "Not implemented for phase type '{}'", type());
    }

    //! @}
    //! @name Properties of the Standard State of the Species in the Solution
    //! @{

    //! Get the array of chemical potentials at unit activity for the species at
    //! their standard states at the current *T* and *P* of the solution.
    /*!
     * These are the standard state chemical potentials @f$ \mu^0_k(T,P)
     * @f$. The values are evaluated at the current temperature and pressure of
     * the solution
     *
     * @param mu      Output vector of chemical potentials.
     *                Length: m_kk.
     */
    virtual void getStandardChemPotentials(double* mu) const {
        throw NotImplementedError("ThermoPhase::getStandardChemPotentials",
                                  "Not implemented for phase type '{}'", type());
    }

    //! Get the nondimensional Enthalpy functions for the species at their
    //! standard states at the current *T* and *P* of the solution.
    /*!
     * @param hrt      Output vector of nondimensional standard state enthalpies.
     *                 Length: m_kk.
     */
    virtual void getEnthalpy_RT(double* hrt) const {
        throw NotImplementedError("ThermoPhase::getEnthalpy_RT",
                                  "Not implemented for phase type '{}'", type());
    }

    //! Get the array of nondimensional Entropy functions for the standard state
    //! species at the current *T* and *P* of the solution.
    /*!
     * @param sr   Output vector of nondimensional standard state entropies.
     *             Length: m_kk.
     */
    virtual void getEntropy_R(double* sr) const {
        throw NotImplementedError("ThermoPhase::getEntropy_R",
                                  "Not implemented for phase type '{}'", type());
    }

    //! Get the nondimensional Gibbs functions for the species in their standard
    //! states at the current *T* and *P* of the solution.
    /*!
     * @param grt  Output vector of nondimensional standard state Gibbs free
     *             energies. Length: m_kk.
     */
    virtual void getGibbs_RT(double* grt) const {
        throw NotImplementedError("ThermoPhase::getGibbs_RT",
                                  "Not implemented for phase type '{}'", type());
    }

    //! Get the Gibbs functions for the standard state of the species at the
    //! current *T* and *P* of the solution
    /*!
     * Units are Joules/kmol
     * @param gpure  Output vector of standard state Gibbs free energies.
     *               Length: m_kk.
     */
    virtual void getPureGibbs(double* gpure) const {
        throw NotImplementedError("ThermoPhase::getPureGibbs",
                                  "Not implemented for phase type '{}'", type());
    }

    //! Returns the vector of nondimensional Internal Energies of the standard
    //! state species at the current *T* and *P* of the solution
    /*!
     * @param urt  output vector of nondimensional standard state internal energies
     *             of the species. Length: m_kk.
     */
    virtual void getIntEnergy_RT(double* urt) const {
        throw NotImplementedError("ThermoPhase::getIntEnergy_RT",
                                  "Not implemented for phase type '{}'", type());
    }

    //! Get the nondimensional Heat Capacities at constant pressure for the
    //! species standard states at the current *T* and *P* of the
    //! solution
    /*!
     * @param cpr   Output vector of nondimensional standard state heat
     *              capacities. Length: m_kk.
     */
    virtual void getCp_R(double* cpr) const {
        throw NotImplementedError("ThermoPhase::getCp_R",
                                  "Not implemented for phase type '{}'", type());
    }

    //!  Get the molar volumes of the species standard states at the current
    //!  *T* and *P* of the solution.
    /*!
     * units = m^3 / kmol
     *
     * @param vol     Output vector containing the standard state volumes.
     *                Length: m_kk.
     */
    virtual void getStandardVolumes(double* vol) const {
        throw NotImplementedError("ThermoPhase::getStandardVolumes",
                                  "Not implemented for phase type '{}'", type());
    }

    //! @}
    //! @name Thermodynamic Values for the Species Reference States
    //! @{

    //! Returns the vector of nondimensional enthalpies of the reference state
    //! at the current temperature of the solution and the reference pressure
    //! for the species.
    /*!
     * @param hrt     Output vector containing the nondimensional reference
     *                state enthalpies. Length: m_kk.
     */
    virtual void getEnthalpy_RT_ref(double* hrt) const {
        throw NotImplementedError("ThermoPhase::getEnthalpy_RT_ref",
                                  "Not implemented for phase type '{}'", type());
    }

    //! Returns the vector of nondimensional Gibbs Free Energies of the
    //! reference state at the current temperature of the solution and the
    //! reference pressure for the species.
    /*!
     * @param grt     Output vector containing the nondimensional reference state
     *                Gibbs Free energies.  Length: m_kk.
     */
    virtual void getGibbs_RT_ref(double* grt) const {
        throw NotImplementedError("ThermoPhase::getGibbs_RT_ref",
                                  "Not implemented for phase type '{}'", type());
    }

    //! Returns the vector of the Gibbs function of the reference state at the
    //! current temperature of the solution and the reference pressure for the
    //! species.
    /*!
     * @param g       Output vector containing the reference state
     *                Gibbs Free energies. Length: m_kk. Units: J/kmol.
     */
    virtual void getGibbs_ref(double* g) const {
        throw NotImplementedError("ThermoPhase::getGibbs_ref",
                                  "Not implemented for phase type '{}'", type());
    }

    //! Returns the vector of nondimensional entropies of the reference state at
    //! the current temperature of the solution and the reference pressure for
    //! each species.
    /*!
     * @param er      Output vector containing the nondimensional reference
     *                state entropies. Length: m_kk.
     */
    virtual void getEntropy_R_ref(double* er) const {
        throw NotImplementedError("ThermoPhase::getEntropy_R_ref",
                                  "Not implemented for phase type '{}'", type());
    }

    //! Returns the vector of nondimensional internal Energies of the reference
    //! state at the current temperature of the solution and the reference
    //! pressure for each species.
    /*!
     * @param urt    Output vector of nondimensional reference state internal
     *               energies of the species. Length: m_kk
     */
    virtual void getIntEnergy_RT_ref(double* urt) const {
        throw NotImplementedError("ThermoPhase::getIntEnergy_RT_ref",
                                  "Not implemented for phase type '{}'", type());
    }

    //! Returns the vector of nondimensional constant pressure heat capacities
    //! of the reference state at the current temperature of the solution and
    //! reference pressure for each species.
    /*!
     * @param cprt   Output vector of nondimensional reference state
     *               heat capacities at constant pressure for the species.
     *               Length: m_kk
     */
    virtual void getCp_R_ref(double* cprt) const {
        throw NotImplementedError("ThermoPhase::getCp_R_ref",
                                  "Not implemented for phase type '{}'", type());
    }

    //! Get the molar volumes of the species reference states at the current
    //! *T* and *P_ref* of the solution.
    /*!
     * units = m^3 / kmol
     *
     * @param vol     Output vector containing the standard state volumes.
     *                Length: m_kk.
     */
    virtual void getStandardVolumes_ref(double* vol) const {
        throw NotImplementedError("ThermoPhase::getStandardVolumes_ref",
                                  "Not implemented for phase type '{}'", type());
    }

    // The methods below are not virtual, and should not be overloaded.

    //! @}
    //! @name Specific Properties
    //! @{

    //! Specific enthalpy. Units: J/kg.
    double enthalpy_mass() const {
        return enthalpy_mole()/meanMolecularWeight();
    }

    //! Specific internal energy. Units: J/kg.
    double intEnergy_mass() const {
        return intEnergy_mole()/meanMolecularWeight();
    }

    //! Specific entropy. Units: J/kg/K.
    double entropy_mass() const {
        return entropy_mole()/meanMolecularWeight();
    }

    //! Specific Gibbs function. Units: J/kg.
    double gibbs_mass() const {
        return gibbs_mole()/meanMolecularWeight();
    }

    //! Specific heat at constant pressure. Units: J/kg/K.
    double cp_mass() const {
        return cp_mole()/meanMolecularWeight();
    }

    //! Specific heat at constant volume. Units: J/kg/K.
    double cv_mass() const {
        return cv_mole()/meanMolecularWeight();
    }
    //! @}

    //! Return the Gas Constant multiplied by the current temperature
    /*!
     *  The units are Joules kmol-1
     */
    double RT() const {
        return temperature() * GasConstant;
    }

    //! @name Setting the State
    //! @anchor sec-thermophase-set-state
    //!
    //! These methods set all or part of the thermodynamic state.
    //! @{

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
    virtual void setState_TPX(double t, double p, const double* x);

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
    virtual void setState_TPX(double t, double p, const Composition& x);

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
    virtual void setState_TPX(double t, double p, const string& x);

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
    virtual void setState_TPY(double t, double p, const double* y);

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
    virtual void setState_TPY(double t, double p, const Composition& y);

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
    virtual void setState_TPY(double t, double p, const string& y);

    //! Set the temperature (K) and pressure (Pa)
    /*!
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param t    Temperature (K)
     * @param p    Pressure (Pa)
     */
    virtual void setState_TP(double t, double p);

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
        throw NotImplementedError("ThermoPhase::setState_ST",
                                  "Not implemented for phase type '{}'", type());
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
        throw NotImplementedError("ThermoPhase::setState_TV",
                                  "Not implemented for phase type '{}'", type());
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
        throw NotImplementedError("ThermoPhase::setState_PV",
                                  "Not implemented for phase type '{}'", type());
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
        throw NotImplementedError("ThermoPhase::setState_UP",
                                  "Not implemented for phase type '{}'", type());
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
        throw NotImplementedError("ThermoPhase::setState_VH",
                                  "Not implemented for phase type '{}'", type());
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
        throw NotImplementedError("ThermoPhase::setState_TH",
                                  "Not implemented for phase type '{}'", type());
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
        throw NotImplementedError("ThermoPhase::setState_SH",
                                  "Not implemented for phase type '{}'", type());
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
     * @since New in %Cantera 3.0.
     */
    virtual void setState_DP(double rho, double p) {
        throw NotImplementedError("ThermoPhase::setState_DP",
                                  "Not implemented for phase type '{}'", type());
    }

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

    //! @}
    //! @name Set Mixture Composition by Mixture Fraction
    //! @{

    //! Set the mixture composition according to the
    //! mixture fraction = kg fuel / (kg oxidizer + kg fuel)
    /*!
     * Fuel and oxidizer compositions are given either as
     * mole fractions or mass fractions (specified by `basis`)
     * and do not need to be normalized. Pressure and temperature are
     * kept constant. Elements C, S, H and O are considered for the oxidation.
     *
     * @param mixFrac    mixture fraction (between 0 and 1)
     * @param fuelComp   composition of the fuel
     * @param oxComp     composition of the oxidizer
     * @param basis      either ThermoPhase::molar or ThermoPhase::mass.
     *                   Fuel and oxidizer composition are interpreted
     *                   as mole or mass fractions (default: molar)
     */
    void setMixtureFraction(double mixFrac, const double* fuelComp,
                            const double* oxComp, ThermoBasis basis=ThermoBasis::molar);
    //! @copydoc ThermoPhase::setMixtureFraction
    void setMixtureFraction(double mixFrac, const string& fuelComp,
            const string& oxComp, ThermoBasis basis=ThermoBasis::molar);
    //! @copydoc ThermoPhase::setMixtureFraction
    void setMixtureFraction(double mixFrac, const Composition& fuelComp,
            const Composition& oxComp, ThermoBasis basis=ThermoBasis::molar);
    //! @}
    //! @name Compute Mixture Fraction
    //! @{

    //! Compute the mixture fraction = kg fuel / (kg oxidizer + kg fuel) for
    //! the current mixture given fuel and oxidizer compositions.
    /*!
     * Fuel and oxidizer compositions are given either as
     * mole fractions or mass fractions (specified by `basis`)
     * and do not need to be normalized.
     * The mixture fraction @f$ Z @f$ can be computed from a single element
     * @f[ Z_m = \frac{Z_{\mathrm{mass},m}-Z_{\mathrm{mass},m,\mathrm{ox}}}
     * {Z_{\mathrm{mass},\mathrm{fuel}}-Z_{\mathrm{mass},m,\mathrm{ox}}} @f] where
     * @f$ Z_{\mathrm{mass},m} @f$ is the elemental mass fraction of element m
     * in the mixture, and @f$ Z_{\mathrm{mass},m,\mathrm{ox}} @f$ and
     * @f$ Z_{\mathrm{mass},m,\mathrm{fuel}} @f$ are the elemental mass fractions
     * of the oxidizer and fuel, or from the Bilger mixture fraction,
     * which considers the elements C, S, H and O @cite bilger1979
     * @f[ Z_{\mathrm{Bilger}} = \frac{\beta-\beta_{\mathrm{ox}}}
     * {\beta_{\mathrm{fuel}}-\beta_{\mathrm{ox}}} @f]
     * with @f$ \beta = 2\frac{Z_C}{M_C}+2\frac{Z_S}{M_S}+\frac{1}{2}\frac{Z_H}{M_H}
     * -\frac{Z_O}{M_O} @f$
     * and @f$ M_m @f$ the atomic weight of element @f$ m @f$.
     *
     * @param fuelComp   composition of the fuel
     * @param oxComp     composition of the oxidizer
     * @param basis      either ThermoBasis::molar or ThermoBasis::mass.
     *                   Fuel and oxidizer composition are interpreted
     *                   as mole or mass fractions (default: molar)
     * @param element    either "Bilger" to compute the mixture fraction
     *                   in terms of the Bilger mixture fraction, or
     *                   an element name, to compute the mixture fraction
     *                   based on a single element (default: "Bilger")
     * @returns          mixture fraction (kg fuel / kg mixture)
     */
    double mixtureFraction(const double* fuelComp, const double* oxComp,
                           ThermoBasis basis=ThermoBasis::molar,
                           const string& element="Bilger") const;
    //! @copydoc ThermoPhase::mixtureFraction
    double mixtureFraction(const string& fuelComp, const string& oxComp,
                           ThermoBasis basis=ThermoBasis::molar,
                           const string& element="Bilger") const;
    //! @copydoc ThermoPhase::mixtureFraction
    double mixtureFraction(const Composition& fuelComp, const Composition& oxComp,
                           ThermoBasis basis=ThermoBasis::molar,
                           const string& element="Bilger") const;
    //! @}
    //! @name Set Mixture Composition by Equivalence Ratio
    //! @{

    //! Set the mixture composition according to the equivalence ratio.
    /*!
     * Fuel and oxidizer compositions are given either as
     * mole fractions or mass fractions (specified by `basis`)
     * and do not need to be normalized. Pressure and temperature are
     * kept constant. Elements C, S, H and O are considered for the oxidation.
     *
     * @param phi        equivalence ratio
     * @param fuelComp   composition of the fuel
     * @param oxComp     composition of the oxidizer
     * @param basis      either ThermoBasis::mole or ThermoBasis::mass.
     *                   Fuel and oxidizer composition are interpreted
     *                   as mole or mass fractions (default: molar)
     */
    void setEquivalenceRatio(double phi, const double* fuelComp, const double* oxComp,
                             ThermoBasis basis=ThermoBasis::molar);
    //! @copydoc ThermoPhase::setEquivalenceRatio
    void setEquivalenceRatio(double phi, const string& fuelComp,
            const string& oxComp, ThermoBasis basis=ThermoBasis::molar);
    //! @copydoc ThermoPhase::setEquivalenceRatio
    void setEquivalenceRatio(double phi, const Composition& fuelComp,
            const Composition& oxComp, ThermoBasis basis=ThermoBasis::molar);
    //! @}

    //! @name Compute Equivalence Ratio
    //! @{

    //! Compute the equivalence ratio for the current mixture
    //! given the compositions of fuel and oxidizer
    /*!
     * The equivalence ratio @f$ \phi @f$ is computed from
     * @f[ \phi = \frac{Z}{1-Z}\frac{1-Z_{\mathrm{st}}}{Z_{\mathrm{st}}} @f]
     * where @f$ Z @f$ is the Bilger mixture fraction @cite bilger1979 of the mixture
     * given the specified fuel and oxidizer compositions
     * @f$ Z_{\mathrm{st}} @f$ is the mixture fraction at stoichiometric
     * conditions. Fuel and oxidizer compositions are given either as
     * mole fractions or mass fractions (specified by `basis`)
     * and do not need to be normalized.
     * Elements C, S, H and O are considered for the oxidation.
     * If fuel and oxidizer composition are unknown or not specified,
     * use the version that takes no arguments.
     *
     * @param fuelComp   composition of the fuel
     * @param oxComp     composition of the oxidizer
     * @param basis      either ThermoPhase::mole or ThermoPhase::mass.
     *                   Fuel and oxidizer composition are interpreted
     *                   as mole or mass fractions (default: molar)
     * @returns          equivalence ratio
     * @see mixtureFraction for the definition of the Bilger mixture fraction
     * @see equivalenceRatio() for the computation of @f$ \phi @f$ without arguments
     */
    double equivalenceRatio(const double* fuelComp, const double* oxComp,
                            ThermoBasis basis=ThermoBasis::molar) const;
    //! @copydoc ThermoPhase::equivalenceRatio
    double equivalenceRatio(const string& fuelComp, const string& oxComp,
                            ThermoBasis basis=ThermoBasis::molar) const;
    //! @copydoc ThermoPhase::equivalenceRatio
    double equivalenceRatio(const Composition& fuelComp,
            const Composition& oxComp, ThermoBasis basis=ThermoBasis::molar) const;
    //! @}

    //! Compute the equivalence ratio for the current mixture
    //! from available oxygen and required oxygen
    /*!
     * Computes the equivalence ratio @f$ \phi @f$ from
     * @f[ \phi =
     * \frac{Z_{\mathrm{mole},C} + Z_{\mathrm{mole},S} + \frac{1}{4}Z_{\mathrm{mole},H}}
     * {\frac{1}{2}Z_{\mathrm{mole},O}} @f]
     * where @f$ Z_{\mathrm{mole},m} @f$ is the elemental mole fraction
     * of element @f$ m @f$. In this special case, the equivalence ratio
     * is independent of a fuel or oxidizer composition because it only
     * considers the locally available oxygen compared to the required oxygen
     * for complete oxidation. It is the same as assuming that the oxidizer
     * only contains O (and inert elements) and the fuel contains only
     * H, C and S (and inert elements). If either of these conditions is
     * not met, use the version of this functions which takes the fuel and
     * oxidizer compositions as input
     *
     * @returns                equivalence ratio
     * @see equivalenceRatio   compute the equivalence ratio from specific
     *                         fuel and oxidizer compositions
     */
    double equivalenceRatio() const;

    //! @name Compute Stoichiometric Air to Fuel Ratio
    //! @{

    //! Compute the stoichiometric air to fuel ratio (kg oxidizer / kg fuel)
    //! given fuel and oxidizer compositions.
    /*!
     * Fuel and oxidizer compositions are given either as
     * mole fractions or mass fractions (specified by `basis`)
     * and do not need to be normalized.
     * Elements C, S, H and O are considered for the oxidation.
     * Note that the stoichiometric air to fuel ratio @f$ \mathit{AFR}_{\mathrm{st}} @f$
     * does not depend on the current mixture composition. The current air to fuel ratio
     * can be computed from @f$ \mathit{AFR} = \mathit{AFR}_{\mathrm{st}}/\phi @f$
     * where @f$ \phi @f$ is the equivalence ratio of the current mixture
     *
     * @param fuelComp   composition of the fuel
     * @param oxComp     composition of the oxidizer
     * @param basis      either ThermoPhase::mole or ThermoPhase::mass.
     *                   Fuel and oxidizer composition are interpreted
     *                   as mole or mass fractions (default: molar)
     * @returns          Stoichiometric Air to Fuel Ratio (kg oxidizer / kg fuel)
     */
    double stoichAirFuelRatio(const double* fuelComp, const double* oxComp,
                              ThermoBasis basis=ThermoBasis::molar) const;
    //! @copydoc ThermoPhase::stoichAirFuelRatio
    double stoichAirFuelRatio(const string& fuelComp, const string& oxComp,
                              ThermoBasis basis=ThermoBasis::molar) const;
    //! @copydoc ThermoPhase::stoichAirFuelRatio
    double stoichAirFuelRatio(const Composition& fuelComp,
            const Composition& oxComp, ThermoBasis basis=ThermoBasis::molar) const;
    //! @}

    //! Return intermediate or model-specific parameters used by particular
    //! derived classes. Specific parameters are described in overidden
    //! methods of classes that derive from the base class.
    virtual AnyMap getAuxiliaryData()
    {
        return AnyMap();
    }

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
    void setState_HPorUV(double h, double p, double tol=1e-9, bool doUV = false);

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
    void setState_conditional_TP(double t, double p, bool set_p);

    //! Helper function for computing the amount of oxygen required for complete
    //! oxidation.
    /*!
     * @param y       array of (possibly non-normalized) mass fractions (length m_kk)
     * @returns       amount of required oxygen in kmol O / kg mixture
     */
    double o2Required(const double* y) const;

    //! Helper function for computing the amount of oxygen
    //! available in the current mixture.
    /*!
     * @param y       array of (possibly non-normalized) mass fractions (length m_kk)
     * @returns       amount of O in kmol O / kg mixture
     */
    double o2Present(const double* y) const;

public:
    //! @name Chemical Equilibrium
    //!
    //! Chemical equilibrium.
    //! @{

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
     * @ingroup equilGroup
     */
    void equilibrate(const string& XY, const string& solver="auto",
                     double rtol=1e-9, int max_steps=50000, int max_iter=100,
                     int estimate_equil=0, int log_level=0);

    //!This method is used by the ChemEquil equilibrium solver.
    /*!
     * It sets the state such that the chemical potentials satisfy
     * @f[ \frac{\mu_k}{\hat R T} = \sum_m A_{k,m}
     * \left(\frac{\lambda_m} {\hat R T}\right) @f] where
     * @f$ \lambda_m @f$ is the element potential of element m. The
     * temperature is unchanged.  Any phase (ideal or not) that
     * implements this method can be equilibrated by ChemEquil.
     *
     * @param mu_RT Input vector of dimensionless chemical potentials
     *                  The length is equal to nSpecies().
     */
    virtual void setToEquilState(const double* mu_RT) {
        throw NotImplementedError("ThermoPhase::setToEquilState");
    }

    //! Indicates whether this phase type can be used with class MultiPhase for
    //! equilibrium calculations. Returns `false` for special phase types which
    //! already represent multi-phase mixtures, namely PureFluidPhase.
    virtual bool compatibleWithMultiPhase() const {
        return true;
    }

    //! @}
    //! @name Critical State Properties
    //!
    //! These methods are only implemented by subclasses that implement
    //! liquid-vapor equations of state.
    //! @{

    //! Critical temperature (K).
    virtual double critTemperature() const {
        throw NotImplementedError("ThermoPhase::critTemperature");
    }

    //! Critical pressure (Pa).
    virtual double critPressure() const {
        throw NotImplementedError("ThermoPhase::critPressure");
    }

    //! Critical volume (m3/kmol).
    virtual double critVolume() const {
        throw NotImplementedError("ThermoPhase::critVolume");
    }

    //! Critical compressibility (unitless).
    virtual double critCompressibility() const {
        throw NotImplementedError("ThermoPhase::critCompressibility");
    }

    //! Critical density (kg/m3).
    virtual double critDensity() const {
        throw NotImplementedError("ThermoPhase::critDensity");
    }

    //! @}
    //! @name Saturation Properties
    //!
    //! These methods are only implemented by subclasses that implement full
    //! liquid-vapor equations of state.
    //! @{

    //! Return the saturation temperature given the pressure
    /*!
     * @param p Pressure (Pa)
     */
    virtual double satTemperature(double p) const {
        throw NotImplementedError("ThermoPhase::satTemperature");
    }

    //! Return the saturation pressure given the temperature
    /*!
     * @param t Temperature (Kelvin)
     */
    virtual double satPressure(double t) {
        throw NotImplementedError("ThermoPhase::satPressure");
    }

    //! Return the fraction of vapor at the current conditions
    virtual double vaporFraction() const {
        throw NotImplementedError("ThermoPhase::vaporFraction");
    }

    //! Set the state to a saturated system at a particular temperature
    /*!
     * @param t  Temperature (kelvin)
     * @param x  Fraction of vapor
     */
    virtual void setState_Tsat(double t, double x) {
        throw NotImplementedError("ThermoPhase::setState_Tsat");
    }

    //! Set the state to a saturated system at a particular pressure
    /*!
     * @param p  Pressure (Pa)
     * @param x  Fraction of vapor
     */
    virtual void setState_Psat(double p, double x) {
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

    //! @}
    //! @name Initialization Methods - For Internal Use (ThermoPhase)
    //!
    //! The following methods are used in the process of constructing
    //! the phase and setting its parameters from a specification in an
    //! input file. They are not normally used in application programs.
    //! To see how they are used, see importPhase().
    //! @{

    bool addSpecies(shared_ptr<Species> spec) override;

    void modifySpecies(size_t k, shared_ptr<Species> spec) override;

    //! Return a changeable reference to the calculation manager for species
    //! reference-state thermodynamic properties
    /*!
     * @param k   Species id. The default is -1, meaning return the default
     */
    virtual MultiSpeciesThermo& speciesThermo(int k = -1);

    virtual const MultiSpeciesThermo& speciesThermo(int k = -1) const;

    /**
     * Initialize a ThermoPhase object using an input file.
     *
     * Used to implement constructors for derived classes which take a
     * file name and phase name as arguments.
     *
     * @param inputFile Input file containing the description of the phase. If blank,
     *     no setup will be performed.
     * @param id  Optional parameter identifying the name of the phase. If
     *            blank, the first phase definition encountered will be used.
     */
    void initThermoFile(const string& inputFile, const string& id);

    //! Initialize the ThermoPhase object after all species have been set up
    /*!
     * This method is provided to allow subclasses to perform any initialization
     * required after all species have been added. For example, it might be used
     * to resize internal work arrays that must have an entry for each species.
     * The base class implementation does nothing, and subclasses that do not
     * require initialization do not need to overload this method. Derived
     * classes which do override this function should call their parent class's
     * implementation of this function as their last action.
     *
     * When importing from an AnyMap phase description (or from a YAML file),
     * setupPhase() adds all the species, stores the input data in #m_input, and then
     * calls this method to set model parameters from the data stored in #m_input.
     */
    virtual void initThermo();

    //! Set equation of state parameters from an AnyMap phase description.
    //! Phases that need additional parameters from the root node should
    //! override this method.
    virtual void setParameters(const AnyMap& phaseNode,
                               const AnyMap& rootNode=AnyMap());

    //! Returns the parameters of a ThermoPhase object such that an identical
    //! one could be reconstructed using the newThermo(AnyMap&) function.
    //! @param withInput  If true, include additional input data fields associated
    //!   with the phase description, such as user-defined fields from a YAML input
    //!   file, as returned by the input() method.
    AnyMap parameters(bool withInput=true) const;

    //! Get phase-specific parameters of a Species object such that an
    //! identical one could be reconstructed and added to this phase.
    /*!
     * @param name         Name of the species
     * @param speciesNode  Mapping to be populated with parameters
     */
    virtual void getSpeciesParameters(const string& name, AnyMap& speciesNode) const {}

    //! Access input data associated with the phase description
    const AnyMap& input() const;
    AnyMap& input();

    void invalidateCache() override;

    //! @}
    //! @name  Derivatives of Thermodynamic Variables needed for Applications
    //!
    //! Derivatives of the activity coefficients are needed to evaluate terms arising
    //! in multicomponent transport models for non-ideal systems. While %Cantera does
    //! not currently implement such models, these derivatives are provided by a few
    //! phase models.
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
    virtual void getdlnActCoeffds(const double dTds, const double* const dXds,
                                  double* dlnActCoeffds) const {
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
    virtual void getdlnActCoeffdlnX_diag(double* dlnActCoeffdlnX_diag) const {
        throw NotImplementedError("ThermoPhase::getdlnActCoeffdlnX_diag");
    }

    //! Get the array of log species mole number derivatives of the log activity
    //! coefficients
    /*!
     * For ideal mixtures  (unity activity coefficients), this can return zero.
     * Implementations should take the derivative of the logarithm of the
     * activity coefficient with respect to the logarithm of the concentration-
     * like variable (for example, moles) that represents the standard state. This
     * quantity is to be used in conjunction with derivatives of that species
     * mole number variable when the derivative of the chemical potential is
     * taken.
     *
     * units = dimensionless
     *
     * @param dlnActCoeffdlnN_diag    Output vector of derivatives of the
     *                                log Activity Coefficients. length = m_kk
     */
    virtual void getdlnActCoeffdlnN_diag(double* dlnActCoeffdlnN_diag) const {
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
     * @f[
     *     \frac{d \ln(\gamma_m) }{d \ln( n_k ) }\Bigg|_{n_i}
     * @f]
     *
     * When implemented, this method is used within the VCS equilibrium solver to
     * calculate the Jacobian elements, which accelerates convergence of the algorithm.
     *
     * @param ld                 Number of rows in the matrix
     * @param dlnActCoeffdlnN    Output vector of derivatives of the
     *                           log Activity Coefficients. length = m_kk * m_kk
     */
    virtual void getdlnActCoeffdlnN(const size_t ld, double* const dlnActCoeffdlnN);

    virtual void getdlnActCoeffdlnN_numderiv(const size_t ld,
                                             double* const dlnActCoeffdlnN);

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
    virtual string report(bool show_thermo=true, double threshold=-1e-14) const;

    //! @}

    //! Set the link to the Solution object that owns this ThermoPhase
    //! @param soln Weak pointer to the parent Solution object
    virtual void setSolution(std::weak_ptr<Solution> soln) {
        m_soln = soln;
    }

protected:
    //! Store the parameters of a ThermoPhase object such that an identical
    //! one could be reconstructed using the newThermo(AnyMap&) function. This
    //! does not include user-defined fields available in input().
    virtual void getParameters(AnyMap& phaseNode) const;

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

    //! Stored value of the electric potential for this phase. Units are Volts.
    double m_phi = 0.0;

    //! Boolean indicating whether a charge neutrality condition is a necessity
    /*!
     * Note, the charge neutrality condition is not a necessity for ideal gas
     * phases. There may be a net charge in those phases, because the NASA
     * polynomials for ionized species in Ideal gases take this condition into
     * account. However, liquid phases usually require charge neutrality in
     * order for their derived thermodynamics to be valid.
     */
    bool m_chargeNeutralityNecessary = false;

    //! Contains the standard state convention
    int m_ssConvention = cSS_CONVENTION_TEMPERATURE;

    //! last value of the temperature processed by reference state
    mutable double m_tlast = 0.0;

    //! reference to Solution
    std::weak_ptr<Solution> m_soln;
};

}

#endif
