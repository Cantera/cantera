/**
 *  @file IdealGasPhase.h
 *   ThermoPhase object for the ideal gas equation of
 * state - workhorse for %Cantera (see \ref thermoprops
 * and class \link Cantera::IdealGasPhase IdealGasPhase\endlink).
 */

//  Copyright 2001 California Institute of Technology
#ifndef CT_IDEALGASPHASE_H
#define CT_IDEALGASPHASE_H

#include "mix_defs.h"
#include "ThermoPhase.h"
#include "SpeciesThermo.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

//!  Class %IdealGasPhase represents low-density gases that obey the
//!  ideal gas equation of state.
/*!
 *
 * %IdealGasPhase derives from class ThermoPhase,
 * and overloads the virtual methods defined there with ones that
 * use expressions appropriate for ideal gas mixtures.
 *
 * The independent unknowns are density, mass fraction, and temperature.
 * the #setPressure() function will calculate the density consistent with
 * the current mass fraction vector and temperature and the desired pressure,
 * and then set the density.
 *
 * <HR>
 * <H2> Specification of Species Standard State Properties </H2>
 * <HR>
 *
 *  It is assumed that the reference state thermodynamics may be
 *  obtained by a pointer to a populated species thermodynamic property
 *  manager class in the base class, ThermoPhase::m_spthermo
 *  (see the base class \link Cantera#SpeciesThermo SpeciesThermo \endlink for a
 *  description of the specification of reference state species thermodynamics functions).
 *  The reference state,
 *  where the pressure is fixed at a single pressure,
 *  is a key species property calculation for the Ideal Gas Equation
 *  of state.
 *
 *  This class is optimized for speed of execution. All calls to thermodynamic functions
 *  first call internal routines (aka #enthalpy_RT_ref()) which return references
 *  the reference state thermodynamics functions. Within these internal reference
 *  state functions, the function #_updateThermo() is called, that first checks to see
 *  whether the temperature has changed. If it has, it updates the internal reference
 *  state thermo functions by calling the SpeciesThermo object.
 *
 *  Functions for the calculation of standard state properties for species
 *  at arbitrary pressure are provided in %IdealGasPhase. However, they
 *  are all derived from their reference state counterparts.
 *
 *  The standard state enthalpy is independent of pressure:
 *
 *       \f[
 *            h^o_k(T,P) = h^{ref}_k(T)
 *       \f]
 *
 *  The standard state constant-pressure heat capacity is independent of pressure:
 *
 *       \f[
 *            Cp^o_k(T,P) = Cp^{ref}_k(T)
 *       \f]
 *
 *  The standard state entropy depends in the following fashion on pressure:
 *
 *       \f[
 *            S^o_k(T,P) = S^{ref}_k(T) -  R \ln(\frac{P}{P_{ref}})
 *       \f]
 *  The standard state gibbs free energy is obtained from the enthalpy and entropy
 *  functions:
 *
 *       \f[
 *            \mu^o_k(T,P) =  h^o_k(T,P) - S^o_k(T,P) T
 *       \f]
 *
 *       \f[
 *            \mu^o_k(T,P) =  \mu^{ref}_k(T) + R T \ln( \frac{P}{P_{ref}})
 *       \f]
 *
 * where
 *       \f[
 *            \mu^{ref}_k(T) =   h^{ref}_k(T)   - T S^{ref}_k(T)
 *       \f]
 *
 *  The standard state internal energy is obtained from the enthalpy function also
 *
 *       \f[
 *            u^o_k(T,P) = h^o_k(T) - R T
 *       \f]
 *
 *  The molar volume of a species is given by the ideal gas law
 *
 *       \f[
 *            V^o_k(T,P) = \frac{R T}{P}
 *       \f]
 *
 *  where R is the molar gas constant. For a complete list of physical constants
 *  used within %Cantera, see \ref physConstants .
 *
 * <HR>
 * <H2> Specification of Solution Thermodynamic Properties </H2>
 * <HR>
 *
 * The activity of a species defined in the phase is given by the ideal gas law:
 *       \f[
 *            a_k = X_k
 *       \f]
 * where \f$ X_k \f$ is the mole fraction of species <I>k</I>.
 * The chemical potential for species <I>k</I> is equal to
 *
 *       \f[
 *            \mu_k(T,P) = \mu^o_k(T, P) + R T \log(X_k)
 *       \f]
 *
 * In terms of the reference state, the above can be rewritten
 *
 *
 *       \f[
 *            \mu_k(T,P) = \mu^{ref}_k(T, P) + R T \log(\frac{P X_k}{P_{ref}})
 *       \f]
 *
 * The partial molar entropy for species <I>k</I> is given by the following relation,
 *
 *       \f[
 *            \tilde{s}_k(T,P) = s^o_k(T,P) - R \log(X_k) = s^{ref}_k(T) - R \log(\frac{P X_k}{P_{ref}})
 *       \f]
 *
 * The partial molar enthalpy for species <I>k</I> is
 *
 *       \f[
 *            \tilde{h}_k(T,P) = h^o_k(T,P) = h^{ref}_k(T)
 *       \f]
 *
 * The partial molar Internal Energy for species <I>k</I> is
 *
 *       \f[
 *            \tilde{u}_k(T,P) = u^o_k(T,P) = u^{ref}_k(T)
 *       \f]
 *
 * The partial molar Heat Capacity for species <I>k</I> is
 *
 *       \f[
 *            \tilde{Cp}_k(T,P) = Cp^o_k(T,P) = Cp^{ref}_k(T)
 *       \f]
 *
 *
 * <HR>
 * <H2> %Application within %Kinetics Managers </H2>
 * <HR>
 *
 *   \f$ C^a_k\f$ are defined such that \f$ a_k = C^a_k /
 *   C^s_k, \f$ where \f$ C^s_k \f$ is a standard concentration
 *   defined below and \f$ a_k \f$ are activities used in the
 *   thermodynamic functions.  These activity (or generalized)
 *   concentrations are used
 *   by kinetics manager classes to compute the forward and
 *   reverse rates of elementary reactions.
 *   The activity concentration,\f$  C^a_k \f$,is given by the following expression.
 *
 *       \f[
 *            C^a_k = C^s_k  X_k  = \frac{P}{R T} X_k
 *       \f]
 *
 * The standard concentration for species <I>k</I> is independent of <I>k</I> and equal to
 *
 *        \f[
 *            C^s_k =  C^s = \frac{P}{R T}
 *        \f]
 *
 * For example, a bulk-phase binary gas reaction between species j and k, producing
 * a new gas species l would have the
 * following equation for its rate of progress variable, \f$ R^1 \f$, which has
 * units of kmol m-3 s-1.
 *
 *   \f[
 *    R^1 = k^1 C_j^a C_k^a =  k^1 (C^s a_j) (C^s a_k)
 *   \f]
 *  where
 *   \f[
 *      C_j^a = C^s a_j \quad \mbox{and} \quad C_k^a = C^s a_k
 *   \f]
 *
 *  \f$ C_j^a \f$ is the activity concentration of species j, and
 *  \f$ C_k^a \f$ is the activity concentration of species k. \f$ C^s \f$
 *  is the standard concentration. \f$ a_j \f$ is
 *  the activity of species j which is equal to the mole fraction of j.
 *
 *  The reverse rate constant can then be obtained from the law of microscopic reversibility
 *  and the equilibrium expression for the system.
 *
 *   \f[
 *         \frac{a_j a_k}{ a_l} = K_a^{o,1} = \exp(\frac{\mu^o_l - \mu^o_j - \mu^o_k}{R T} )
 *   \f]
 *
 *  \f$  K_a^{o,1} \f$ is the dimensionless form of the equilibrium constant, associated with
 *  the pressure dependent standard states \f$ \mu^o_l(T,P) \f$ and their associated activities,
 *  \f$ a_l \f$, repeated here:
 *
 *       \f[
 *            \mu_l(T,P) = \mu^o_l(T, P) + R T \log(a_l)
 *       \f]
 *
 *  We can switch over to expressing the equilibrium constant in terms of the reference
 *  state chemical potentials
 *
 *   \f[
 *       K_a^{o,1} = \exp(\frac{\mu^{ref}_l - \mu^{ref}_j - \mu^{ref}_k}{R T} ) * \frac{P_{ref}}{P}
 *   \f]
 *
 *   The concentration equilibrium constant, \f$ K_c \f$, may be obtained by changing over
 *   to activity concentrations. When this is done:
 *
 *   \f[
 *         \frac{C^a_j C^a_k}{ C^a_l} = C^o K_a^{o,1} = K_c^1 =
 *             \exp(\frac{\mu^{ref}_l - \mu^{ref}_j - \mu^{ref}_k}{R T} ) * \frac{P_{ref}}{RT}
 *   \f]
 *
 *    %Kinetics managers will calculate the concentration equilibrium constant, \f$ K_c \f$,
 *    using the second and third part of the above expression as a definition for the concentration
 *    equilibrium constant.
 *
 *    For completeness, the pressure equilibrium constant may be obtained as well
 *
 *   \f[
 *         \frac{P_j P_k}{ P_l P_{ref}} = K_p^1 =
               \exp\left(\frac{\mu^{ref}_l - \mu^{ref}_j - \mu^{ref}_k}{R T} \right)
 *   \f]
 *
 *   \f$ K_p \f$ is the simplest form of the equilibrium constant for ideal gases. However, it isn't
 *   necessarily the simplest form of the equilibrium constant for other types of phases; \f$ K_c \f$ is
 *   used instead because it is completely general.
 *
 *   The reverse rate of progress may be written down as
 *   \f[
 *    R^{-1} = k^{-1} C_l^a =  k^{-1} (C^o a_l)
 *   \f]
 *
 *  where we can use the concept of microscopic reversibility to
 *  write the reverse rate constant in terms of the
 *  forward rate constant and the concentration equilibrium
 *   constant, \f$ K_c \f$.
 *
 *    \f[
 *       k^{-1} =  k^1 K^1_c
 *    \f]
 *
 *  \f$k^{-1} \f$ has units of s-1.
 *
 * <HR>
 * <H2> Instantiation of the Class </H2>
 * <HR>
 *
 * The constructor for this phase is located in the default ThermoFactory
 * for %Cantera. A new %IdealGasPhase may be created by the following code
 * snippet:
 *
 * @code
 *    XML_Node *xc = get_XML_File("silane.xml");
 *    XML_Node * const xs = xc->findNameID("phase", "silane");
 *    ThermoPhase *silane_tp = newPhase(*xs);
 *    IdealGasPhase *silaneGas = dynamic_cast <IdealGasPhase *>(silane_tp);
 * @endcode
 *
 * or by the following constructor:
 *
 * @code
 *    XML_Node *xc = get_XML_File("silane.xml");
 *    XML_Node * const xs = xc->findNameID("phase", "silane");
 *    IdealGasPhase *silaneGas = new IdealGasPhase(*xs);
 * @endcode
 *
 * <HR>
 * <H2> XML Example </H2>
 * <HR>
 *   An example of an XML Element named phase setting up a IdealGasPhase
 *   object named silane is given below.
 *
 * @code
 * <!--     phase silane      -->
 * <phase dim="3" id="silane">
 *   <elementArray datasrc="elements.xml"> Si  H  He </elementArray>
 *   <speciesArray datasrc="#species_data">
 *     H2  H  HE  SIH4  SI  SIH  SIH2  SIH3  H3SISIH  SI2H6
 *     H2SISIH2  SI3H8  SI2  SI3
 *   </speciesArray>
 *   <reactionArray datasrc="#reaction_data"/>
 *   <thermo model="IdealGas"/>
 *   <kinetics model="GasKinetics"/>
 *   <transport model="None"/>
 * </phase>
 * @endcode
 *
 *   The model attribute "IdealGas" of the thermo XML element identifies the phase as
 *   being of the type handled by the IdealGasPhase object.
 *
 *    @ingroup thermoprops
 *
 */
class IdealGasPhase: public ThermoPhase
{
public:
    //! Default empty Constructor
    IdealGasPhase();

    //! Construct and initialize an IdealGasPhase ThermoPhase object
    //! directly from an ASCII input file
    /*!
     * @param inputFile Name of the input file containing the phase XML data
     *                  to set up the object
     * @param id        ID of the phase in the input file. Defaults to the
     *                  empty string.
     */
    IdealGasPhase(const std::string& inputFile, const std::string& id = "");

    //! Construct and initialize an IdealGasPhase ThermoPhase object
    //! directly from an XML database
    /*!
     *  @param phaseRef XML phase node containing the description of the phase
     *  @param id     id attribute containing the name of the phase.
     *                (default is the empty string)
     */
    IdealGasPhase(XML_Node& phaseRef, const std::string& id = "");

    //! Copy Constructor
    /*!
     * Copy constructor for the object. Constructed
     * object will be a clone of this object, but will
     * also own all of its data.
     * This is a wrapper around the assignment operator
     *
     * @param right Object to be copied.
     */
    IdealGasPhase(const IdealGasPhase& right);

    //! Assignment operator
    /*!
     * Assignment operator for the object. Constructed
     * object will be a clone of this object, but will
     * also own all of its data.
     *
     * @param right Object to be copied.
     */
    IdealGasPhase& operator=(const IdealGasPhase& right);

    //! Duplicator from the ThermoPhase parent class
    /*!
     * Given a pointer to a ThermoPhase object, this function will
     * duplicate the ThermoPhase object and all underlying structures.
     * This is basically a wrapper around the inherited copy constructor.
     *
     * @return returns a pointer to a %ThermoPhase object, containing
     *      a copy of the current object
     */
    ThermoPhase* duplMyselfAsThermoPhase() const;

    //! Equation of state flag.
    /*!
     *  Returns the value cIdealGas, defined  in mix_defs.h.
     */
    virtual int eosType() const {
        return cIdealGas;
    }

    //! @name Molar Thermodynamic Properties of the Solution
    //! @{

    //! Return the Molar enthalpy. Units: J/kmol.
    /*!
     * For an ideal gas mixture,
     * \f[
     * \hat h(T) = \sum_k X_k \hat h^0_k(T),
     * \f]
     * and is a function only of temperature.
     * The standard-state pure-species enthalpies
     * \f$ \hat h^0_k(T) \f$ are computed by the species thermodynamic
     * property manager.
     *
     * \see SpeciesThermo
     */
    virtual doublereal enthalpy_mole() const {
        return GasConstant * temperature() * mean_X(&enthalpy_RT_ref()[0]);
    }

    /**
     * Molar internal energy. J/kmol. For an ideal gas mixture,
     * \f[
     * \hat u(T) = \sum_k X_k \hat h^0_k(T) - \hat R T,
     * \f]
     * and is a function only of temperature.
     * The reference-state pure-species enthalpies
     * \f$ \hat h^0_k(T) \f$ are computed by the species thermodynamic
     * property manager.
     * @see SpeciesThermo
     */
    virtual doublereal intEnergy_mole() const;

    /**
     * Molar entropy. Units: J/kmol/K.
     * For an ideal gas mixture,
     * \f[
     * \hat s(T, P) = \sum_k X_k \hat s^0_k(T) - \hat R \log (P/P^0).
     * \f]
     * The reference-state pure-species entropies
     * \f$ \hat s^0_k(T) \f$ are computed by the species thermodynamic
     * property manager.
     * @see SpeciesThermo
     */
    virtual doublereal entropy_mole() const;

    /**
     * Molar Gibbs free Energy for an ideal gas.
     * Units =  J/kmol.
     */
    virtual doublereal gibbs_mole() const;

    /**
     * Molar heat capacity at constant pressure. Units: J/kmol/K.
     * For an ideal gas mixture,
     * \f[
     * \hat c_p(t) = \sum_k \hat c^0_{p,k}(T).
     * \f]
     * The reference-state pure-species heat capacities
     * \f$ \hat c^0_{p,k}(T) \f$ are computed by the species thermodynamic
     * property manager.
     * @see SpeciesThermo
     */
    virtual doublereal cp_mole() const;

    /**
     * Molar heat capacity at constant volume. Units: J/kmol/K.
     * For an ideal gas mixture,
     * \f[ \hat c_v = \hat c_p - \hat R. \f]
     */
    virtual doublereal cv_mole() const;

    /**
     * @returns species translational/rotational specific heat at
     * constant volume.  Inferred from the species gas
     * constant and number of translational/rotational
     * degrees of freedom.  The translational/rotational
     * modes are assumed to be fully populated, and are
     * given by
     * \f[
     *   C^{tr}_{v,s} \equiv \frac{\partial e^{tr}_s}{\partial T} = \frac{5}{2} R_s
     * \f]
     * for diatomic molecules and
     * \f[
     *   C^{tr}_{v,s} \equiv \frac{\partial e^{tr}_s}{\partial T} = \frac{3}{2} R_s
     * \f]
     * for atoms.
     */
    virtual doublereal cv_tr(doublereal) const;

    /**
     * @returns species translational specific heat at constant volume.
     * Since the translational modes are assumed to be fully populated
     * this is simply
     * \f[
     *   C^{trans}_{v,s} \equiv \frac{\partial e^{trans}_s}{\partial T} = \frac{3}{2} R_s
     * \f]
     */
    virtual doublereal cv_trans() const;

    /**
     * @returns species rotational specific heat at constant volume.
     * By convention, we lump the translational/rotational components
     * \f[
     *   C^{tr}_{v,s} \equiv C^{trans}_{v,s} + C^{rot}_{v,s}
     * \f]
     * so then
     * \f[
     *   C^{rot}_{v,s} \equiv C^{tr}_{v,s} - C^{trans}_{v,s}
     * \f]
     */
    virtual doublereal cv_rot(double atomicity) const;

    /**
     * @returns species vibrational specific heat at
     * constant volume,
     * \f[
     *     C^{vib}_{v,s} = \frac{\partial e^{vib}_{v,s} }{\partial T}
     * \f]
     * where the species vibration energy \f$ e^{vib}_{v,s} \f$ is
     * - atom:
     *   0
     * - Diatomic:
     *   \f[ \frac{R_s \theta_{v,s}}{e^{\theta_{v,s}/T}-1} \f]
     * - General Molecule:
     *   \f[
     *       \sum_i \frac{R_s \theta_{v,s,i}}{e^{\theta_{v,s,i}/T}-1}
     *   \f]
     */
    virtual doublereal cv_vib(int k, doublereal T) const;

    //! @}
    //! @name Mechanical Equation of State
    //! @{

    /**
     * Pressure. Units: Pa.
     * For an ideal gas mixture,
     * \f[ P = n \hat R T. \f]
     */
    virtual doublereal pressure() const {
        return GasConstant * molarDensity() * temperature();
    }

    //! Set the pressure at constant temperature and composition.
    /*!
     *  Units: Pa.
     *   This method is implemented by setting the mass density to
     * \f[
     * \rho = \frac{P \overline W}{\hat R T }.
     * \f]
     *
     * @param p Pressure (Pa)
     */
    virtual void setPressure(doublereal p) {
        setDensity(p * meanMolecularWeight() / (GasConstant * temperature()));
    }

    //! Returns  the isothermal compressibility. Units: 1/Pa.
    /**
     * The isothermal compressibility is defined as
     * \f[
     * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
     * \f]
     *  For ideal gases it's equal to the inverse of the pressure
     */
    virtual doublereal isothermalCompressibility() const {
        return 1.0 / pressure();
    }

    //! Return the volumetric thermal expansion coefficient. Units: 1/K.
    /*!
     * The thermal expansion coefficient is defined as
     * \f[
     * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
     * \f]
     * For ideal gases, it's equal to the inverse of the temperature.
     */
    virtual doublereal thermalExpansionCoeff() const {
        return 1.0 / temperature();
    }

    //@}

    /**
     * @name Chemical Potentials and Activities
     *
     * The activity \f$a_k\f$ of a species in solution is
     * related to the chemical potential by
     * \f[
     *  \mu_k(T,P,X_k) = \mu_k^0(T,P)
     * + \hat R T \log a_k.
     *  \f]
     * The quantity \f$\mu_k^0(T,P)\f$ is
     * the standard state chemical potential at unit activity.
     * It may depend on the pressure and the temperature. However,
     * it may not depend on the mole fractions of the species
     * in the solution.
     *
     * The activities are related to the generalized
     * concentrations, \f$\tilde C_k\f$, and standard
     * concentrations, \f$C^0_k\f$, by the following formula:
     *
     *  \f[
     *  a_k = \frac{\tilde C_k}{C^0_k}
     *  \f]
     * The generalized concentrations are used in the kinetics classes
     * to describe the rates of progress of reactions involving the
     * species. Their formulation depends upon the specification
     * of the rate constants for reaction, especially the units used
     * in specifying the rate constants. The bridge between the
     * thermodynamic equilibrium expressions that use a_k and the
     * kinetics expressions which use the generalized concentrations
     * is provided by the multiplicative factor of the
     * standard concentrations.
     * @{
     */

    //! This method returns the array of generalized concentrations.
    /*!
     *  For an ideal gas mixture, these are simply the actual concentrations.
     *
     * @param c Output array of generalized concentrations. The
     *           units depend upon the implementation of the
     *           reaction rate expressions within the phase.
     */
    virtual void getActivityConcentrations(doublereal* c) const {
        getConcentrations(c);
    }

    //! Returns the standard concentration \f$ C^0_k \f$, which is used to normalize
    //! the generalized concentration.
    /*!
     * This is defined as the concentration  by which the generalized
     * concentration is normalized to produce the activity.
     * In many cases, this quantity will be the same for all species in a phase.
     * Since the activity for an ideal gas mixture is
     * simply the mole fraction, for an ideal gas \f$ C^0_k = P/\hat R T \f$.
     *
     * @param k Optional parameter indicating the species. The default
     *          is to assume this refers to species 0.
     * @return
     *   Returns the standard Concentration in units of m3 kmol-1.
     */
    virtual doublereal standardConcentration(size_t k = 0) const;

    //! Returns the natural logarithm of the standard
    //! concentration of the kth species
    /*!
     * @param k    index of the species. (defaults to zero)
     */
    virtual doublereal logStandardConc(size_t k = 0) const;

    //! Get the array of non-dimensional activity coefficients at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     *  For ideal gases, the activity coefficients are all equal to one.
     *
     * @param ac Output vector of activity coefficients. Length: m_kk.
     */
    virtual void getActivityCoefficients(doublereal* ac) const;

    //@}
    /// @name Partial Molar Properties of the Solution
    //@{

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

    //!  Get the species partial molar enthalpies. Units: J/kmol.
    /*!
     * @param hbar    Output vector of species partial molar enthalpies.
     *                Length: m_kk. units are J/kmol.
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;

    //! Get the species partial molar entropies. Units: J/kmol/K.
    /*!
     * @param sbar    Output vector of species partial molar entropies.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;

    //! Get the species partial molar enthalpies. Units: J/kmol.
    /*!
     * @param ubar    Output vector of species partial molar internal energies.
     *                Length = m_kk. units are J/kmol.
     */
    virtual void getPartialMolarIntEnergies(doublereal* ubar) const;

    //! Get the partial molar heat capacities Units: J/kmol/K
    /*!
     * @param cpbar   Output vector of species partial molar heat capacities at constant pressure.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const;

    //! Get the species partial molar volumes. Units: m^3/kmol.
    /*!
     *  @param vbar   Output vector of species partial molar volumes.
     *                Length = m_kk. units are m^3/kmol.
     */
    virtual void getPartialMolarVolumes(doublereal* vbar) const;

    //@}
    /// @name  Properties of the Standard State of the Species in the Solution
    //@{

    //! Get the array of chemical potentials at unit activity for the
    //! species standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * These are the standard state chemical potentials \f$ \mu^0_k(T,P)
     * \f$. The values are evaluated at the current
     * temperature and pressure of the solution
     *
     * @param mu  Output vector of chemical potentials.
     *            Length: m_kk.
     */
    virtual void getStandardChemPotentials(doublereal* mu) const;

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

    //! Get the nondimensional Gibbs functions for the species
    //! standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param grt  Output vector of nondimensional standard state gibbs free energies
     *             Length: m_kk.
     */
    virtual void getGibbs_RT(doublereal* grt) const;

    //! Get the Gibbs functions for the standard
    //! state of the species at the current <I>T</I> and <I>P</I> of the solution
    /*!
     * Units are Joules/kmol
     * @param gpure  Output vector of  standard state gibbs free energies
     *               Length: m_kk.
     */
    virtual void getPureGibbs(doublereal* gpure) const;

    //!  Returns the vector of nondimensional Internal Energies  of the standard
    //!  state species at the current <I>T</I> and <I>P</I> of the solution
    /*!
     * @param urt  output vector of nondimensional standard state internal energies
     *             of the species. Length: m_kk.
     */
    virtual void getIntEnergy_RT(doublereal* urt) const;

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

    //@}
    /// @name Thermodynamic Values for the Species Reference States
    //@{

    //!  Returns the vector of nondimensional
    //!  enthalpies of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    /*!
     * @param hrt     Output vector containing the nondimensional reference state
     *                enthalpies.  Length: m_kk.
     */
    virtual void getEnthalpy_RT_ref(doublereal* hrt) const;

#ifdef H298MODIFY_CAPABILITY

    virtual void modifyOneHf298SS(const size_t& k, const doublereal Hf298New) {
        m_spthermo->modifyOneHf298(k, Hf298New);
        m_tlast += 0.0001234;
    }
#endif
    //!  Returns the vector of nondimensional
    //!  Gibbs Free Energies of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    /*!
     * @param grt     Output vector containing the nondimensional reference state
     *                Gibbs Free energies.  Length: m_kk.
     */
    virtual void getGibbs_RT_ref(doublereal* grt) const;

    //!  Returns the vector of the
    //!  gibbs function of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    /*!
     *  units = J/kmol
     *
     * @param g       Output vector containing the  reference state
     *                Gibbs Free energies.  Length: m_kk. Units: J/kmol.
     */
    virtual void getGibbs_ref(doublereal* g) const;

    //!  Returns the vector of nondimensional
    //!  entropies of the reference state at the current temperature
    //!  of the solution and the reference pressure for each species.
    /*!
     * @param er      Output vector containing the nondimensional reference state
     *                entropies.  Length: m_kk.
     */
    virtual void getEntropy_R_ref(doublereal* er) const;

    //! Returns the vector of nondimensional
    //!  internal Energies of the reference state at the current temperature
    //!  of the solution and the reference pressure for each species.
    /*!
     * @param urt    Output vector of nondimensional reference state
     *               internal energies of the species.
     *               Length: m_kk
     */
    virtual void getIntEnergy_RT_ref(doublereal* urt) const;

    //!  Returns the vector of nondimensional
    //!  constant pressure heat capacities of the reference state
    //!  at the current temperature of the solution
    //!  and reference pressure for each species.
    /*!
     * @param cprt   Output vector of nondimensional reference state
     *               heat capacities at constant pressure for the species.
     *               Length: m_kk
     */
    virtual void getCp_R_ref(doublereal* cprt) const;

    //!  Get the molar volumes of the species standard states at the current
    //!  <I>T</I> and <I>P_ref</I> of the solution.
    /*!
     * units = m^3 / kmol
     *
     * @param vol     Output vector containing the standard state volumes.
     *                Length: m_kk.
     */
    virtual void getStandardVolumes_ref(doublereal* vol) const;

    //@}
    /// @name NonVirtual Internal methods to Return References to Reference State Thermo
    //@{

    //! Returns a reference to the dimensionless reference state enthalpy vector.
    /*!
     * This function is part of the layer that checks/recalculates the reference
     * state thermo functions.
     */
    const vector_fp& enthalpy_RT_ref() const {
        _updateThermo();
        return m_h0_RT;
    }

    //! Returns a reference to the dimensionless reference state Gibbs free energy vector.
    /*!
     * This function is part of the layer that checks/recalculates the reference
     * state thermo functions.
     */
    const vector_fp& gibbs_RT_ref() const {
        _updateThermo();
        return m_g0_RT;
    }

    //! Returns a reference to the dimensionless reference state Entropy vector.
    /*!
     * This function is part of the layer that checks/recalculates the reference
     * state thermo functions.
     */
    const vector_fp& entropy_R_ref() const {
        _updateThermo();
        return m_s0_R;
    }

    //! Returns a reference to the dimensionless reference state Heat Capacity vector.
    /*!
     * This function is part of the layer that checks/recalculates the reference
     * state thermo functions.
     */
    const vector_fp& cp_R_ref() const {
        _updateThermo();
        return m_cp0_R;
    }

    //@}

    //! Initialize the ThermoPhase object after all species have been set up
    /*!
     * @internal Initialize.
     *
     * This method performs any initialization required after all
     * species have been added. For example, it is used to
     * resize internal work arrays that must have an entry for
     * each species.
     * This method is called from ThermoPhase::initThermoXML(),
     * which is called from importPhase(),
     * just prior to returning from the function, importPhase().
     *
     * @see importCTML.cpp
     */
    virtual void initThermo();

    //! Method used by the ChemEquil equilibrium solver.
    /*!
     * @internal
     *
     * Set mixture to an equilibrium state consistent with specified
     * element potentials and temperature.
     * It sets the state such that the chemical potentials satisfy
     * \f[ \frac{\mu_k}{\hat R T} = \sum_m A_{k,m}
     * \left(\frac{\lambda_m} {\hat R T}\right) \f] where
     * \f$ \lambda_m \f$ is the element potential of element m. The
     * temperature is unchanged.  Any phase (ideal or not) that
     * implements this method can be equilibrated by ChemEquil.
     *
     * @param lambda_RT vector of non-dimensional element potentials
     *                  \f[ \lambda_m/RT \f].
     */
    virtual void setToEquilState(const doublereal* lambda_RT);

protected:
    //! Reference state pressure
    /*!
     *  Value of the reference state pressure in Pascals.
     *  All species must have the same reference state pressure.
     */
    doublereal m_p0;

    //! last value of the temperature processed by reference state
    mutable doublereal m_tlast;

    //! Temporary storage for log of p/RT
    mutable doublereal m_logc0;

    //! Temporary storage for dimensionless reference state enthalpies
    mutable vector_fp m_h0_RT;

    //! Temporary storage for dimensionless reference state heat capacities
    mutable vector_fp m_cp0_R;

    //! Temporary storage for dimensionless reference state gibbs energies
    mutable vector_fp m_g0_RT;

    //! Temporary storage for dimensionless reference state entropies
    mutable vector_fp m_s0_R;

    mutable vector_fp m_expg0_RT;

    //! Temporary array containing internally calculated partial pressures
    mutable vector_fp m_pp;

private:
    //! Update the species reference state thermodynamic functions
    /*!
     *  This method is called each time a thermodynamic property is requested,
     *  to check whether the internal species properties within the object
     *  need to be updated. Currently, this updates the species thermo
     *  polynomial values for the current value of the temperature. A check is
     *  made to see if the temperature has changed since the last evaluation.
     *  This object does not contain any persistent data that depends on the
     *  concentration, that needs to be updated. The state object modifies its
     *  concentration dependent information at the time the setMoleFractions()
     *  (or equivalent) call is made.
     */
    void _updateThermo() const;
};
}

#endif
