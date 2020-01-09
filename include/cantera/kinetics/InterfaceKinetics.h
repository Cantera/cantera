/**
 * @file InterfaceKinetics.h
 * @ingroup chemkinetics
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IFACEKINETICS_H
#define CT_IFACEKINETICS_H

#include "Kinetics.h"
#include "Reaction.h"
#include "RateCoeffMgr.h"

namespace Cantera
{

// forward declarations
class SurfPhase;
class ImplicitSurfChem;

//! A kinetics manager for heterogeneous reaction mechanisms. The reactions are
//! assumed to occur at a 2D interface between two 3D phases.
/*!
 * There are some important additions to the behavior of the kinetics class due
 * to the presence of multiple phases and a heterogeneous interface. If a
 * reactant phase doesn't exists, i.e., has a mole number of zero, a
 * heterogeneous reaction can not proceed from reactants to products. Note it
 * could perhaps proceed from products to reactants if all of the product phases
 * exist.
 *
 * In order to make the determination of whether a phase exists or not actually
 * involves the specification of additional information to the kinetics object.,
 * which heretofore has only had access to intrinsic field information about the
 * phases (i.e., temperature pressure, and mole fraction).
 *
 * The extrinsic specification of whether a phase exists or not must be
 * specified on top of the intrinsic calculation of the reaction rate. This
 * class carries a set of booleans indicating whether a phase in the
 * heterogeneous mechanism exists or not.
 *
 * Additionally, the class carries a set of booleans around indicating whether a
 * product phase is stable or not. If a phase is not thermodynamically stable,
 * it may be the case that a particular reaction in a heterogeneous mechanism
 * will create a product species in the unstable phase. However, other reactions
 * in the mechanism will destruct that species. This may cause oscillations in
 * the formation of the unstable phase from time step to time step within a ODE
 * solver, in practice. In order to avoid this situation, a set of booleans is
 * tracked which sets the stability of a phase. If a phase is deemed to be
 * unstable, then species in that phase will not be allowed to be birthed by the
 * kinetics operator. Nonexistent phases are deemed to be unstable by default,
 * but this can be changed.
 *
 *  @ingroup chemkinetics
 */
class InterfaceKinetics : public Kinetics
{
public:
    //! Constructor
    /*!
     * @param thermo The optional parameter may be used to initialize the object
     *               with one ThermoPhase object.
     *               HKM Note -> Since the interface kinetics object will
     *               probably require multiple ThermoPhase objects, this is
     *               probably not a good idea to have this parameter.
     */
    InterfaceKinetics(thermo_t* thermo = 0);

    virtual ~InterfaceKinetics();

    virtual std::string kineticsType() const {
        return "Surf";
    }

    //! Set the electric potential in the nth phase
    /*!
     * @param n phase Index in this kinetics object.
     * @param V Electric potential (volts)
     */
    void setElectricPotential(int n, doublereal V);

    //! @name Reaction Rates Of Progress
    //! @{

    //! Equilibrium constant for all reactions including the voltage term
    /*!
     *   Kc = exp(deltaG/RT)
     *
     *   where deltaG is the electrochemical potential difference between
     *   products minus reactants.
     */
    virtual void getEquilibriumConstants(doublereal* kc);

    //! values needed to convert from exchange current density to surface
    //! reaction rate.
    void updateExchangeCurrentQuantities();

    virtual void getDeltaGibbs(doublereal* deltaG);

    virtual void getDeltaElectrochemPotentials(doublereal* deltaM);
    virtual void getDeltaEnthalpy(doublereal* deltaH);
    virtual void getDeltaEntropy(doublereal* deltaS);

    virtual void getDeltaSSGibbs(doublereal* deltaG);
    virtual void getDeltaSSEnthalpy(doublereal* deltaH);
    virtual void getDeltaSSEntropy(doublereal* deltaS);

    //! @}
    //! @name Reaction Mechanism Informational Query Routines
    //! @{

    virtual void getActivityConcentrations(doublereal* const conc);

    //! Return the charge transfer rxn Beta parameter for the ith reaction
    /*!
     *  Returns the beta parameter for a charge transfer reaction. This
     *  parameter is not important for non-charge transfer reactions.
     *  Note, the parameter defaults to zero. However, a value of 0.5
     *  should be supplied for every charge transfer reaction if
     *  no information is known, as a value of 0.5 pertains to a
     *  symmetric transition state. The value can vary between 0 to 1.
     *
     *  @param irxn Reaction number in the kinetics mechanism
     *
     *  @return Beta parameter. This defaults to zero, even for charge
     *    transfer reactions.
     */
    doublereal electrochem_beta(size_t irxn) const;

    virtual bool isReversible(size_t i) {
        if (std::find(m_revindex.begin(), m_revindex.end(), i)
                < m_revindex.end()) {
            return true;
        } else {
            return false;
        }
    }

    virtual void getFwdRateConstants(doublereal* kfwd);
    virtual void getRevRateConstants(doublereal* krev,
                                     bool doIrreversible = false);

    //! Return effective preexponent for the specified reaction
    /*!
     *  Returns effective preexponent, accounting for surface coverage
     *  dependencies.
     *
     *  @param irxn Reaction number in the kinetics mechanism
     *  @return Effective preexponent
     */
    double effectivePreExponentialFactor(size_t irxn) {
        return m_rates.effectivePreExponentialFactor(irxn);
    }

    //! Return effective activation energy for the specified reaction
    /*!
     *  Returns effective activation energy, accounting for surface coverage
     *  dependencies.
     *
     *  @param irxn Reaction number in the kinetics mechanism
     *  @return Effective activation energy divided by the gas constant
     */
    double effectiveActivationEnergy_R(size_t irxn) {
       return m_rates.effectiveActivationEnergy_R(irxn);
    }

    //! Return effective temperature exponent for the specified reaction
    /*!
     *  Returns effective temperature exponent, accounting for surface coverage
     *  dependencies. Current parameterization in SurfaceArrhenius does not
     *  change this parameter with the change in surface coverages.
     *
     *  @param irxn Reaction number in the kinetics mechanism
     *  @return Effective temperature exponent
     */
    double effectiveTemperatureExponent(size_t irxn) {
       return m_rates.effectiveTemperatureExponent(irxn);
    }

    //! @}
    //! @name Reaction Mechanism Construction
    //! @{

    //!  Add a phase to the kinetics manager object.
    /*!
     * This must be done before the function init() is called or
     * before any reactions are input.
     *
     * This function calls Kinetics::addPhase(). It also sets the following
     * fields:
     *
     *        m_phaseExists[]
     *
     * @param thermo    Reference to the ThermoPhase to be added.
     */
    virtual void addPhase(thermo_t& thermo);

    virtual void init();
    virtual void resizeSpecies();
    virtual bool addReaction(shared_ptr<Reaction> r);
    virtual void modifyReaction(size_t i, shared_ptr<Reaction> rNew);
    //! @}

    //! Internal routine that updates the Rates of Progress of the reactions
    /*!
     *  This is actually the guts of the functionality of the object
     */
    virtual void updateROP();

    //! Update properties that depend on temperature
    /*!
     *  Current objects that this function updates:
     *       m_kdata->m_logtemp
     *       m_kdata->m_rfn
     *       m_rates.
     *       updateKc();
     */
    void _update_rates_T();

    //! Update properties that depend on the electric potential
    void _update_rates_phi();

    //! Update properties that depend on the species mole fractions and/or
    //! concentration,
    /*!
     * This method fills out the array of generalized concentrations by calling
     * method getActivityConcentrations for each phase, which classes
     * representing phases should overload to return the appropriate quantities.
     */
    void _update_rates_C();

    //! Advance the surface coverages in time
    /*!
     * This method carries out a time-accurate advancement of the
     * surface coverages for a specified amount of time.
     *
     *  \f[
     *    \dot {\theta}_k = \dot s_k (\sigma_k / s_0)
     *  \f]
     *
     * @param tstep  Time value to advance the surface coverages
     * @param rtol   The relative tolerance for the integrator
     * @param atol   The absolute tolerance for the integrator
     * @param maxStepSize   The maximum step-size the integrator is allowed to take.
     *                      If zero, this option is disabled.
     * @param maxSteps   The maximum number of time-steps the integrator can take.
     *                   If not supplied, uses the default value in CVodeIntegrator (20000).
     * @param maxErrTestFails   the maximum permissible number of error test failures
     *                           If not supplied, uses the default value in CVODES (7).
     */
    void advanceCoverages(doublereal tstep, double rtol=1.e-7,
                          double atol=1.e-14, double maxStepSize=0,
                          size_t maxSteps=20000, size_t maxErrTestFails=7);

    //! Solve for the pseudo steady-state of the surface problem
    /*!
     * This is the same thing as the advanceCoverages() function,
     * but at infinite times.
     *
     * Note, a direct solve is carried out under the hood here,
     * to reduce the computational time.
     *
     * @param ifuncOverride One of the values defined in @ref solvesp_methods.
     *         The default is -1, which means that the program will decide.
     * @param timeScaleOverride When a pseudo transient is selected this value
     *             can be used to override the default time scale for
     *             integration which is one. When SFLUX_TRANSIENT is used, this
     *             is equal to the time over which the equations are integrated.
     *             When SFLUX_INITIALIZE is used, this is equal to the time used
     *             in the initial transient algorithm, before the equation
     *             system is solved directly.
     */
    void solvePseudoSteadyStateProblem(int ifuncOverride = -1,
                                       doublereal timeScaleOverride = 1.0);

    void setIOFlag(int ioFlag);

    //! Update the standard state chemical potentials and species equilibrium
    //! constant entries
    /*!
     *  Virtual because it is overridden when dealing with experimental open
     *  circuit voltage overrides
     */
    virtual void updateMu0();

    //! Update the equilibrium constants and stored electrochemical potentials
    //! in molar units for all reversible reactions and for all species.
    /*!
     *  Irreversible reactions have their equilibrium constant set
     *  to zero. For reactions involving charged species the equilibrium
     *  constant is adjusted according to the electrostatic potential.
     */
    void updateKc();

    //! Apply modifications for the forward reaction rate for interfacial charge
    //! transfer reactions
    /*!
     * For reactions that transfer charge across a potential difference,
     * the activation energies are modified by the potential difference.
     * (see, for example, ...). This method applies this correction.
     *
     * @param kfwd  Vector of forward reaction rate constants on which to have
     *              the voltage correction applied
     */
    void applyVoltageKfwdCorrection(doublereal* const kfwd);

    //! When an electrode reaction rate is optionally specified in terms of its
    //! exchange current density, adjust kfwd to the standard reaction rate
    //! constant form and units. When the BV reaction types are used, keep the
    //! exchange current density form.
    /*!
     *  For a reaction rate constant that was given in units of Amps/m2
     *  (exchange current density formulation with iECDFormulation == true),
     *  convert the rate to kmoles/m2/s.
     *
     *  For a reaction rate constant that was given in units of kmol/m2/sec when
     *  the reaction type is a Butler-Volmer form, convert it to exchange
     *  current density form (amps/m2).
     *
     * @param kfwd  Vector of forward reaction rate constants, given in either
     *              normal form or in exchange current density form.
     */
    void convertExchangeCurrentDensityFormulation(doublereal* const kfwd);

    //! Set the existence of a phase in the reaction object
    /*!
     *  Tell the kinetics object whether a phase in the object exists. This is
     *  actually an extrinsic specification that must be carried out on top of
     *  the intrinsic calculation of the reaction rate. The routine will also
     *  flip the IsStable boolean within the kinetics object as well.
     *
     *  @param iphase  Index of the phase. This is the order within the
     *      internal thermo vector object
     *  @param exists  Boolean indicating whether the phase exists or not
     */
    void setPhaseExistence(const size_t iphase, const int exists);

    //! Set the stability of a phase in the reaction object
    /*!
     *  Tell the kinetics object whether a phase in the object is stable.
     *  Species in an unstable phase will not be allowed to have a positive
     *  rate of formation from this kinetics object. This is actually an
     *  extrinsic specification that must be carried out on top of the
     *  intrinsic calculation of the reaction rate.
     *
     *  While conceptually not needed since kinetics is consistent with thermo
     *  when taken as a whole, in practice it has found to be very useful to
     *  turn off the creation of phases which shouldn't be forming. Typically
     *  this can reduce the oscillations in phase formation and destruction
     *  which are observed.
     *
     *  @param iphase  Index of the phase. This is the order within the
     *      internal thermo vector object
     *  @param isStable Flag indicating whether the phase is stable or not
     */
    void setPhaseStability(const size_t iphase, const int isStable);

    //! Gets the phase existence int for the ith phase
    /*!
     * @param iphase  Phase Id
     * @return The int specifying whether the kinetics object thinks the phase
     *         exists or not. If it exists, then species in that phase can be
     *         a reactant in reactions.
     */
    int phaseExistence(const size_t iphase) const;

    //! Gets the phase stability int for the ith phase
    /*!
     * @param iphase  Phase Id
     * @return The int specifying whether the kinetics object thinks the phase
     *         is stable with nonzero mole numbers. If it stable, then the
     *         kinetics object will allow for rates of production of of
     *         species in that phase that are positive.
     */
    int phaseStability(const size_t iphase) const;

    //! @deprecated To be removed after Cantera 2.5.
    virtual void determineFwdOrdersBV(ElectrochemicalReaction& r, vector_fp& fwdFullorders);

protected:
    //! Build a SurfaceArrhenius object from a Reaction, taking into account
    //! the possible sticking coefficient form and coverage dependencies
    //! @param i  Reaction number
    //! @param r  Reaction object containing rate coefficient parameters
    //! @param replace  True if replacing an existing reaction
    SurfaceArrhenius buildSurfaceArrhenius(size_t i, InterfaceReaction& r,
                                           bool replace);

    //! Temporary work vector of length m_kk
    vector_fp m_grt;

    //! List of reactions numbers which are reversible reactions
    /*!
     *  This is a vector of reaction numbers. Each reaction in the list is
     *  reversible. Length = number of reversible reactions
     */
    std::vector<size_t> m_revindex;

    //! Templated class containing the vector of reactions for this interface
    /*!
     *  The templated class is described in RateCoeffMgr.h
     *  The class SurfaceArrhenius is described in RxnRates.h
     */
    Rate1<SurfaceArrhenius> m_rates;

    bool m_redo_rates;

    //! Vector of irreversible reaction numbers
    /*!
     * vector containing the reaction numbers of irreversible reactions.
     */
    std::vector<size_t> m_irrev;

    //! Array of concentrations for each species in the kinetics mechanism
    /*!
     * An array of generalized concentrations \f$ C_k \f$ that are defined
     * such that \f$ a_k = C_k / C^0_k, \f$ where \f$ C^0_k \f$ is a standard
     * concentration/ These generalized concentrations are used by this
     * kinetics manager class to compute the forward and reverse rates of
     * elementary reactions. The "units" for the concentrations of each phase
     * depend upon the implementation of kinetics within that phase. The order
     * of the species within the vector is based on the order of listed
     * ThermoPhase objects in the class, and the order of the species within
     * each ThermoPhase class.
     */
    vector_fp m_conc;

    //! Array of activity concentrations for each species in the kinetics object
    /*!
     * An array of activity concentrations \f$ Ca_k \f$ that are defined
     * such that \f$ a_k = Ca_k / C^0_k, \f$ where \f$ C^0_k \f$ is a standard
     * concentration. These activity concentrations are used by this
     * kinetics manager class to compute the forward and reverse rates of
     * elementary reactions. The "units" for the concentrations of each phase
     * depend upon the implementation of kinetics within that phase. The order
     * of the species within the vector is based on the order of listed
     * ThermoPhase objects in the class, and the order of the species within
     * each ThermoPhase class.
     */
    vector_fp m_actConc;

    //! Vector of standard state chemical potentials for all species
    /*!
     * This vector contains a temporary vector of standard state chemical
     * potentials for all of the species in the kinetics object
     *
     * Length = m_kk. Units = J/kmol.
     */
    vector_fp m_mu0;

    //! Vector of chemical potentials for all species
    /*!
     * This vector contains a vector of chemical potentials for all of the
     * species in the kinetics object
     *
     * Length = m_kk. Units = J/kmol.
     */
    vector_fp m_mu;

    //! Vector of standard state electrochemical potentials modified by a
    //! standard concentration term.
    /*!
     * This vector contains a temporary vector of standard state electrochemical
     * potentials + RTln(Cs) for all of the species in the kinetics object
     *
     * In order to get the units correct for the concentration equilibrium
     * constant, each species needs to have an RT ln(Cs)  added to its
     * contribution to the equilibrium constant Cs is the standard concentration
     * for the species. Frequently, for solid species, Cs is equal to 1.
     * However, for gases Cs is P/RT. Length = m_kk. Units = J/kmol.
     */
    vector_fp m_mu0_Kc;

    //! Vector of phase electric potentials
    /*!
     * Temporary vector containing the potential of each phase in the kinetics
     * object. length = number of phases. Units = Volts.
     */
    vector_fp m_phi;

    //! Vector of potential energies due to Voltages
    /*!
     * Length is the number of species in kinetics mech. It's used to store the
     * potential energy due to the voltage.
     */
    vector_fp m_pot;

    //! Storage for the net electric energy change due to reaction.
    /*!
     * Length is number of reactions. It's used to store the net electric
     * potential energy change due to the reaction.
     *
     *  deltaElectricEnergy_[jrxn] = sum_i ( F V_i z_i nu_ij)
     */
    vector_fp deltaElectricEnergy_;

    //! Pointer to the single surface phase
    SurfPhase* m_surf;

    //! Pointer to the Implicit surface chemistry object
    /*!
     * Note this object is owned by this InterfaceKinetics object. It may only
     * be used to solve this single InterfaceKinetics object's surface problem
     * uncoupled from other surface phases.
     */
    ImplicitSurfChem* m_integrator;

    //! Electrochemical transfer coefficient for the forward direction
    /*!
     *   Electrochemical transfer coefficient for all reactions that have
     *   transfer reactions the reaction is given by m_ctrxn[i]
     */
    vector_fp m_beta;

    //! Vector of reaction indexes specifying the id of the charge transfer
    //! reactions in the mechanism
    /*!
     *  Vector of reaction indices which involve charge transfers. This provides
     *  an index into the m_beta and m_ctrxn_BVform array.
     *
     *        irxn = m_ctrxn[i]
     */
    std::vector<size_t> m_ctrxn;

    //! Vector of Reactions which follow the Butler-Volmer methodology for specifying the
    //! exchange current density first. Then, the other forms are specified based on this form.
    /*!
     *     Length is equal to the number of reactions with charge transfer coefficients, m_ctrxn[]
     *
     *    m_ctrxn_BVform[i] = 0;  This means that the irxn reaction is calculated via the standard forward
     *                            and reverse reaction rates
     *    m_ctrxn_BVform[i] = 1;  This means that the irxn reaction is calculated via the BV format
     *                            directly.
     *    m_ctrxn_BVform[i] = 2;  this means that the irxn reaction is calculated via the BV format
     *                            directly, using concentrations instead of activity concentrations.
     * @deprecated To be removed after Cantera 2.5.
     */
    std::vector<size_t> m_ctrxn_BVform;

    //! Vector of booleans indicating whether the charge transfer reaction rate constant
    //! is described by an exchange current density rate constant expression
    /*!
     *   Length is equal to the number of reactions with charge transfer coefficients, m_ctrxn[]
     *
     *   m_ctrxn_ecdf[irxn] = 0   This means that the rate coefficient calculator will calculate
     *                            the rate constant as a chemical forward rate constant, a standard format.
     *   m_ctrxn_ecdf[irxn] = 1   this means that the rate coefficient calculator will calculate
     *                            the rate constant as an exchange current density rate constant expression.
     */
    vector_int m_ctrxn_ecdf;

    //! Vector of standard concentrations
    /*!
     *   Length number of kinetic species
     *   units depend on the definition of the standard concentration within each phase
     */
    vector_fp m_StandardConc;

    //!  Vector of delta G^0, the standard state Gibbs free energies for each reaction
    /*!
     *    Length is the number of reactions
     *    units are Joule kmol-1
     */
    vector_fp m_deltaG0;

    //! Vector of deltaG[] of reaction, the delta Gibbs free energies for each reaction
    /*!
     *    Length is the number of reactions
     *    units are Joule kmol-1
     */
    vector_fp m_deltaG;

    //! Vector of the products of the standard concentrations of the reactants
    /*!
     *   Units vary wrt what the units of the standard concentrations are
     *   Length = number of reactions.
     */
    vector_fp m_ProdStanConcReac;

    bool m_ROP_ok;

    //! Current temperature of the data
    doublereal m_temp;

    //! Current log of the temperature
    doublereal m_logtemp;

    //! Boolean flag indicating whether any reaction in the mechanism
    //! has a coverage dependent forward reaction rate
    /*!
     *   If this is true, then the coverage dependence is multiplied into
     *   the forward reaction rates constant
     */
    bool m_has_coverage_dependence;

    //! Boolean flag indicating whether any reaction in the mechanism
    //! has a beta electrochemical parameter.
    /*!
     *  If this is true, the Butler-Volmer correction is applied
     *  to the forward reaction rate for those reactions.
     *
     *    fac = exp ( - beta * (delta_phi))
     */
    bool m_has_electrochem_rxns;

    //! Boolean flag indicating whether any reaction in the mechanism
    //! is described by an exchange current density expression
    /*!
     *  If this is true, the standard state Gibbs free energy of the reaction
     *  and the product of the reactant standard concentrations must be
     *  precalculated in order to calculate the rate constant.
     */
    bool m_has_exchange_current_density_formulation;

    //! Int flag to indicate that some phases in the kinetics mechanism are
    //! non-existent.
    /*!
     *  We change the ROP vectors to make sure that non-existent phases are
     *  treated correctly in the kinetics operator. The value of this is equal
     *  to the number of phases which don't exist.
     */
    int m_phaseExistsCheck;

    //!  Vector of booleans indicating whether phases exist or not
    /*!
     *  Vector of booleans indicating whether a phase exists or not. We use this
     *  to set the ROP's so that unphysical things don't happen. For example, a
     *  reaction can't go in the forwards direction if a phase in which a
     *  reactant is present doesn't exist. Because InterfaceKinetics deals with
     *  intrinsic quantities only normally, nowhere else is this extrinsic
     *  concept introduced except here.
     *
     *  length = number of phases in the object. By default all phases exist.
     */
    std::vector<bool> m_phaseExists;

    //! Vector of int indicating whether phases are stable or not
    /*!
     *  Vector of booleans indicating whether a phase is stable or not under
     *  the current conditions. We use this to set the ROP's so that
     *  unphysical things don't happen.
     *
     *  length = number of phases in the object. By default all phases are stable.
     */
    vector_int m_phaseIsStable;

    //! Vector of vector of booleans indicating whether a phase participates in
    //! a reaction as a reactant
    /*!
     *  m_rxnPhaseIsReactant[j][p] indicates whether a species in phase p
     *  participates in reaction j as a reactant.
     */
    std::vector<std::vector<bool> > m_rxnPhaseIsReactant;

    //! Vector of vector of booleans indicating whether a phase participates in a
    //! reaction as a product
    /*!
     *  m_rxnPhaseIsReactant[j][p] indicates whether a species in phase p
     *  participates in reaction j as a product.
     */
    std::vector<std::vector<bool> > m_rxnPhaseIsProduct;

    //! Values used for converting sticking coefficients into rate constants
    struct StickData {
        size_t index; //!< index of the sticking reaction in the full reaction list
        double order; //!< exponent applied to site density term
        double multiplier; //!< multiplicative factor in rate expression
        bool use_motz_wise; //!< 'true' if Motz & Wise correction is being used
    };

    //! Data for sticking reactions
    std::vector<StickData> m_stickingData;

    void applyStickingCorrection(double T, double* kf);

    int m_ioFlag;

    //! Number of dimensions of reacting phase (2 for InterfaceKinetics, 1 for
    //! EdgeKinetics)
    size_t m_nDim;
};

}

#endif
