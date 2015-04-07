/**
 * @file Kinetics.h
 *  Base class for kinetics managers and also contains the kineticsmgr
 *  module documentation (see \ref  kineticsmgr and class
 *  \link Cantera::Kinetics Kinetics\endlink).
 */

// Copyright 2001-2004  California Institute of Technology

#ifndef CT_KINETICS_H
#define CT_KINETICS_H

#include "cantera/base/ctexceptions.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/mix_defs.h"

namespace Cantera
{

// forward references
class ReactionData;

/**
 * @defgroup chemkinetics Chemical Kinetics
 */

/// @defgroup kineticsmgr Kinetics Managers
/// @section kinmodman Models and Managers
///
/// A kinetics manager is a C++ class that implements a kinetics
/// model; a kinetics model is a set of mathematical equation
/// describing how various kinetic quantities are to be computed --
/// reaction rates, species production rates, etc. Many different
/// kinetics models might be defined to handle different types of
/// kinetic processes. For example, one kinetics model might use
/// expressions valid for elementary reactions in ideal gas
/// mixtures. It might, for example, require the reaction orders
/// to be integral and equal to the forward stoichiometric
/// coefficients, require that each reaction be reversible with a
/// reverse rate satisfying detailed balance, include
/// pressure-dependent unimolecular reactions, etc. Another
/// kinetics model might be designed for heterogeneous chemistry
/// at interfaces, and might allow empirical reaction orders,
/// coverage-dependent activation energies, irreversible
/// reactions, and include effects of potential differences across
/// the interface on reaction rates.
///
/// A kinetics manager implements a kinetics model. Since the
/// model equations may be complex and expensive to evaluate, a
/// kinetics manager may adopt various strategies to 'manage' the
/// computation and evaluate the expressions efficiently. For
/// example, if there are rate coefficients or other quantities
/// that depend only on temperature, a manager class may choose to
/// store these quantities internally, and re-evaluate them only
/// when the temperature has actually changed. Or a manager
/// designed for use with reaction mechanisms with a few repeated
/// activation energies might precompute the terms \f$ exp(-E/RT)
/// \f$, instead of evaluating the exponential repeatedly for each
/// reaction. There are many other possible 'management styles',
/// each of which might be better suited to some reaction
/// mechanisms than others.
///
/// But however a manager structures the internal computation, the
/// tasks the manager class must perform are, for the most part,
/// the same. It must be able to compute reaction rates, species
/// production rates, equilibrium constants, etc. Therefore, all
/// kinetics manager classes should have a common set of public
/// methods, but differ in how they implement these methods.
///
/// A kinetics manager computes reaction rates of progress,
/// species production rates, equilibrium constants, and similar
/// quantities for a reaction mechanism. All kinetics manager
/// classes derive from class Kinetics, which defines a common
/// public interface for all kinetics managers. Each derived class
/// overloads the virtual methods of Kinetics to implement a
/// particular kinetics model.
///
/// For example, class GasKinetics implements reaction rate
/// expressions appropriate for homogeneous reactions in ideal gas
/// mixtures, and class InterfaceKinetics implements expressions
/// appropriate for heterogeneous mechanisms at interfaces,
/// including how to handle reactions involving charged species of
/// phases with different electric potentials --- something that
/// class GasKinetics doesn't deal with at all.
///
/// Many of the methods of class Kinetics write into arrays the
/// values of some quantity for each species, for example the net
/// production rate.  These methods always write the results into
/// flat arrays, ordered by phase in the order the phase was
/// added, and within a phase in the order the species were added
/// to the phase (which is the same ordering as in the input
/// file). Example: suppose a heterogeneous mechanism involves
/// three phases -- a bulk phase 'a', another bulk phase 'b', and
/// the surface phase 'a:b' at the a/b interface. Phase 'a'
/// contains 12 species, phase 'b' contains 3, and at the
/// interface there are 5 adsorbed species defined in phase
/// 'a:b'. Then methods like getNetProductionRates(doublereal* net)
/// will write and output array of length 20, beginning at the location
/// pointed to by 'net'. The first 12 values will be the net production
/// rates for all 12 species of phase 'a' (even if some do not participate
/// in the reactions), the next 3 will be for phase 'b', and finally the
/// net production rates for the surface species will occupy the last
/// 5 locations.
/// @ingroup chemkinetics


//! Public interface for kinetics managers.
/*!
 * This class serves as a base class to derive 'kinetics
 * managers', which are classes that manage homogeneous chemistry
 * within one phase, or heterogeneous chemistry at one
 * interface. The virtual methods of this class are meant to be
 * overloaded in subclasses. The non-virtual methods perform
 * generic functions and are implemented in Kinetics. They should
 * not be overloaded. Only those methods required by a subclass
 * need to be overloaded; the rest will throw exceptions if
 * called.
 *
 * When the nomenclature "kinetics species index" is used below,
 * this means that the species index ranges over all species in
 * all phases handled by the kinetics manager.
 *
 *  @ingroup kineticsmgr
 */
class Kinetics
{

public:
    /**
     * @name Constructors and General Information about Mechanism
     */
    //@{

    /// Default constructor.
    Kinetics();

    /// Destructor.
    virtual ~Kinetics();

    //!Copy Constructor for the Kinetics object.
    Kinetics(const Kinetics&);

    //! Assignment operator
    /*!
     * @param right    Reference to Kinetics object to be copied into the
     *                 current one.
     */
    Kinetics& operator=(const Kinetics& right);

    //! Duplication routine for objects which inherit from Kinetics
    /*!
     *  This function can be used to duplicate objects derived from Kinetics
     *  even if the application only has a pointer to Kinetics to work with.
     *
     *  These routines are basically wrappers around the derived copy
     *  constructor.
     *
     * @param  tpVector Vector of pointers to ThermoPhase objects. this is the
     *                  #m_thermo vector within this object
     */
    virtual Kinetics* duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const;

    //! Reassign the pointers within the Kinetics object
    /*!
     *  This type or routine is necessary because the Kinetics object doesn't
     *  own the ThermoPhase objects. After a duplication, we need to point to
     *  different ThermoPhase objects.
     *
     *  We check that the ThermoPhase objects are aligned in the same order and have
     *  the following identical properties to the ones that they are replacing:
     *
     *  - ThermoPhase::id()
     *  - ThermoPhase::eosType()
     *  - ThermoPhase::nSpecies()
     *
     *  @param tpVector Vector of pointers to ThermoPhase objects. this is the
     *         #m_thermo vector within this object
     */
    virtual void assignShallowPointers(const std::vector<thermo_t*> & tpVector);

    //!  Identifies the kinetics manager type.
    /*!
     *   Each class derived from Kinetics should overload this method to
     *   return a unique integer. Standard values are defined in file
     *   mix_defs.h.
     */
    virtual int type() const;

    //! Number of reactions in the reaction mechanism.
    size_t nReactions() const {
        return m_ii;
    }

    //! Check that the specified reaction index is in range
    //! Throws an exception if i is greater than nReactions()
    void checkReactionIndex(size_t m) const;

    //! Check that an array size is at least nReactions()
    //! Throws an exception if ii is less than nReactions(). Used before calls
    //! which take an array pointer.
    void checkReactionArraySize(size_t ii) const;

    //! Check that the specified species index is in range
    //! Throws an exception if k is greater than nSpecies()-1
    void checkSpeciesIndex(size_t k) const;

    //! Check that an array size is at least nSpecies()
    //! Throws an exception if kk is less than nSpecies(). Used before calls
    //! which take an array pointer.
    void checkSpeciesArraySize(size_t mm) const;

    //@}
    //! @name Information/Lookup Functions about Phases and Species
    //@{

    /**
     * The number of phases participating in the reaction
     * mechanism. For a homogeneous reaction mechanism, this will
     * always return 1, but for a heterogeneous mechanism it will
     * return the total number of phases in the mechanism.
     */
    size_t nPhases() const {
        return m_thermo.size();
    }

    //! Check that the specified phase index is in range
    //! Throws an exception if m is greater than nPhases()
    void checkPhaseIndex(size_t m) const;

    //! Check that an array size is at least nPhases()
    //! Throws an exception if mm is less than nPhases(). Used before calls
    //! which take an array pointer.
    void checkPhaseArraySize(size_t mm) const;

    /**
     * Return the phase index of a phase in the list of phases
     * defined within the object.
     *
     *  @param ph std::string name of the phase
     *
     * If a -1 is returned, then the phase is not defined in
     * the Kinetics object.
     */
    size_t phaseIndex(const std::string& ph) {
        if (m_phaseindex.find(ph) == m_phaseindex.end()) {
            return npos;
        } else {
            return m_phaseindex[ph] - 1;
        }
    }

    /**
     * This returns the integer index of the phase which has ThermoPhase type
     * cSurf. For heterogeneous mechanisms, this identifies the one surface
     * phase. For homogeneous mechanisms, this returns -1.
     */
    size_t surfacePhaseIndex() {
        return m_surfphase;
    }

    /**
     * Phase where the reactions occur. For heterogeneous mechanisms, one of
     * the phases in the list of phases represents the 2D interface or 1D edge
     * at which the reactions take place. This method returns the index of the
     * phase with the smallest spatial dimension (1, 2, or 3) among the list
     * of phases.  If there is more than one, the index of the first one is
     * returned. For homogeneous mechanisms, the value 0 is returned.
     */
    size_t reactionPhaseIndex() {
        return m_rxnphase;
    }

    /**
     * This method returns a reference to the nth ThermoPhase object defined
     * in this kinetics mechanism.  It is typically used so that member
     * functions of the ThermoPhase object may be called. For homogeneous
     * mechanisms, there is only one object, and this method can be called
     * without an argument to access it.
     *
     * @param n Index of the ThermoPhase being sought.
     */
    thermo_t& thermo(size_t n=0) {
        return *m_thermo[n];
    }
    const thermo_t& thermo(size_t n=0) const {
        return *m_thermo[n];
    }

    /**
     * The total number of species in all phases participating in the kinetics
     * mechanism. This is useful to dimension arrays for use in calls to
     * methods that return the species production rates, for example.
     */
    size_t nTotalSpecies() const {
        return m_kk;
    }

    /**
     * The location of species k of phase n in species arrays.
     * Kinetics manager classes return species production rates in
     * flat arrays, with the species of each phases following one
     * another, in the order the phases were added.  This method
     * is useful to find the value for a particular species of a
     * particular phase in arrays returned from methods like
     * getCreationRates that return an array of species-specific
     * quantities.
     *
     * Example: suppose a heterogeneous mechanism involves three
     * phases.  The first contains 12 species, the second 26, and
     * the third 3.  Then species arrays must have size at least
     * 41, and positions 0 - 11 are the values for the species in
     * the first phase, positions 12 - 37 are the values for the
     * species in the second phase, etc.  Then
     * kineticsSpeciesIndex(7, 0) = 7, kineticsSpeciesIndex(4, 1)
     * = 16, and kineticsSpeciesIndex(2, 2) = 40.
     *
     * @param k species index
     * @param n phase index for the species
     */
    size_t kineticsSpeciesIndex(size_t k, size_t n) const {
        return m_start[n] + k;
    }

    //! Return the name of the kth species in the kinetics manager.
    /*!
     *  k is an integer from 0 to ktot - 1, where ktot is the number of
     * species in the kinetics manager, which is the sum of the number of
     * species in all phases participating in the kinetics manager.  If k is
     * out of bounds, the string "<unknown>" is returned.
     *
     * @param k species index
     */
    std::string kineticsSpeciesName(size_t k) const;

    /**
     * This routine will look up a species number based on the input
     * std::string nm. The lookup of species will occur for all phases
     * listed in the kinetics object.
     *
     *  return
     *   - If a match is found, the position in the species list is returned.
     *   - If no match is found, the value -1 is returned.
     *
     * @param nm   Input string name of the species
     */
    size_t kineticsSpeciesIndex(const std::string& nm) const;

    /**
     * This routine will look up a species number based on the input
     * std::string nm. The lookup of species will occur in the specified
     * phase of the object, or all phases if ph is "<any>".
     *
     *  return
     *   - If a match is found, the position in the species list is returned.
     *   - If no match is found, the value npos (-1) is returned.
     *
     * @param nm   Input string name of the species
     * @param ph   Input string name of the phase.
     */
    size_t kineticsSpeciesIndex(const std::string& nm,
                                const std::string& ph) const;

    /**
     * This function looks up the name of a species and returns a
     * reference to the ThermoPhase object of the phase where the species
     * resides. Will throw an error if the species doesn't match.
     *
     * @param nm   String containing the name of the species.
     */
    thermo_t& speciesPhase(const std::string& nm);

    /**
     * This function takes as an argument the kineticsSpecies index
     * (i.e., the list index in the list of species in the kinetics
     * manager) and returns the species' owning ThermoPhase object.
     *
     * @param k          Species index
     */
    thermo_t& speciesPhase(size_t k) {
        return thermo(speciesPhaseIndex(k));
    }

    /**
     * This function takes as an argument the kineticsSpecies index (i.e., the
     * list index in the list of species in the kinetics manager) and returns
     * the index of the phase owning the species.
     *
     * @param k          Species index
     */
    size_t speciesPhaseIndex(size_t k);

    //! @}
    //! @name Reaction Rates Of Progress
    //! @{

    //!  Return the forward rates of progress of the reactions
    /*!
     * Forward rates of progress.  Return the forward rates of
     * progress in array fwdROP, which must be dimensioned at
     * least as large as the total number of reactions.
     *
     * @param fwdROP  Output vector containing forward rates
     *                of progress of the reactions. Length: m_ii.
     */
    virtual void getFwdRatesOfProgress(doublereal* fwdROP) {
        throw NotImplementedError("Kinetics::getFwdRatesOfProgress");
    }

    //!  Return the Reverse rates of progress of the reactions
    /*!
     * Return the reverse rates of progress in array revROP, which must be
     * dimensioned at least as large as the total number of reactions.
     *
     * @param revROP  Output vector containing reverse rates
     *                of progress of the reactions. Length: m_ii.
     */
    virtual void getRevRatesOfProgress(doublereal* revROP) {
        throw NotImplementedError("Kinetics::getRevRatesOfProgress");
    }

    /**
     * Net rates of progress.  Return the net (forward - reverse) rates of
     * progress in array netROP, which must be dimensioned at least as large
     * as the total number of reactions.
     *
     * @param netROP  Output vector of the net ROP. Length: m_ii.
     */
    virtual void getNetRatesOfProgress(doublereal* netROP) {
        throw NotImplementedError("Kinetics::getNetRatesOfProgress");
    }

    //! Return a vector of Equilibrium constants.
    /*!
     *  Return the equilibrium constants of the reactions in concentration
     *  units in array kc, which must be dimensioned at least as large as the
     *  total number of reactions.
     *
     * @param kc   Output vector containing the equilibrium constants.
     *             Length: m_ii.
     */
    virtual void getEquilibriumConstants(doublereal* kc) {
        throw NotImplementedError("Kinetics::getEquilibriumConstants");
    }

    /**
     * Change in species properties. Given an array of molar species
     * property values \f$ z_k, k = 1, \dots, K \f$, return the
     * array of reaction values
     * \f[
     * \Delta Z_i = \sum_k \nu_{k,i} z_k, i = 1, \dots, I.
     * \f]
     * For example, if this method is called with the array of
     * standard-state molar Gibbs free energies for the species,
     * then the values returned in array \c deltaProperty would be
     * the standard-state Gibbs free energies of reaction for each
     * reaction.
     *
     * @param property Input vector of property value. Length: m_kk.
     * @param deltaProperty Output vector of deltaRxn. Length: m_ii.
     */
    virtual void getReactionDelta(const doublereal* property,
                                  doublereal* deltaProperty) {
        throw NotImplementedError("Kinetics::getReactionDelta");
    }

    //! Return the vector of values for the reaction gibbs free energy change.
    /*!
     * These values depend upon the concentration of the solution.
     *
     *  units = J kmol-1
     *
     * @param deltaG  Output vector of  deltaG's for reactions Length: m_ii.
     */
    virtual void getDeltaGibbs(doublereal* deltaG) {
        throw NotImplementedError("Kinetics::getDeltaGibbs");
    }

    //! Return the vector of values for the reaction electrochemical free
    //! energy change.
    /*!
     * These values depend upon the concentration of the solution and the
     * voltage of the phases
     *
     *  units = J kmol-1
     *
     * @param deltaM  Output vector of  deltaM's for reactions Length: m_ii.
     */
    virtual void getDeltaElectrochemPotentials(doublereal* deltaM) {
        throw NotImplementedError("Kinetics::getDeltaElectrochemPotentials");
    }

    /**
     * Return the vector of values for the reactions change in enthalpy.
     * These values depend upon the concentration of the solution.
     *
     *  units = J kmol-1
     *
     * @param deltaH  Output vector of deltaH's for reactions Length: m_ii.
     */
    virtual void getDeltaEnthalpy(doublereal* deltaH) {
        throw NotImplementedError("Kinetics::getDeltaEnthalpy");
    }

    /**
     * Return the vector of values for the reactions change in entropy.  These
     * values depend upon the concentration of the solution.
     *
     *  units = J kmol-1 Kelvin-1
     *
     * @param deltaS  Output vector of deltaS's for reactions Length: m_ii.
     */
    virtual void getDeltaEntropy(doublereal* deltaS) {
        throw NotImplementedError("Kinetics::getDeltaEntropy");
    }

    /**
     * Return the vector of values for the reaction standard state
     * gibbs free energy change.  These values don't depend upon
     * the concentration of the solution.
     *
     *  units = J kmol-1
     *
     * @param deltaG  Output vector of ss deltaG's for reactions Length: m_ii.
     */
    virtual void getDeltaSSGibbs(doublereal* deltaG) {
        throw NotImplementedError("Kinetics::getDeltaSSGibbs");
    }

    /**
     * Return the vector of values for the change in the standard
     * state enthalpies of reaction.  These values don't depend
     * upon the concentration of the solution.
     *
     *  units = J kmol-1
     *
     * @param deltaH  Output vector of ss deltaH's for reactions Length: m_ii.
     */
    virtual void getDeltaSSEnthalpy(doublereal* deltaH) {
        throw NotImplementedError("Kinetics::getDeltaSSEnthalpy");
    }

    /**
     * Return the vector of values for the change in the standard
     * state entropies for each reaction.  These values don't
     * depend upon the concentration of the solution.
     *
     *  units = J kmol-1 Kelvin-1
     *
     * @param deltaS  Output vector of ss deltaS's for reactions Length: m_ii.
     */
    virtual void getDeltaSSEntropy(doublereal* deltaS) {
        throw NotImplementedError("Kinetics::getDeltaSSEntropy");
    }

    //! @}
    //! @name Species Production Rates
    //! @{

    /**
     * Species creation rates [kmol/m^3/s or kmol/m^2/s]. Return the species
     * creation rates in array cdot, which must be dimensioned at least as
     * large as the total number of species in all phases. @see nTotalSpecies.
     *
     * @param cdot   Output vector of creation rates. Length: m_kk.
     */
    virtual void getCreationRates(doublereal* cdot) {
        throw NotImplementedError("Kinetics::getCreationRates");
    }

    /**
     * Species destruction rates [kmol/m^3/s or kmol/m^2/s]. Return the
     * species destruction rates in array ddot, which must be dimensioned at
     * least as large as the total number of species. @see nTotalSpecies.
     *
     * @param ddot   Output vector of destruction rates. Length: m_kk.
     */
    virtual void getDestructionRates(doublereal* ddot) {
        throw NotImplementedError("Kinetics::getDestructionRates");
    }

    /**
     * Species net production rates [kmol/m^3/s or kmol/m^2/s]. Return
     * the species net production rates (creation - destruction)
     * in array wdot, which must be dimensioned at least as large
     * as the total number of species. @see nTotalSpecies.
     *
     * @param wdot   Output vector of net production rates. Length: m_kk.
     */
    virtual void getNetProductionRates(doublereal* wdot) {
        throw NotImplementedError("Kinetics::getNetProductionRates");
    }

    //! @}
    //! @name Reaction Mechanism Informational Query Routines
    //! @{

    /**
     * Stoichiometric coefficient of species k as a reactant in reaction i.
     *
     * @param k   kinetic species index
     * @param i   reaction index
     */
    virtual doublereal reactantStoichCoeff(size_t k, size_t i) const {
        throw NotImplementedError("Kinetics::reactantStoichCoeff");
    }

    /**
     * Stoichiometric coefficient of species k as a product in reaction i.
     *
     * @param k   kinetic species index
     * @param i   reaction index
     */
    virtual doublereal productStoichCoeff(size_t k, size_t i) const {
        throw NotImplementedError("Kinetics::productStoichCoeff");
    }

    //! Reactant order of species k in reaction i.
    /*!
     * This is the nominal order of the activity concentration in
     * determining the forward rate of progress of the reaction
     *
     * @param k   kinetic species index
     * @param i   reaction index
     */
    virtual doublereal reactantOrder(size_t k, size_t i) const {
        throw NotImplementedError("Kinetics::reactantOrder");
    }

    //! product Order of species k in reaction i.
    /*!
     * This is the nominal order of the activity concentration of species k in
     * determining the reverse rate of progress of the reaction i
     *
     * For irreversible reactions, this will all be zero.
     *
     * @param k   kinetic species index
     * @param i   reaction index
     */
    virtual doublereal productOrder(int k, int i) const {
        throw NotImplementedError("Kinetics::productOrder");
    }

    //! Get the vector of activity concentrations used in the kinetics object
    /*!
     *  @param[out] conc  Vector of activity concentrations. Length is equal
     *               to the number of species in the kinetics object
     */
    virtual void getActivityConcentrations(doublereal* const conc) {
        throw NotImplementedError("Kinetics::getActivityConcentrations");
    }

    /**
     * Returns a read-only reference to the vector of reactant
     * index numbers for reaction i.
     *
     * @param i  reaction index
     */
    virtual const std::vector<size_t>& reactants(size_t i) const {
        return m_reactants[i];
    }

    /**
     * Returns a read-only reference to the vector of product
     * index numbers for reaction i.
     *
     * @param i reaction index
     */
    virtual const std::vector<size_t>& products(size_t i) const {
        return m_products[i];
    }

    /**
     * Flag specifying the type of reaction. The legal values and
     * their meaning are specific to the particular kinetics
     * manager.
     *
     * @param i   reaction index
     */
    virtual int reactionType(size_t i) const {
        throw NotImplementedError("Kinetics::reactionType");
    }

    /**
     * True if reaction i has been declared to be reversible. If
     * isReversible(i) is false, then the reverse rate of progress
     * for reaction i is always zero.
     *
     * @param i   reaction index
     */
    virtual bool isReversible(size_t i) {
        throw NotImplementedError("Kinetics::isReversible");
    }

    /**
     * Return a string representing the reaction.
     *
     * @param i   reaction index
     */
    virtual std::string reactionString(size_t i) const {
        throw NotImplementedError("Kinetics::reactionStd::String");
    }

    //! Returns a string containing the reactants side of the reaction equation.
    virtual std::string reactantString(size_t i) const {
        throw NotImplementedError("Kinetics::reactionString");
    }

    //! Returns a string containing the products side of the reaction equation.
    virtual std::string productString(size_t i) const {
        throw NotImplementedError("Kinetics::productString");
    }

    /**
     * Return the forward rate constants
     *
     * length is the number of reactions. units depends on many issues.
     *
     * @param kfwd    Output vector containing the forward reaction rate
     *                constants. Length: m_ii.
     */
    virtual void getFwdRateConstants(doublereal* kfwd) {
        throw NotImplementedError("Kinetics::getFwdRateConstants");
    }

    /**
     * Return the reverse rate constants.
     *
     * length is the number of reactions. units depends on many issues. Note,
     * this routine will return rate constants for irreversible reactions if
     * the default for doIrreversible is overridden.
     *
     * @param krev   Output vector of reverse rate constants.
     * @param doIrreversible boolean indicating whether irreversible reactions
     *                       should be included.
     */
    virtual void getRevRateConstants(doublereal* krev,
                                     bool doIrreversible = false) {
        throw NotImplementedError("Kinetics::getFwdRateConstants");
    }

    //! @}
    //! @name Reaction Mechanism Construction
    //! @{

    //!  Add a phase to the kinetics manager object.
    /*!
     * This must be done before the function init() is called or before any
     * reactions are input. The following fields are updated:
     *
     *  - #m_start -> vector of integers, containing the starting position of
     *    the species for each phase in the kinetics mechanism.
     *  - #m_surfphase -> index of the surface phase.
     *  - #m_thermo -> vector of pointers to ThermoPhase phases that
     *    participate in the kinetics mechanism.
     *  - #m_phaseindex -> map containing the std::string id of each
     *    ThermoPhase phase as a key and the index of the phase within the
     *    kinetics manager object as the value.
     *
     * @param thermo    Reference to the ThermoPhase to be added.
     */
    virtual void addPhase(thermo_t& thermo);

    /**
     * Prepare the class for the addition of reactions. This method is called
     * by importKinetics() after all phases have been added but before any
     * reactions have been. The base class method does nothing, but derived
     * classes may use this to perform any initialization (allocating arrays,
     * etc.) that requires knowing the phases and species, but before any
     * reactions are added.
     */
    virtual void init() {}

    /**
     * Finish adding reactions and prepare for use. This method is called by
     * importKinetics() after all reactions have been entered into the
     * mechanism and before the mechanism is used to calculate reaction rates.
     * The base class method does nothing, but derived classes may use this to
     * perform any initialization (allocating arrays, etc.) that must be done
     * after the reactions are entered.
     */
    virtual void finalize();

    /**
     * Add a single reaction to the mechanism. This routine
     * must be called after init() and before finalize().
     *
     * @param r      Reference to the ReactionData object for the reaction
     *               to be added.
     */
    virtual void addReaction(ReactionData& r) {
        throw NotImplementedError("Kinetics::addReaction");
    }

    virtual const std::vector<grouplist_t>& reactantGroups(size_t i) {
        return m_dummygroups;
    }

    virtual const std::vector<grouplist_t>& productGroups(size_t i) {
        return m_dummygroups;
    }

    //@}
    //! @name Altering Reaction Rates
    /*!
     * These methods alter reaction rates. They are designed primarily for
     * carrying out sensitivity analysis, but may be used for any purpose
     * requiring dynamic alteration of rate constants.  For each reaction, a
     * real-valued multiplier may be defined that multiplies the reaction rate
     * coefficient. The multiplier may be set to zero to completely remove a
     * reaction from the mechanism.
     */
    //@{

    //! The current value of the multiplier for reaction i.
    /*!
     * @param i index of the reaction
     */
    doublereal multiplier(size_t i) const {
        return m_perturb[i];
    }

    //! Set the multiplier for reaction i to f.
    /*!
     *  @param i  index of the reaction
     *  @param f  value of the multiplier.
     */
    void setMultiplier(size_t i, doublereal f) {
        m_perturb[i] = f;
    }

    //@}

    /**
     * Increment the number of reactions in the mechanism by one.
     * @todo Should be protected?
     */
    void incrementRxnCount() {
        m_ii++;
        m_perturb.push_back(1.0);
    }

    /**
     * Returns true if the kinetics manager has been properly
     * initialized and finalized.
     */
    virtual bool ready() const {
        return false;
    }

    /*!
     * Takes as input an array of properties for all species in the mechanism
     * and copies those values belonging to a particular phase to the output
     * array.
     * @param data Input data array.
     * @param phase Pointer to one of the phase objects participating in this
     *     reaction mechanism
     * @param phase_data Output array where the values for the the specified
     *     phase are to be written.
     */
    void selectPhase(const doublereal* data, const thermo_t* phase,
                     doublereal* phase_data);

protected:
    //! Number of reactions in the mechanism
    size_t m_ii;

    //! The number of species in all of the phases
    //! that participate in this kinetics mechanism.
    size_t m_kk;

    /// Vector of perturbation factors for each reaction's rate of
    /// progress vector. It is initialized to one.
    vector_fp m_perturb;

    /**
     * This is a vector of vectors containing the reactants for
     * each reaction. The outer vector is over the number of
     * reactions, m_ii.  The inner vector is a list of species
     * indices. If the stoichiometric coefficient for a reactant
     * is greater than one, then the reactant is listed
     * contiguously in the vector a number of times equal to its
     * stoichiometric coefficient.
     * NOTE: These vectors will be wrong if there are real
     *       stoichiometric coefficients in the expression.
     */
    std::vector<std::vector<size_t> > m_reactants;

    /**
     * This is a vector of vectors containing the products for
     * each reaction. The outer vector is over the number of
     * reactions, m_ii.  The inner vector is a list of species
     * indices. If the stoichiometric coefficient for a product is
     * greater than one, then the reactant is listed contiguously
     * in the vector a number of times equal to its stoichiometric
     * coefficient.
     * NOTE: These vectors will be wrong if there are real
     *       stoichiometric coefficients in the expression.
     */
    std::vector<std::vector<size_t> > m_products;

    //! m_thermo is a vector of pointers to ThermoPhase objects that are
    //! involved with this kinetics operator
    /*!
     * For homogeneous kinetics applications, this vector
     * will only have one entry. For interfacial reactions, this
     * vector will consist of multiple entries; some of them will
     * be surface phases, and the other ones will be bulk phases.
     * The order that the objects are listed determines the order
     * in which the species comprising each phase are listed in
     * the source term vector, originating from the reaction
     * mechanism.
     *
     * Note that this kinetics object doesn't own these ThermoPhase objects
     * and is not responsible for creating or deleting them.
     */
    std::vector<thermo_t*> m_thermo;

    /**
     * m_start is a vector of integers specifying the beginning position
     * for the species vector for the n'th phase in the kinetics
     * class.
     */
    std::vector<size_t>  m_start;

    /**
     * Mapping of the phase id, i.e., the id attribute in the xml
     * phase element to the position of the phase within the
     * kinetics object.  Positions start with the value of 1. The
     * member function, phaseIndex() decrements by one before
     * returning the index value, so that missing phases return
     * -1.
     */
    std::map<std::string, size_t> m_phaseindex;

    //! Index in the list of phases of the one surface phase.
    size_t m_surfphase;

    //! Phase Index where reactions are assumed to be taking place
    /*!
     *  We calculate this by assuming that the phase with the lowest
     *  dimensionality is the phase where reactions are taking place.
     */
    size_t m_rxnphase;

    //! number of spatial dimensions of lowest-dimensional phase.
    size_t m_mindim;

private:
    //! Vector of group lists
    std::vector<grouplist_t> m_dummygroups;
};

}

#endif
