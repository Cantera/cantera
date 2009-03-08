/**
 * @file Kinetics.h
 *
 *  $Author: hkmoffa $
 *  $Date: 2006/04/30 18:01:43 $
 *  $Revision: 1.24 $
 */

// Copyright 2001-2004  California Institute of Technology

/**
 * @defgroup kineticsGroup Kinetics
 */

#ifndef CT_KINETICS_H
#define CT_KINETICS_H

#include "ctexceptions.h"
#include "ThermoPhase.h"
#include "mix_defs.h"

namespace Cantera {

    // forward references
    class ReactionData;

    /// @defgroup kineticsmgr Kinetics Managers 
    /// @section kinmodman Models and Managers
    /// 
    /// A kinetics manager is a C++ class that implements a kinetics
    /// model; a kinetics model is a set of mathematical equation
    /// describing how various kinetic quanities are to be computed --
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
    /// Kinetics managers may be also created that hard-wire a
    /// particular reaction mechanism in C++ code. This can often
    /// result in faster performance. An example of this is the
    /// kinetics manager GRI30_Kinetics that hard-wires the rate
    /// expressions for the natural gas combustion mechanism GRI-3.0.
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


    /// Public interface for kinetics managers. This class serves as a
    /// base class to derive 'kinetics managers', which are classes
    /// that manage homogeneous chemistry within one phase, or
    /// heterogeneous chemistry at one interface. The virtual methods
    /// of this class are meant to be overloaded in subclasses. The
    /// non-virtual methods perform generic functions and are
    /// implemented in Kinetics. They should not be overloaded. Only
    /// those methods required by a subclass need to be overloaded;
    /// the rest will throw exceptions if called.  @ingroup kinetics
    /// @ingroup kineticsmgr

    class Kinetics {

    public:

        // typedefs
        typedef ThermoPhase thermo_t;

	/**
	 * @name Constructors and General Information about Mechanism
	 */
	//@{

        /// Default constructor.
        Kinetics();

        /// This constructor initializes with a starting phase.
        /// @deprecated
        //        Kinetics(thermo_t* thermo);

        /// Destructor. 
        virtual ~Kinetics();

        ///  Identifies the kinetics manager type.  Each class derived
        ///  from Kinetics should overload this method to return a
        ///  unique integer. Standard values are defined in file
        ///  mix_defs.h.
        virtual int type() { return 0; }

        /// Number of reactions in the reaction mechanism.
        int nReactions() const {return m_ii;}

	//@}


        /**
         * @name Information/Lookup Functions about Phases and Species
         */
        //@{

        /**
         * The number of phases participating in the reaction
         * mechanism. For a homogeneous reaction mechanism, this will
         * always return 1, but for a heterogeneous mechanism it will
         * return the total number of phases in the mechanism.
         */  
        int nPhases() const { return static_cast<int>(m_thermo.size()); }

	/**
	 * Return the phase index of a phase in the list of phases
	 * defined within the object.
	 *
	 *  @param ph string name of the phase
	 *
	 * If a -1 is returned, then the phase is not defined in
	 * the Kinetics object.
	 */
        int phaseIndex(string ph) { 
            if (m_phaseindex.find(ph) == m_phaseindex.end()) {
                return -1;
            }
            else {
                return m_phaseindex[ph] - 1;
            }
        }

	/**
	 * This returns the integer index of the phase which has
	 * ThermoPhase type cSurf. For heterogeneous mechanisms, this
	 * identifies the one surface phase. For homogeneous
	 * mechanisms, this reurns -1.
	 */
        int surfacePhaseIndex() { return m_surfphase; }

        /**
         * Phase where the reactions occur. For heterogeneous
         * mechanisms, one of the phases in the list of phases
         * represents the 2D interface or 1D edge at which the
         * reactions take place. This method returns the index of the
         * phase with the smallest spatial dimension (1, 2, or 3)
         * among the list of phases.  If there is more than one, the
         * index of the first one is returned. For homogeneous
         * mechanisms, the value 0 is returned.
         */
        int reactionPhaseIndex() { return m_rxnphase; }


	/**
	 * This method returns a reference to the nth ThermoPhase
	 * object defined in this kinetics mechanism.  It is typically
	 * used so that member functions of the ThermoPhase object may
	 * be called. For homogeneous mechanisms, there is only one
	 * object, and this method can be called without an argument
	 * to access it.
	 */
        thermo_t& thermo(int n=0) { return *m_thermo[n]; }
        const thermo_t& thermo(int n=0) const { return *m_thermo[n]; }

	/**
	 * This method returns a reference to the nth ThermoPhase
	 * defined in this kinetics mechanism.
	 * It is typically used so that member functions of the
	 * ThermoPhase may be called. @deprecated This method is redundant.
	 */
        thermo_t& phase(int n=0) { 
            deprecatedMethod("Kinetics","phase","thermo");
            return *m_thermo[n]; 
        }
        const thermo_t& phase(int n=0) const { 
            deprecatedMethod("Kinetics","phase","thermo");
            return *m_thermo[n]; 
        }

        /**
	 * The total number of species in all phases participating in
	 * the kinetics mechanism. This is useful to dimension arrays
	 * for use in calls to methods that return the species
	 * production rates, for example.
	 */
        int nTotalSpecies() const {
            int n=0, np;
            np = nPhases();
            for (int p = 0; p < np; p++) n += thermo(p).nSpecies();
            return n;
        }

	/**
	 * Returns the starting index of the species in the nth phase
	 * associated with the reaction mechanism.
	 *
	 * @param n Return the index of first species in the nth phase
	 *          associated with the reaction mechanism.
	 */
        int start(int n) { 
            deprecatedMethod("Kinetics","start","kineticsSpeciesIndex(0,n)");
            return m_start[n]; 
        }


	/**
	 * The location of species k of phase n in species arrays.
         * Kinetics manager classes return species production rates in
         * flat arrays, with the species of each phases following one
         * another, in the order the phases were added.  This method
         * is useful to find the value for a particular species of a
         * particular phase in arrrays returned from methods like
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
        int kineticsSpeciesIndex(int k, int n) const {
            return m_start[n] + k;
        }

	/**
	 * Return the string name of the kth species in the kinetics
	 * manager. k is an integer from 0 to ktot - 1, where ktot is
	 * the number of species in the kinetics manager, which is the
	 * sum of the number of species in all phases participating in
	 * the kinetics manager.  If k is out of bounds, the string
	 * "<unknown>" is returned.
	 */
        string kineticsSpeciesName(int k) const;

	/**
	 * This routine will look up a species number based on
	 * the input string nm. The lookup of species will
	 * occur for all phases listed in the kinetics object,
	 * unless the string ph refers to a specific phase of
	 * the object. 
	 *
	 *  return
	 *   - If a match is found, the position in the species list
	 *   is returned. 
         *   - If a specific phase is specified and no match is found,
         *   the value -1 is returned.
	 *   - If no match is found in any phase, the value -2 is returned.
	 */
        int kineticsSpeciesIndex(string nm, string ph = "<any>") const;

	/**
	 * This function looks up the string name of a species and
	 * returns a reference to the ThermoPhase object of the
	 * phase where the species resides.
	 * Will throw an error if the species string doesn't match.
	 */
        thermo_t& speciesPhase(string nm);

	/**
	 * This function takes as an argument the kineticsSpecies index
	 * (i.e., the list index in the list of species in the kinetics
	 * manager) and returns the species' owning ThermoPhase object.
	 */
        thermo_t& speciesPhase(int k) {
            return thermo(speciesPhaseIndex(k));                
        }

	/**
	 * This function takes as an argument the kineticsSpecies index
	 * (i.e., the list index in the list of species in the kinetics
	 * manager) and returns the index of the phase owning the 
	 * species.
	 */
        int speciesPhaseIndex(int k);

	//@}



        /**
         * @name Reaction Rates Of Progress
         */
        //@{

        /**
         * Forward rates of progress.  Return the forward rates of
         * progress in array fwdROP, which must be dimensioned at
         * least as large as the total number of reactions.
         */
        virtual void getFwdRatesOfProgress(doublereal* fwdROP) {
            err("getFwdRatesOfProgress");
        }
 
        /**
         * Reverse rates of progress.  Return the reverse rates of
         * progress in array revROP, which must be dimensioned at
         * least as large as the total number of reactions.
         */
        virtual void getRevRatesOfProgress(doublereal* revROP) {
            err("getRevRatesOfProgress");
        } 

        /**
         * Net rates of progress.  Return the net (forward - reverse)
         * rates of progress in array netROP, which must be
         * dimensioned at least as large as the total number of
         * reactions.
         */
        virtual void getNetRatesOfProgress(doublereal* netROP) {
            err("getNetRatesOfProgress");
        }


        /**
         * Equilibrium constants. Return the equilibrium constants of
         * the reactions in concentration units in array kc, which
         * must be dimensioned at least as large as the total number
         * of reactions. 
         */  
         virtual void getEquilibriumConstants(doublereal* kc) {
             err("getEquilibriumConstants");
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
         */
        virtual void getReactionDelta(const doublereal* property,
            doublereal* deltaProperty) {
            err("getReactionDelta");
        }

	/**
	 * Return the vector of values for the reaction gibbs free
	 * energy change.  These values depend upon the concentration
	 * of the solution.
	 *
	 *  units = J kmol-1
	 */
	virtual void getDeltaGibbs( doublereal* deltaG) {
	    err("getDeltaGibbs");
	}

	/**
	 * Return the vector of values for the reactions change in
	 * enthalpy.  These values depend upon the concentration of
	 * the solution.
	 *
	 *  units = J kmol-1
	 */
	virtual void getDeltaEnthalpy( doublereal* deltaH) {
	    err("getDeltaEnthalpy");
	}

	/**
	 * Return the vector of values for the reactions change in
	 * entropy.  These values depend upon the concentration of the
	 * solution.
	 *
	 *  units = J kmol-1 Kelvin-1
	 */
	virtual void getDeltaEntropy( doublereal* deltaS) {
	    err("getDeltaEntropy");
	}

	/**
	 * Return the vector of values for the reaction standard state
	 * gibbs free energy change.  These values don't depend upon
	 * the concentration of the solution.
	 *
	 *  units = J kmol-1
	 */
	virtual void getDeltaSSGibbs( doublereal* deltaG) {
	    err("getDeltaSSGibbs");
	}

	/**
	 * Return the vector of values for the change in the standard
	 * state enthalpies of reaction.  These values don't depend
	 * upon the concentration of the solution.
	 *
	 *  units = J kmol-1
	 */
	virtual void getDeltaSSEnthalpy( doublereal* deltaH) {
	    err("getDeltaSSEnthalpy");
	}

	/**
	 * Return the vector of values for the change in the standard
	 * state entropies for each reaction.  These values don't
	 * depend upon the concentration of the solution.
	 *
	 *  units = J kmol-1 Kelvin-1
	 */
	virtual void getDeltaSSEntropy( doublereal* deltaS) {
	    err("getDeltaSSEntropy");
	}


	//@}
        /**
         * @name Species Production Rates
         */
        //@{

        /**
         * Species creation rates [kmol/m^3/s or kmol/m^2/s]. Return the
         * species creation rates in array cdot, which must be
         * dimensioned at least as large as the total number of
         * species in all phases. @see nTotalSpecies.
         *  
         */ 
        virtual void getCreationRates(doublereal* cdot) {
            err("getCreationRates");
        }

        /**
         * Species destruction rates [kmol/m^3/s or kmol/m^2/s]. Return
         * the species destruction rates in array ddot, which must be
         * dimensioned at least as large as the total number of
         * species. @see nTotalSpecies.
         *  
         */ 
        virtual void getDestructionRates(doublereal* ddot) {
            err("getDestructionRates");
        }

        /**
         * Species net production rates [kmol/m^3/s or kmol/m^2/s]. Return
         * the species net production rates (creation - destruction)
         * in array wdot, which must be dimensioned at least as large
         * as the total number of species. @see nTotalSpecies.
         */ 
        virtual void getNetProductionRates(doublereal* wdot) {
            err("getNetProductionRates");
        }

        //@}


        /**
         * @name Reaction Mechanism Informational Query Routines
         */
        //@{

        /**
         * Stoichiometric coefficient of species k as a reactant in
         * reaction i.  
         */
        virtual doublereal reactantStoichCoeff(int k, int i) const { 
            err("reactantStoichCoeff");
            return -1.0;
        }
        /**
         * Stoichiometric coefficient of species k as a product in
         * reaction i.  
         */
        virtual doublereal productStoichCoeff(int k, int i) const { 
            err("productStoichCoeff");
            return -1.0;
        }

        virtual doublereal reactantOrder(int k, int i) const {
            err("reactantOrder");
            return -1.0;
        }

        /**
         * Returns a read-only reference to the vector of reactant
         * index numbers for reaction i.
         */ 
        virtual const vector_int& reactants(int i) const { 
            return m_reactants[i]; 
        }

        /**
         * Returns a read-only reference to the vector of product
         * index numbers for reaction i.
         */ 
        virtual const vector_int& products(int i) const { 
            return m_products[i]; 
        }

        /**
         * Flag specifying the type of reaction. The legal values and
         * their meaning are specific to the particular kinetics
         * manager.
         */
        virtual int reactionType(int i) const { 
            err("reactionType"); 
            return -1; 
        }

	/**
         * True if reaction i has been declared to be reversible. If
         * isReversible(i) is false, then the reverse rate of progress
         * for reaction i is always zero.
         */
        virtual bool isReversible(int i){
	    err("isReversible");
	    return false;
	}

	/**
	 * Return a string representing the reaction.
	 */
        virtual string reactionString(int i) const {
            err("reactionString"); return "<null>";
        }

	/**
	 * Return the forward rate constants
	 *
	 * length is the number of reactions. units depends
	 * on many issues. @todo DGG: recommend changing name to 
         * getFwdRateCoefficients.
	 */
	virtual void getFwdRateConstants(doublereal *kfwd) {
         err("getFwdRateConstants");
	}

	/**
	 * Return the reverse rate constants.
	 *
	 * length is the number of reactions. units depends
	 * on many issues. Note, this routine will return rate constants
	 * for irreversible reactions if the default for
	 * doIrreversible is overridden. @todo DGG: recommend changing name to 
         * getRevRateCoefficients.
	 */
	virtual void getRevRateConstants(doublereal *krev, 
					 bool doIrreversible = false) {
	    err("getFwdRateConstants");
	}


	/**
	 * Return the activation energies in Kelvin.
	 *
	 * length is the number of reactions
	 */
	virtual void getActivationEnergies(doublereal *E) {
            err("getActivationEnergies");
	}


	//@}
        /**
         * @name Reaction Mechanism Construction
         */
        //@{

        /**
         * Add a phase to the kinetics manager object. This must
	 * be done before the function init() is called or 
	 * before any reactions are input.
	 * The following fields are updated:
	 *  m_start -> vector of integers, containing the
	 *             starting position of the species for
	 *             each phase in the kinetics mechanism.
	 *  m_surfphase -> index of the surface phase.
	 *  m_thermo -> vector of pointers to ThermoPhase phases
	 *              that participate in the kinetics 
	 *              mechanism.
	 *  m_phaseindex -> map containing the string id of each
	 *              ThermoPhase phase as a key and the
	 *              index of the phase within the kinetics
	 *              manager object as the value.
         */
        void addPhase(thermo_t& thermo);

        /**
         * Prepare the class for the addition of reactions. This
	 * method is called by function importKinetics after all
	 * phases have been added but before any reactions have
	 * been. The base class method does nothing, but derived
	 * classes may use this to perform any initialization
	 * (allocating arrays, etc.) that requires knowing the phases
	 * and species, but before any reactions are added.
         */
        virtual void init() {}

        /**
	 * Finish adding reactions and prepare for use. This method is
	 * called by function importKinetics after all reactions have
	 * been entered into the mechanism and before the mechanism is
	 * used to calculate reaction rates. The base class method
	 * does nothing, but derived classes may use this to perform
	 * any initialization (allocating arrays, etc.) that must be
	 * done after the reactions are entered.
	 */
        virtual void finalize() {}

	/**
	 * Add a single reaction to the mechanism. This routine
	 * must be called after init() and before finalize().
	 */
        virtual void addReaction(const ReactionData& r) {
	    err("addReaction");
	}

        virtual const vector<grouplist_t>& reactantGroups(int i) { 
	    //err("reactantGroups"); 
	    return m_dummygroups;
	}

        virtual const vector<grouplist_t>& productGroups(int i) {
	    //err("productGroups"); 
	    return m_dummygroups;
	}


	//@}
        /**
         * @name Altering Reaction Rates
         *
         * These methods alter reaction rates. They are designed
         * primarily for carrying out sensitivity analysis, but may be
         * used for any purpose requiring dynamic alteration of rate
         * constants.  For each reaction, a real-valued multiplier may
         * be defined that multiplies the reaction rate
         * coefficient. The multiplier may be set to zero to
         * completely remove a reaction from the mechanism.
         */
        //@{

        /// The current value of the multiplier for reaction i.
        doublereal multiplier(int i) const {return m_perturb[i];}

        /// Set the multiplier for reaction i to f.
        void setMultiplier(int i, doublereal f) {m_perturb[i] = f;}
        
        //@}

	/**
	 * Increment the number of reactions in the mechanism by one.
         * @todo Should be protected?
	 */
        void incrementRxnCount() { m_ii++; m_perturb.push_back(1.0); }

	/**
	 * Returns true if the kinetics manager has been properly
	 * initialized and finalized. 
	 */
        virtual bool ready() const {
	    return false;
	}


        /** 
         * Extract from array \c data the portion pertaining to phase \c phase.
         */
        void selectPhase(const doublereal* data, const thermo_t* phase,
            doublereal* phase_data);


        /// For internal use. May be removed in a future release.
        int index(){ return m_index; }
        void setIndex(int index) { m_index = index; }


    protected:

	
        /// Number of reactions in the mechanism
        int m_ii;
	
        /// Vector of perturbation factors for each reaction's rate of
        /// progress vector. It is initialized to one.
        ///
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
        vector<vector_int> m_reactants;

	/**
	 * This is a vector of vectors containing the products for
	 * each reaction. The outer vector is over the number of
	 * reactions, m_ii.  The inner vector is a list of species
	 * indeces. If the stoichiometric coefficient for a product is
	 * greater than one, then the reactant is listed contiguously
	 * in the vector a number of times equal to its stoichiometric
	 * coefficient.
         * NOTE: These vectors will be wrong if there are real 
         *       stoichiometric coefficients in the expression.
	 */
        vector<vector_int> m_products;

	/**
	 * m_thermo is a vector of pointers to ThermoPhase
	 * objects. For homogeneous kinetics applications, this vector
	 * will only have one entry. For interfacial reactions, this
	 * vector will consist of multiple entries; some of them will
	 * be surface phases, and the other ones will be bulk phases.
	 * The order that the objects are listed determines the order
	 * in which the species comprising each phase are listed in
	 * the source term vector, originating from the reaction
	 * mechanism.
	 */
        vector<thermo_t*> m_thermo;

	/**
	 * m_start is a vector of integers specifying the beginning position
	 * for the species vector for the n'th phase in the kinetics
	 * class.
	 */
        vector_int  m_start;

        /**
	 * Mapping of the phase id, i.e., the id attribute in the xml
	 * phase element to the position of the phase within the
	 * kinetics object.  Positions start with the value of 1. The
	 * member function, phaseIndex() decrements by one before
	 * returning the index value, so that missing phases return
	 * -1.
	 */
        map<string, int> m_phaseindex;
        int m_index;

	/**
	 * Index in the list of phases of the one surface phase. 
	 */
        int m_surfphase;

        /**
         * Index in the list of phases of the one phase where the reactions
         * occur.
         */
        int m_rxnphase;

        /// number of spatial dimensions of lowest-dimensional phase.
        int m_mindim;

    private:

        vector<grouplist_t> m_dummygroups;
        void err(string m) const;

    };


    typedef Kinetics kinetics_t;

}



#endif
