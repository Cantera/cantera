/**
 * @file Kinetics.h
 *
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology

/**
 * @defgroup kineticsGroup Kinetics
 */

#ifndef CT_KINETICS_H
#define CT_KINETICS_H

#include "ctexceptions.h"
//#include "Phase.h"
#include "ThermoPhase.h"

namespace Cantera {

    class ReactionData;

    /**
     * Public interface for kinetics managers. This class serves as a
     * base class to derive 'kinetics managers', which are classes
     * that manage homogeneous chemistry within one phase.
     */
    class Kinetics {

    public:

        // typedefs
        typedef ThermoPhase thermo_t;

        /// Constructor.
        Kinetics() : m_ii(0), m_thermo(0), m_index(-1) {}
        Kinetics(thermo_t* thermo) 
            : m_ii(0), m_index(-1) {
            if (thermo) {
                m_start.push_back(0);
                m_thermo.push_back(thermo);
            }
        }

        /// Destructor. Does nothing.
        virtual ~Kinetics() {} // delete m_xml; }

        int index(){ return m_index; }
        void setIndex(int index) { m_index = index; }

        //XML_Node& xml() { return *m_xml; }

        /// Identifies subclass.
        virtual int type() { return 0; }

        int start(int n) { return m_start[n]; }

        /// Number of reactions
        int nReactions() const {return m_ii;}

        /// Number of species
        int nTotalSpecies() const {
            int n=0, np;
            np = nPhases();
            for (int p = 0; p < np; p++) n += thermo(p).nSpecies();
            return n;
        }

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


        /**
         * Returns a read-only reference to the vector of reactant
         * index numbers for reaction i.
         */ 
        virtual const vector_int& reactants(int i) const { 
            return m_reactants[i]; 
        }
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
         * @name Reaction Rates Of Progress
         */
        //@{

        /**
         * Forward rates of progress.
         * Return the forward rates of progress in array fwdROP, which
         * must be dimensioned at least as large as the total number
         * of reactions.
         */
        virtual void getFwdRatesOfProgress(doublereal* fwdROP) {
            err("getFwdRatesOfProgress");
        }
 
        /**
         * Reverse rates of progress.
         * Return the reverse rates of progress in array revROP, which
         * must be dimensioned at least as large as the total number
         * of reactions.
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
         * True if reaction i has been declared to be reversible. If
         * isReversible(i) is false, then the reverse rate of progress
         * for reaction i is always zero.
         */
        virtual bool isReversible(int i){return false;}

        /**
         * Species creation rates [kmol/m^3]. Return the species 
         * creation rates in array cdot, which must be
         * dimensioned at least as large as the total number of
         * species.
         *  
         */ 
        virtual void getCreationRates(doublereal* cdot) {
            err("getCreationRates");
        }

        /**
         * Species destruction rates [kmol/m^3]. Return the species 
         * destruction rates in array ddot, which must be
         * dimensioned at least as large as the total number of
         * species.
         *  
         */ 
        virtual void getDestructionRates(doublereal* ddot) {
            err("getDestructionRates");
        }

        /**
         * Species net production rates [kmol/m^3]. Return the species
         * net production rates (creation - destruction) in array
         * wdot, which must be dimensioned at least as large as the
         * total number of species.
         */ 
        virtual void getNetProductionRates(doublereal* wdot) {
            err("getNetProductionRates");
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
        //@}


        /**
         * @name Reaction Mechanism Construction
         */
        //@{


        /**
         * Return the number of phases defined within the kinetics
	 * object.
         */  
        int nPhases() const { return m_thermo.size(); }
	/**
	 * Return the phase index of a phase in the list of phases
	 * defined within the object.
	 *
	 *  Input
	 * ----------
	 *  ph = string name of the phase
	 *
	 * If a -1 is returned, then the phase is not defined in
	 * the Kinetics object.
	 * (HKM -> unfound object will create another entry in the
	 *  map, suggest rewriting this function)
	 */
        int phaseIndex(string ph) { return m_phaseindex[ph] - 1; }

        /**
         * Add a phase.
         */
        void addPhase(thermo_t& thermo) { 
            if (m_thermo.size() > 0) {
                m_start.push_back(m_start.back() 
                    + m_thermo.back()->nSpecies());
            }
            else {
                m_start.push_back(0);
            }
            m_thermo.push_back(&thermo);
            m_phaseindex[m_thermo.back()->id()] = nPhases();
        }

	/**
	 * This method returns a reference to the nth ThermoPhase
	 * defined in this kinetics mechanism.
	 * It is typically used so that member functions of the the
	 * ThermoPhase may be called.
	 */
        thermo_t& thermo(int n=0) { return *m_thermo[n]; }
        const thermo_t& thermo(int n=0) const { return *m_thermo[n]; }

	/**
	 * This method returns a reference to the nth thermophase
	 * defined in this kinetics mechanism.
	 * It is typically used so that member functions of the
	 * ThermoPhase may be called.
	 */
        thermo_t& phase(int n=0) { return *m_thermo[n]; }
        const thermo_t& phase(int n=0) const { return *m_thermo[n]; }

	/**
	 * This method returns the index of a species in the source
	 * term vector for this kinetics object.
	 *
	 * @param k species index 
	 * @param n phase index for the species
	 */
        int kineticsSpeciesIndex(int k, int n) {
            return m_start[n] + k;
        }
	/**
	 * This routine will look up a species number based on
	 * the input string nm. The lookup of species will
	 * occur for all phases listed in the kinetics obect,
	 * unless the string ph refers to a specific phase of
	 * the object. 
	 *
	 *  return
	 *   If a match is found, the position in the species list
	 *   is returned. 
	 *   If no match is found, the value -2 is returned.
	 */
        int kineticsSpeciesIndex(string nm, string ph = "<any>") {
	  int np = m_thermo.size();
	  int k;
	  string id;
	  for (int n = 0; n < np; n++) {
	    id = thermo(n).id();
	    if (ph == id) {
	      k = thermo(n).speciesIndex(nm);
	      if (k < 0) return -1;
	      return k + m_start[n];
	    }
	    else if (ph == "<any>") {
	      /*
	       * Call the speciesIndex() member function of the
	       * ThermoPhase object to find a match.
	       */
	      k = thermo(n).speciesIndex(nm);
	      if (k >= 0) return k + m_start[n];
	    }                    
	  }
	  return -2;
        }

        thermo_t& speciesPhase(string nm) {
            int np = m_thermo.size();
            int k;
            string id;
            for (int n = 0; n < np; n++) {
                k = thermo(n).speciesIndex(nm);
                if (k >= 0) return thermo(n);
            }
            throw CanteraError("speciesPhase", "unknown species "+nm);
        }


        /**
         * Prepare the class for the addition of reactions. This function
	 * must be called after instantiation of the class, but before
	 * any reactions are actually added to the mechanism.
         */
        virtual void init() {err("init");}

        /**
	 * Finish adding reactions and prepare for use. This function
	 * must be called after all reactions are entered into the mechanism
	 * and before the mechanism is used to calculate reaction rates. 
	 */
        virtual void finalize() {err("finalize");}

        virtual void addReaction(const ReactionData& r) {err("addReaction");}

        virtual string reactionString(int i) const {
            err("reactionString"); return "<null>";
        }


        virtual const vector<grouplist_t>& reactantGroups(int i) { 
	    err("reactantGroups"); 
	    return m_dummygroups;
	}

        virtual const vector<grouplist_t>& productGroups(int i) {
	    err("productGroups"); 
	    return m_dummygroups;
	}

	//@}
        /**
         * @name Altering Reaction Rates
         *
         * These methods alter reaction rates. They are designed
         * primarily for carrying out sensitivity analysis.
         */
        //@{

        /// The current value of the multiplier for reaction i.
        doublereal multiplier(int i) const {return m_perturb[i];}

        /// Set the multiplier for reaction i to f.
        void setMultiplier(int i, doublereal f) {m_perturb[i] = f;}
        
        //@}
	/**
	 * Increment the number of reactions in the mechanism by one.
	 */
        void incrementRxnCount() { m_ii++; m_perturb.push_back(1.0); }

	/**
	 *
	 */
        virtual bool ready() const {
	    err("ready()");
	    return false;
	}

    protected:

	/**
	 * m_ii is the number of reactions in the mechanism
	 */
        int m_ii;

	/**
	 * m_perturb is a vector of perturbation factors for each
	 * reaction's rate of progress vector. It is initialized to one.
	 */
        vector_fp m_perturb;

	/**
	 * This is a vector of vectors containing the reactants for each
	 * reaction. The outer vector is over the number of reactions, m_ii.
	 * The inner vector is a list of species indecises. If the stoichiometric
	 * coefficient for a reactant is greater than one, then the
	 * reactant is listed contiguously in the vector a number of 
	 * times equal to its stoichiometric coefficient.
	 */
        vector<vector_int> m_reactants;

	/**
	 * This is a vector of vectors containing the products for each
	 * reaction. The outer vector is over the number of reactions, m_ii.
	 * The inner vector is a list of species indecises. If the stoichiometric
	 * coefficient for a product is greater than one, then the
	 * reactant is listed contiguously in the vector a number of 
	 * times equal to its stoichiometric coefficient.
	 */

        vector<vector_int> m_products;

	/**
	 * m_thermo is a vector of pointers to ThermoPhase objects. For 
	 * homogeneous kinetics applications, this vector will only
	 * consist of one entry. For interfacial reactions, this vector
	 * will consist of multiple entries; some of them will be surface
	 * phases, and the other ones will be bulk phases. 
	 * The order that the objects are listed determines the order in
	 * which the species comprising each phase are listed in the
	 * source term vector, originating from the reaction mechanism.
	 */
        vector<thermo_t*> m_thermo;

	/**
	 * m_start is a vector of integers specifying the beginning position
	 * for the species vector for the n'th phase in the kinetics
	 * class.
	 */
        vector_int  m_start;

        //        XML_Node* m_xml;
        map<string, int> m_phaseindex;
        int m_index;

    private:

        vector<grouplist_t> m_dummygroups;
        void err(string m) const {
            throw CanteraError("Kinetics::" + m, 
			       "The default Base class method was called, when "
		               "the inherited class's method should have been called");
        }
    };

    typedef Kinetics kinetics_t;
}



#endif
