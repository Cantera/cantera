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
         * Get the nth Phase object.
         */
        //phase_t& phase(int n=0) { return *m_phase[n]; }
        //const phase_t& phase(int n=0) const { return *m_phase[n]; }
        int nPhases() const { return m_thermo.size(); }
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
        thermo_t& thermo(int n=0) { return *m_thermo[n]; }
        const thermo_t& thermo(int n=0) const { return *m_thermo[n]; }
        thermo_t& phase(int n=0) { return *m_thermo[n]; }
        const thermo_t& phase(int n=0) const { return *m_thermo[n]; }

        int kineticsSpeciesIndex(int k, int n) {
            return m_start[n] + k;
        }

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
                    k = thermo(n).speciesIndex(nm);
                    if (k >= 0) return k + m_start[n];
                }                    
            }
            return -2;
        }

        /**
         * Prepare to add reactions.
         */
        virtual void init() {err("init");}

        /// Finished adding reactions. Prepare for use.
        virtual void finalize() {err("finalize");}

        virtual void addReaction(const ReactionData& r) {err("addReaction");}

        virtual string reactionString(int i) const {
            err("reactionString"); return "<null>";
        }


        virtual const vector<grouplist_t>& reactantGroups(int i)  
            { err("reactantGroups"); return m_dummygroups; }

        virtual const vector<grouplist_t>& productGroups(int i) 
            { err("productGroups"); return m_dummygroups; }

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

        void incrementRxnCount() { m_ii++; m_perturb.push_back(1.0); }

        virtual bool ready() const {return false;}

    protected:

        int m_ii;
        vector_fp m_perturb;
        vector<vector_int> m_reactants;
        vector<vector_int> m_products;
        vector<thermo_t*> m_thermo;
        vector_int  m_start;
        //        XML_Node* m_xml;
        map<string, int> m_phaseindex;
        int m_index;

    private:

        vector<grouplist_t> m_dummygroups;
        void err(string m) const {
            throw CanteraError("Kinetics::"+m,"Base class method called.");
        }
    };

    typedef Kinetics kinetics_t;
}



#endif
