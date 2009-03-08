/**
 * @file EdgeKinetics.h
 *
 * $Author: dggoodwin $
 * $Revision: 1.6 $
 * $Date: 2006/04/28 17:22:23 $
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_EDGEKINETICS_H
#define CT_EDGEKINETICS_H

#include <fstream>
#include <math.h>
#include <map>
#include <stdlib.h>

#include "mix_defs.h"
#include "Kinetics.h"

#include "utilities.h"
#include "RateCoeffMgr.h"
#include "StoichManager.h"

namespace Cantera {

    // forward references

    class ReactionData;
    class EdgeKineticsData;
    class ThermoPhase;
    class SurfPhase;
    class ImplicitSurfChem;

    /**
     * Holds mechanism-specific data.
     */
    class EdgeKineticsData {
    public:
        EdgeKineticsData() :
            m_ROP_ok(false), 
            m_temp(0.0), m_logtemp(0.0)
            {}
        virtual ~EdgeKineticsData(){}

        doublereal m_logp0, m_logc0;
        array_fp m_ropf, m_ropr, m_ropnet;
        //array_fp m_rfn_low, m_rfn_high;
        bool m_ROP_ok;

        doublereal m_temp, m_logtemp;
        vector_fp m_rfn;
        vector_fp m_rkcn;
    };


    class EdgeKinetics : public Kinetics {

    public:

        /**
	 * Constructor 
	 *
	 */
        EdgeKinetics();

        /// Destructor.
        virtual ~EdgeKinetics();

	/**
	 * Identifies the subclass of the Kinetics manager type.
	 * These are listed in mix_defs.h.
	 */
        virtual int ID() { return cEdgeKinetics; }

	/**
	 * Identifies the subclass of the Kinetics manager type.
	 * These are listed in mix_defs.h.
	 */
        virtual int type() { return cEdgeKinetics; }

	/**
	 * Set the electric potential in the nth phase
	 *
	 * @param n phase Index in this kinetics object. 
	 * @param V Electric potential (volts)
	 */
        void setElectricPotential(int n, doublereal V) {
            thermo(n).setElectricPotential(V);
            m_redo_rates = true;
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
	 *  Units are kmol/m2/s
         */
        virtual void getFwdRatesOfProgress(doublereal* fwdROP) { 
            updateROP(); 
            copy(m_kdata->m_ropf.begin(), m_kdata->m_ropf.end(), fwdROP);
        }

	/**
         * Reverse rates of progress.
         * Return the reverse rates of progress in array revROP, which
         * must be dimensioned at least as large as the total number
         * of reactions.
	 *  Units are kmol/m2/s
         */
        virtual void getRevRatesOfProgress(doublereal* revROP) { 
            updateROP(); 
            copy(m_kdata->m_ropr.begin(), m_kdata->m_ropr.end(), revROP);
        }

	/**
         * Net rates of progress.  Return the net (forward - reverse)
         * rates of progress in array netROP, which must be
         * dimensioned at least as large as the total number of
         * reactions.
	 *  Units are kmol/m2/s
         */
        virtual void getNetRatesOfProgress(doublereal* netROP) { 
            updateROP(); 
            copy(m_kdata->m_ropnet.begin(), m_kdata->m_ropnet.end(), netROP);
        }

        /**
         * Equilibrium constants. Return the equilibrium constants of
         * the reactions in concentration units in array kc, which
         * must be dimensioned at least as large as the total number
         * of reactions.
         */
	virtual void getEquilibriumConstants(doublereal* kc);


	//@}
        /**
         * @name Species Production Rates
         */
        //@{

        /**
         * Species creation rates [kmol/m^2/s]. Return the species 
         * creation rates in array cdot, which must be
         * dimensioned at least as large as the total number of
         * species in all phases of the kinetics
	 * model
         *  
         */ 
        virtual void getCreationRates(doublereal* cdot) {
            updateROP();
            fill(cdot, cdot + m_kk, 0.0);
            m_revProductStoich.incrementSpecies(
                &m_kdata->m_ropf[0], cdot);
            m_irrevProductStoich.incrementSpecies(
                &m_kdata->m_ropf[0], cdot);
            m_reactantStoich.incrementSpecies(
                &m_kdata->m_ropr[0], cdot);
        }

        /**
         * Species destruction rates [kmol/m^2/s]. Return the species 
         * destruction rates in array ddot, which must be
         * dimensioned at least as large as the total number of
         * species in all phases of the kinetics
	 * model
         *  
         */ 
        virtual void getDestructionRates(doublereal* ddot) {
            updateROP();
            fill(ddot, ddot + m_kk, 0.0);
            m_revProductStoich.incrementSpecies(
                &m_kdata->m_ropr[0], ddot);
            m_reactantStoich.incrementSpecies(
                &m_kdata->m_ropf[0], ddot);
        }

	/**
         * Species net production rates [kmol/m^2/s]. Return the species
         * net production rates (creation - destruction) in array
         * wdot, which must be dimensioned at least as large as the
         * total number of species in all phases of the kinetics
	 * model
         */ 
        virtual void getNetProductionRates(doublereal* net) {
            updateROP();
            fill(net, net + m_kk, 0.0);
            m_revProductStoich.incrementSpecies(
                &m_kdata->m_ropnet[0], net);
            m_irrevProductStoich.incrementSpecies(
                &m_kdata->m_ropnet[0], net);
            m_reactantStoich.decrementSpecies(
                &m_kdata->m_ropnet[0], net);
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
            return m_rrxn[k][i];
        }

        /**
         * Stoichiometric coefficient of species k as a product in
         * reaction i.  
         */
        virtual doublereal productStoichCoeff(int k, int i) const {
            return m_prxn[k][i];
        }

        /**
         * Flag specifying the type of reaction. The legal values and
         * their meaning are specific to the particular kinetics
         * manager.
         */
	virtual int reactionType(int i) const {
            return m_index[i].first;
        }

	/**
         * True if reaction i has been declared to be reversible. If
         * isReversible(i) is false, then the reverse rate of progress
         * for reaction i is always zero.
         */
        virtual bool isReversible(int i) {
            if (find(m_revindex.begin(), m_revindex.end(), i) 
                < m_revindex.end()) return true;
            else return false;
        }

	/**
	 * Return a string representing the reaction.
	 */
	virtual string reactionString(int i) const {
            return m_rxneqn[i];
        }

	//@}
        /**
         * @name Reaction Mechanism Construction
         */
        //@{

	/**
	 * Prepare the class for the addition of reactions. This function
	 * must be called after instantiation of the class, but before
	 * any reactions are actually added to the mechanism.
	 * This function calculates m_kk the number of species in all
	 * phases participating in the reaction mechanism. We don't know
	 * m_kk previously, before all phases have been added. 
         */
        virtual void init();

        /**
	 *  Add a single reaction to the mechanism.
	 */
        virtual void addReaction(const ReactionData& r);

        /**
	 * Finish adding reactions and prepare for use. This function
	 * must be called after all reactions are entered into the mechanism
	 * and before the mechanism is used to calculate reaction rates. 
	 */
        virtual void finalize();
        virtual bool ready() const;


        void updateROP();


        const vector<grouplist_t>& reactantGroups(int i)
            { return m_rgroups[i]; }
        const vector<grouplist_t>& productGroups(int i)
            { return m_pgroups[i]; }

        void _update_rates_T();
        void _update_rates_phi();
        void _update_rates_C();
        void checkPartialEquil();

    protected:
	/**
	 * m_kk here is the number of species in all of the phases
	 * that participate in the kinetics mechanism.
	 */
        int                                 m_kk;

        Rate1<SurfaceArrhenius>                    m_rates;        
        //Rate1<Arrhenius>                    m_rates;        
        bool                                m_redo_rates;

	/**
	 * Vector of information about reactions in the
	 * mechanism.
	 * The key is the reaction index (0 < i < m_ii).
	 * The first pair is the reactionType of the reaction.
	 * The second pair is ...
	 */
        mutable map<int, pair<int, int> >   m_index;

        vector<int> m_irrev;

        StoichManagerN                      m_reactantStoich;
        StoichManagerN                      m_revProductStoich;
        StoichManagerN                      m_irrevProductStoich;

        StoichManagerN                      m_globalReactantStoich;

        int m_nirrev;

	/**
	 * Number of reversible reactions in the mechanism
	 */
        int m_nrev;

        map<int, vector<grouplist_t> >      m_rgroups;
        map<int, vector<grouplist_t> >      m_pgroups;

        vector<int>                         m_rxntype;

        mutable vector<map<int, doublereal> >     m_rrxn;
        mutable vector<map<int, doublereal> >     m_prxn;

        vector_int m_revindex;
        vector<string> m_rxneqn;

	/**
	 * Temporary data storage used in calculating the rates of
	 * of reactions.
	 */
        EdgeKineticsData* m_kdata;

	/**
	 * An array of generalized concentrations
         * \f$ C_k \f$ that are defined such that \f$ a_k = C_k /
         * C^0_k, \f$ where \f$ C^0_k \f$ is a standard concentration/
         * These generalized concentrations are used
         * by this kinetics manager class to compute the forward and
         * reverse rates of elementary reactions. The "units" for the
	 * concentrations of each phase depend  upon the implementation
	 * of kinetics within that phase.
	 * The order of the species within the vector is based on
	 * the order of listed ThermoPhase objects in the class, and the
	 * order of the species within each ThermoPhase class.
	 */
        vector_fp m_conc;

        vector_fp m_mu0;
        vector_fp m_phi;
        vector_fp m_pot;
        vector_fp m_rwork;
        vector_fp m_E;
        vector_fp m_beta;
        vector_int m_ctrxn;

    private:

        int reactionNumber(){ return m_ii;}
        void addElementaryReaction(const ReactionData& r);
        void addGlobalReaction(const ReactionData& r);
        void installReagents(const ReactionData& r);

        void installGroups(int irxn, const vector<grouplist_t>& r,
            const vector<grouplist_t>& p);
        void updateKc();

        void registerReaction(int rxnNumber, int type, int loc) {
            m_index[rxnNumber] = pair<int, int>(type, loc);
        }
        void applyButlerVolmerCorrection(doublereal* kf);
        bool m_finalized;
        bool m_has_electrochem_rxns;
    };
}

#endif
