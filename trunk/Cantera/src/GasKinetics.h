/**
 * @file GasKinetics.h
 *
 * $Author: hkmoffa $
 * $Revision: 1.9 $
 * $Date: 2006/04/30 18:01:42 $
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_GASKINETICS_H
#define CT_GASKINETICS_H

#include <fstream>
#include <math.h>
#include <map>
#include <stdlib.h>

#include "mix_defs.h"
#include "Kinetics.h"

#include "utilities.h"

#include "ReactionStoichMgr.h"
#include "ThirdBodyMgr.h"
#include "FalloffMgr.h"
#include "RateCoeffMgr.h"

void get_wdot(const doublereal* rop, doublereal* wdot);

namespace Cantera {

    // forward references

    class Enhanced3BConc;
    class ReactionData;
    class GasKineticsData;
    class Thermo;

    /**
     * Holds mechanism-specific data.
     */
    class GasKineticsData {
    public:
        GasKineticsData() :
	    m_logp_ref(0.0),
	    m_logc_ref(0.0),
	    m_logStandConc(0.0),
            m_ROP_ok(false), 
            m_temp(0.0)
            {}
        virtual ~GasKineticsData(){}

        doublereal m_logp_ref, m_logc_ref, m_logStandConc;
        array_fp m_ropf, m_ropr, m_ropnet;
        array_fp m_rfn_low, m_rfn_high;
        bool m_ROP_ok;

        doublereal m_temp;
        vector_fp  m_rfn;
        vector_fp falloff_work;
        vector_fp concm_3b_values;
        vector_fp concm_falloff_values;
        vector_fp m_rkcn;
    };


    /**
     * Kinetics manager for elementary gas-phase chemistry. This
     * kinetics manager implements standard mass-action reaction rate
     * expressions for low-density gases.
     * @ingroup kinetics
     */

    class GasKinetics : public Kinetics {

    public:
	/**
	 * @name Constructors and General Information about Mechanism
	 */
	//@{
        /// Constructor.
        GasKinetics(thermo_t* thermo = 0);

        /// Destructor.
        virtual ~GasKinetics();

        virtual int ID() { return cGasKinetics; }

        virtual doublereal reactantStoichCoeff(int k, int i) const {
            return m_rrxn[k][i];
        }

        virtual doublereal productStoichCoeff(int k, int i) const {
            return m_prxn[k][i];
        }

	//@}
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
            updateROP(); 
            copy(m_kdata->m_ropf.begin(), m_kdata->m_ropf.end(), fwdROP);
        }
        /**
         * Reverse rates of progress.
         * Return the reverse rates of progress in array revROP, which
         * must be dimensioned at least as large as the total number
         * of reactions.
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

	/**
	 * Return the vector of values for the reaction gibbs free energy
	 * change.
	 * These values depend upon the concentration
	 * of the solution.
	 *
	 *  units = J kmol-1
	 */
	virtual void getDeltaGibbs( doublereal* deltaG);

	/**
	 * Return the vector of values for the reactions change in
	 * enthalpy.
	 * These values depend upon the concentration
	 * of the solution.
	 *
	 *  units = J kmol-1
	 */
	virtual void getDeltaEnthalpy( doublereal* deltaH);

	/**
	 * Return the vector of values for the reactions change in
	 * entropy.
	 * These values depend upon the concentration
	 * of the solution.
	 *
	 *  units = J kmol-1 Kelvin-1
	 */
	virtual void getDeltaEntropy(doublereal* deltaS);

	/**
	 * Return the vector of values for the reaction 
	 * standard state gibbs free energy change.
	 * These values don't depend upon the concentration
	 * of the solution.
	 *
	 *  units = J kmol-1
	 */
	virtual void getDeltaSSGibbs(doublereal* deltaG);

	/**
	 * Return the vector of values for the change in the
	 * standard state enthalpies of reaction.
	 * These values don't depend upon the concentration
	 * of the solution.
	 *
	 *  units = J kmol-1
	 */
	virtual void getDeltaSSEnthalpy(doublereal* deltaH);

	/**
	 * Return the vector of values for the change in the
	 * standard state entropies for each reaction.
	 * These values don't depend upon the concentration
	 * of the solution.
	 *
	 *  units = J kmol-1 Kelvin-1
	 */
	virtual void getDeltaSSEntropy(doublereal* deltaS);

	//@}
        /**
         * @name Species Production Rates
         */
        //@{

        /**
         * Species net production rates [kmol/m^3]. Return the species
         * net production rates (creation - destruction) in array
         * wdot, which must be dimensioned at least as large as the
         * total number of species.
         */ 
        virtual void getNetProductionRates(doublereal* net) {
            updateROP();
#ifdef HWMECH
            get_wdot(&m_kdata->m_ropnet[0], net);
#else
            m_rxnstoich->getNetProductionRates(m_kk, &m_kdata->m_ropnet[0], net); 
            //fill(net, net + m_kk, 0.0);
            //m_revProductStoich.incrementSpecies(
            //    m_kdata->m_ropnet.begin(), net);
            //m_irrevProductStoich.incrementSpecies(
            //    m_kdata->m_ropnet.begin(), net);
            //m_reactantStoich.decrementSpecies(
            //    m_kdata->m_ropnet.begin(), net);
#endif
        }

	/**
         * Species creation rates [kmol/m^3]. Return the species 
         * creation rates in array cdot, which must be
         * dimensioned at least as large as the total number of
         * species.
         *  
         */ 
        virtual void getCreationRates(doublereal* cdot) {
            updateROP();
            m_rxnstoich->getCreationRates(m_kk, &m_kdata->m_ropf[0], 
                &m_kdata->m_ropr[0], cdot); 
            //fill(cdot, cdot + m_kk, 0.0);
            //m_revProductStoich.incrementSpecies(
            //    m_kdata->m_ropf.begin(), cdot);
            //m_irrevProductStoich.incrementSpecies(
            //    m_kdata->m_ropf.begin(), cdot);
            //m_reactantStoich.incrementSpecies(
            //    m_kdata->m_ropr.begin(), cdot);
        }

	/**
         * Species destruction rates [kmol/m^3]. Return the species 
         * destruction rates in array ddot, which must be
         * dimensioned at least as large as the total number of
         * species.
         *  
         */ 
        virtual void getDestructionRates(doublereal* ddot) {
            updateROP();
            m_rxnstoich->getDestructionRates(m_kk, &m_kdata->m_ropf[0], 
                &m_kdata->m_ropr[0], ddot); 
            //            fill(ddot, ddot + m_kk, 0.0);
            //m_revProductStoich.incrementSpecies(
            //    m_kdata->m_ropr.begin(), ddot);
            //m_reactantStoich.incrementSpecies(
            //    m_kdata->m_ropf.begin(), ddot);
        }

        //@}
        /**
         * @name Reaction Mechanism Informational Query Routines
         */
        //@{

        /**
         * Flag specifying the type of reaction. The legal values and
         * their meaning are specific to the particular kinetics
         * manager.
         */
	virtual int reactionType(int i) const {
            return m_index[i].first;
        }

        virtual string reactionString(int i) const {
            return m_rxneqn[i];
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
	 * Return the forward rate constants
	 *
	 * length is the number of reactions. units depends
	 * on many issues.
	 */
	virtual void getFwdRateConstants(doublereal *kfwd);

	/**
	 * Return the reverse rate constants.
	 *
	 * length is the number of reactions. units depends
	 * on many issues. Note, this routine will return rate constants
	 * for irreversible reactions if the default for
	 * doIrreversible is overridden.
	 */
	virtual void getRevRateConstants(doublereal *krev, 
					 bool doIrreversible = false);

        //@}
        /**
         * @name Reaction Mechanism Setup Routines
         */
        //@{


	/**
         * Set delta T threshold for updating temperature-dependent
         * rates.
         */
        void setRateUpdateThreshold(doublereal dt) {
            m_dt_threshold = dt;
        }

        virtual void init();

        ///  Add a reaction to the mechanism. 
        void addReaction(const ReactionData& r);

        virtual void finalize();
        virtual bool ready() const;

        virtual void update_T();
        virtual void update_C();

        void updateROP();


        const vector<grouplist_t>& reactantGroups(int i)
            { return m_rgroups[i]; }
        const vector<grouplist_t>& productGroups(int i)
            { return m_pgroups[i]; }


        void _update_rates_T();
        void _update_rates_C();

      //@}

    protected:

        int                                 m_kk, m_nfall;

        vector_int                          m_fallindx;
        doublereal                          m_dt_threshold;

        Rate1<Arrhenius>                    m_falloff_low_rates;
        Rate1<Arrhenius>                    m_falloff_high_rates;        
        Rate1<Arrhenius>                    m_rates;        
        
        mutable map<int, pair<int, int> >   m_index;
        
        FalloffMgr                          m_falloffn;
        
        ThirdBodyMgr<Enhanced3BConc>        m_3b_concm;
        ThirdBodyMgr<Enhanced3BConc>        m_falloff_concm;
        
        vector<int> m_irrev;

        //StoichManagerN                      m_reactantStoich;
        //StoichManagerN                      m_revProductStoich;
        //StoichManagerN                      m_irrevProductStoich;

        ReactionStoichMgr*                   m_rxnstoich;

        vector<int>                         m_fwdOrder;

        int m_nirrev;
        int m_nrev;

        map<int, vector<grouplist_t> >      m_rgroups;
        map<int, vector<grouplist_t> >      m_pgroups;

        vector<int>                         m_rxntype;

        mutable vector<map<int, doublereal> >     m_rrxn;
        mutable vector<map<int, doublereal> >     m_prxn;

	/**
         * Difference between the input global reactants order
         * and the input global products order. Changed to a double
         * to account for the fact that we can have real-valued
         * stoichiometries.
         */
        vector_fp  m_dn;
        vector_int m_revindex;

        vector<string> m_rxneqn;

        GasKineticsData* m_kdata;

        vector_fp m_conc;
        void processFalloffReactions();
        vector_fp m_grt;


    private:

        int reactionNumber(){ return m_ii;}
        vector<map<int, doublereal> > m_stoich;

        void addElementaryReaction(const ReactionData& r);
        void addThreeBodyReaction(const ReactionData& r);
        void addFalloffReaction(const ReactionData& r);
        
        void installReagents(const ReactionData& r); 
        //const vector_int& r,
        //    const vector_int& p, bool reversible);

        void installGroups(int irxn, const vector<grouplist_t>& r,
            const vector<grouplist_t>& p);
        void updateKc();

        void registerReaction(int rxnNumber, int type, int loc) {
            m_index[rxnNumber] = pair<int, int>(type, loc);
        }
        bool m_finalized;
    };
}

#endif
