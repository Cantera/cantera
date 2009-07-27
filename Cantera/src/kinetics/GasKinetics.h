/**
 * @file GasKinetics.h
 *
 * @ingroup chemkinetics
 *
 * $Author: hkmoffa $
 * $Revision: 1.8 $
 * $Date: 2008/12/29 21:34:08 $
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_GASKINETICS_H
#define CT_GASKINETICS_H

#include <fstream>
#include <map>

#include "mix_defs.h"
#include "Kinetics.h"

#include "utilities.h"

#include "ReactionStoichMgr.h"
#include "ThirdBodyMgr.h"
#include "FalloffMgr.h"
#include "RateCoeffMgr.h"

#include <cmath>
#include <cstdlib>

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
        array_fp  m_rfn;
        array_fp falloff_work;
        array_fp concm_3b_values;
        array_fp concm_falloff_values;
        array_fp m_rkcn;
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
	 * @name Constructors and General Information
	 */
	//@{
        /// Constructor.
        GasKinetics(thermo_t* thermo = 0);

        /// Destructor.
        virtual ~GasKinetics();

        virtual int ID() const { return cGasKinetics; }
        virtual int type() const { return cGasKinetics; }

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
            std::copy(m_kdata->m_ropf.begin(), m_kdata->m_ropf.end(), fwdROP);
        }

        /**
         * Reverse rates of progress.
         * Return the reverse rates of progress in array revROP, which
         * must be dimensioned at least as large as the total number
         * of reactions.
         */
        virtual void getRevRatesOfProgress(doublereal* revROP) {
            updateROP();
            std::copy(m_kdata->m_ropr.begin(), m_kdata->m_ropr.end(), revROP);
        }

        /**
         * Net rates of progress.  Return the net (forward - reverse)
         * rates of progress in array netROP, which must be
         * dimensioned at least as large as the total number of
         * reactions.
         */
        virtual void getNetRatesOfProgress(doublereal* netROP) {
            updateROP();
            std::copy(m_kdata->m_ropnet.begin(), m_kdata->m_ropnet.end(), netROP);
        }


        /**
         * Equilibrium constants. Return the equilibrium constants of
         * the reactions in concentration units in array kc, which
         * must be dimensioned at least as large as the total number
         * of reactions.
         */
	virtual void getEquilibriumConstants(doublereal* kc);

	/**
	 * Return the array of values for the reaction gibbs free energy
	 * change.
	 * These values depend on the species concentrations.
	 *
	 *  units = J kmol-1
	 */
	virtual void getDeltaGibbs( doublereal* deltaG);

	/**
	 * Return the array of values for the reaction enthalpy change.
	 * These values depend upon the species concentrations.
	 *
	 *  units = J kmol-1
	 */
	virtual void getDeltaEnthalpy( doublereal* deltaH);

	/**
	 * Return the array of values for the reactions change in
	 * entropy.
	 * These values depend upon the concentration
	 * of the solution.
	 *
	 *  units = J kmol-1 Kelvin-1
	 */
	virtual void getDeltaEntropy(doublereal* deltaS);

	/**
	 * Return the array of values for the reaction
	 * standard state Gibbs free energy change.
	 * These values do not depend on the species 
         * concentrations.
	 *
	 *  units = J kmol-1
	 */
	virtual void getDeltaSSGibbs(doublereal* deltaG);

	/**
	 * Return the array of values for the change in the
	 * standard state enthalpies of reaction.
	 * These values do not depend upon the concentration
	 * of the solution.
	 *
	 *  units = J kmol-1
	 */
	virtual void getDeltaSSEnthalpy(doublereal* deltaH);

	/**
	 * Return the array of values for the change in the
	 * standard state entropies for each reaction.
	 * These values do not depend upon the concentration
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

        //! Return the species net production rates
        /*!
         * Species net production rates [kmol/m^3/s]. Return the species
         * net production rates (creation - destruction) in array
         * wdot, which must be dimensioned at least as large as the
         * total number of species.
         *
         *  @param net  Array of species production rates.
         *             units kmol m-3 s-1
         */
        virtual void getNetProductionRates(doublereal* net) {
            updateROP();
            //#ifdef HWMECH
                //get_wdot(&m_kdata->m_ropnet[0], net);
            //#else
            m_rxnstoich->getNetProductionRates(m_kk, 
                &m_kdata->m_ropnet[0], net);
            //#endif
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

        virtual std::string reactionString(int i) const {
            return m_rxneqn[i];
        }

	/**
         * True if reaction i has been declared to be reversible. If
         * isReversible(i) is false, then the reverse rate of progress
         * for reaction i is always zero.
         */
        virtual bool isReversible(int i) {
          if (std::find(m_revindex.begin(), m_revindex.end(), i)
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

        virtual void init();

        ///  Add a reaction to the mechanism.
        virtual void addReaction(const ReactionData& r);

        virtual void finalize();
        virtual bool ready() const;

        virtual void update_T();
        virtual void update_C();

        void updateROP();


        const std::vector<grouplist_t>& reactantGroups(int i)
            { return m_rgroups[i]; }
        const std::vector<grouplist_t>& productGroups(int i)
            { return m_pgroups[i]; }


        void _update_rates_T();
        void _update_rates_C();

      //@}

    protected:

        int                                 m_kk, m_nfall;

        array_int                          m_fallindx;

        Rate1<Arrhenius>                    m_falloff_low_rates;
        Rate1<Arrhenius>                    m_falloff_high_rates;
        Rate1<Arrhenius>                    m_rates;

        mutable std::map<int, std::pair<int, int> >   m_index;

        FalloffMgr                          m_falloffn;

        ThirdBodyMgr<Enhanced3BConc>        m_3b_concm;
        ThirdBodyMgr<Enhanced3BConc>        m_falloff_concm;

        std::vector<int> m_irrev;

        ReactionStoichMgr*                   m_rxnstoich;

        std::vector<int>                         m_fwdOrder;

        int m_nirrev;
        int m_nrev;

        std::map<int, std::vector<grouplist_t> >      m_rgroups;
        std::map<int, std::vector<grouplist_t> >      m_pgroups;

        std::vector<int>                         m_rxntype;

        mutable std::vector<std::map<int, doublereal> >     m_rrxn;
        mutable std::vector<std::map<int, doublereal> >     m_prxn;

	/**
         * Difference between the input global reactants order
         * and the input global products order. Changed to a double
         * to account for the fact that we can have real-valued
         * stoichiometries.
         */
        array_fp  m_dn;
        array_int m_revindex;

        std::vector<std::string> m_rxneqn;

        GasKineticsData* m_kdata;

        array_fp m_conc;
        void processFalloffReactions();
        array_fp m_grt;


    private:

        int reactionNumber(){ return m_ii;}
        std::vector<std::map<int, doublereal> > m_stoich;

        void addElementaryReaction(const ReactionData& r);
        void addThreeBodyReaction(const ReactionData& r);
        void addFalloffReaction(const ReactionData& r);

        void installReagents(const ReactionData& r);

        void installGroups(int irxn, const std::vector<grouplist_t>& r,
            const std::vector<grouplist_t>& p);
        void updateKc();

        void registerReaction(int rxnNumber, int type, int loc) {
            m_index[rxnNumber] = std::pair<int, int>(type, loc);
        }
        bool m_finalized;
    };
}

#endif
