/**
 * @file GasKinetics.h
 *
 * $Author$
 * $Revision$
 * $Date$
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
#include "StoichManager.h"
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
            m_ROP_ok(false), 
            m_temp(0.0)
            {}
        virtual ~GasKineticsData(){}

        doublereal m_logp0, m_logc0;
        array_fp m_ropf, m_ropr, m_ropnet;
        array_fp m_rfn_low, m_rfn_high;
        bool m_ROP_ok;

        doublereal m_temp;
        vector_fp m_rfn;
        vector_fp falloff_work;
        vector_fp concm_3b_values;
        vector_fp concm_falloff_values;
        vector_fp m_rkcn;
    };


    /**
     * Kinetics manager for elementary gas-phase chemistry. This
     * kinetics manager implements standard mass-action reaction rate
     * expressions for low-density gases. It assumes that all
     * stoichiometric coefficients are integers.
     */

    class GasKinetics : public Kinetics {

    public:

        /// Constructor.
        GasKinetics(thermo_t* thermo = 0);

        /// Destructor.
        virtual ~GasKinetics(){delete m_kdata;}

        virtual int ID() { return cGasKinetics; }

        virtual doublereal reactantStoichCoeff(int k, int i) const {
            return m_rrxn[k][i];
        }

        virtual doublereal productStoichCoeff(int k, int i) const {
            return m_prxn[k][i];
        }

        virtual void getFwdRatesOfProgress(doublereal* fwdROP) { 
            updateROP(); 
            copy(m_kdata->m_ropf.begin(), m_kdata->m_ropf.end(), fwdROP);
        }

        virtual void getRevRatesOfProgress(doublereal* revROP) { 
            updateROP(); 
            copy(m_kdata->m_ropr.begin(), m_kdata->m_ropr.end(), revROP);
        }

        virtual void getNetRatesOfProgress(doublereal* netROP) { 
            updateROP(); 
            copy(m_kdata->m_ropnet.begin(), m_kdata->m_ropnet.end(), netROP);
        }

        virtual void getNetProductionRates(doublereal* net) {
            updateROP();
#ifdef HWMECH
            get_wdot(m_kdata->m_ropnet.begin(), net);
#else
            fill(net, net + m_kk, 0.0);
            m_revProductStoich.incrementSpecies(
                m_kdata->m_ropnet.begin(), net);
            m_irrevProductStoich.incrementSpecies(
                m_kdata->m_ropnet.begin(), net);
            m_reactantStoich.decrementSpecies(
                m_kdata->m_ropnet.begin(), net);
#endif
        }

        virtual void getCreationRates(doublereal* cdot) {
            updateROP();
            fill(cdot, cdot + m_kk, 0.0);
            m_revProductStoich.incrementSpecies(
                m_kdata->m_ropf.begin(), cdot);
            m_irrevProductStoich.incrementSpecies(
                m_kdata->m_ropf.begin(), cdot);
            m_reactantStoich.incrementSpecies(
                m_kdata->m_ropr.begin(), cdot);
        }

        virtual void getDestructionRates(doublereal* ddot) {
            updateROP();
            fill(ddot, ddot + m_kk, 0.0);
            m_revProductStoich.incrementSpecies(
                m_kdata->m_ropr.begin(), ddot);
            m_reactantStoich.incrementSpecies(
                m_kdata->m_ropf.begin(), ddot);
        }

        virtual void getEquilibriumConstants(doublereal* kc);

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

        virtual int reactionType(int i) const {
            return m_index[i].first;
        }

        virtual string reactionString(int i) const {
            return m_rxneqn[i];
        }

        const vector<grouplist_t>& reactantGroups(int i)
            { return m_rgroups[i]; }
        const vector<grouplist_t>& productGroups(int i)
            { return m_pgroups[i]; }

        virtual bool isReversible(int i) {
            if (find(m_revindex.begin(), m_revindex.end(), i) 
                < m_revindex.end()) return true;
            else return false;
        }

        void _update_rates_T();
        void _update_rates_C();

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

        StoichManagerN                      m_reactantStoich;
        StoichManagerN                      m_revProductStoich;
        StoichManagerN                      m_irrevProductStoich;

        vector<int>                         m_fwdOrder;

        int m_nirrev;
        int m_nrev;

        map<int, vector<grouplist_t> >      m_rgroups;
        map<int, vector<grouplist_t> >      m_pgroups;

        vector<int>                         m_rxntype;

        mutable vector<map<int, doublereal> >     m_rrxn;
        mutable vector<map<int, doublereal> >     m_prxn;

        vector_int m_dn;
        vector_int m_revindex;

        map<int, map<int, doublereal> >  m_rstoich;
        map<int, map<int, doublereal> >  m_pstoich;

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
