/**
 * @file InterfaceKinetics.h
 *
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_IFACEKINETICS_H
#define CT_IFACEKINETICS_H

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
    class InterfaceKineticsData;
    class ThermoPhase;

    /**
     * Holds mechanism-specific data.
     */
    class InterfaceKineticsData {
    public:
        InterfaceKineticsData() :
            m_ROP_ok(false), 
            m_temp(0.0)
            {}
        virtual ~InterfaceKineticsData(){}

        doublereal m_logp0, m_logc0;
        array_fp m_ropf, m_ropr, m_ropnet;
        array_fp m_rfn_low, m_rfn_high;
        bool m_ROP_ok;

        doublereal m_temp;
        vector_fp m_rfn;
        vector_fp m_rkcn;
    };


    class InterfaceKinetics : public Kinetics {

    public:

        /// Constructor.
        InterfaceKinetics(thermo_t* thermo = 0);

        /// Destructor.
        virtual ~InterfaceKinetics(){delete m_kdata;}

        virtual int ID() { return cInterfaceKinetics; }

        void setElectricPotential(int n, doublereal V) {
            m_phi[n] = V;
            m_redo_rates = true;
        }

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
            fill(net, net + m_kk, 0.0);
            m_revProductStoich.incrementSpecies(
                m_kdata->m_ropnet.begin(), net);
            m_irrevProductStoich.incrementSpecies(
                m_kdata->m_ropnet.begin(), net);
            m_reactantStoich.decrementSpecies(
                m_kdata->m_ropnet.begin(), net);
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

        virtual void init();

        ///  Add a reaction to the mechanism. 
        virtual void addReaction(const ReactionData& r);

        virtual void finalize();
        virtual bool ready() const;

        //virtual void update_T();
        //virtual void update_C();

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
        void correctElectronTransferRates(doublereal* kf);
        void _update_rates_T();
        void _update_rates_C();

    protected:

        int                                 m_kk;

        Rate1<Arrhenius>                    m_rates;        
        bool                                m_redo_rates;

        mutable map<int, pair<int, int> >   m_index;

        vector<int> m_irrev;

        StoichManagerN                      m_reactantStoich;
        StoichManagerN                      m_revProductStoich;
        StoichManagerN                      m_irrevProductStoich;

        StoichManagerN                      m_globalReactantStoich;

        int m_nirrev;
        int m_nrev;

        map<int, vector<grouplist_t> >      m_rgroups;
        map<int, vector<grouplist_t> >      m_pgroups;

        vector<int>                         m_rxntype;

        mutable vector<map<int, doublereal> >     m_rrxn;
        mutable vector<map<int, doublereal> >     m_prxn;

        vector_int m_revindex;
        vector<string> m_rxneqn;

        InterfaceKineticsData* m_kdata;

        vector_fp m_conc;
        vector_fp m_mu0;
        vector_fp m_phi;
        vector_fp m_pot;
        vector_fp m_rwork;

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
        bool m_finalized;
    };
}

#endif
