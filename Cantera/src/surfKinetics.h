/**
 * @file SurfKinetics.h
 */


/*
 * $Author$
 * $Revision$
 * $Date$
 *
 * Copyright 2001  California Institute of Technology
 */

#ifndef CT_SURFKINETICS_H
#define CT_SURFKINETICS_H

#include <fstream>
#include <math.h>
#include <map>
#include <stdlib.h>

#include "mix_defs.h"
#include "Kinetics.h"

#include "utilities.h"
#include "RateCoeffMgr.h"
#include "surfacePhase.h"

namespace Cantera {

    // forward references
    class ImplicitSurfChem;
    class ReactionData;

    /**
     * Holds mechanism-specific data.
     * @ingroup kineticsGroup
     */
    class SurfKineticsData {
    public:
        SurfKineticsData() :
            m_ROP_ok(false), 
            m_temp(0.0) 
            {}
        virtual ~SurfKineticsData(){}

        vector_fp m_ropf;
        bool m_ROP_ok;

        doublereal m_temp;
        vector_fp m_rfn;
        doublereal m_s0;
    };


    /**
     * A kinetics manager for elementary surface chemistry. 
     */
    class SurfKinetics : public Kinetics {

    public:

        /// Constructor.
        SurfKinetics() : Kinetics() {}
        SurfKinetics(SurfacePhase* surfphase, thermo_t* th1, 
            thermo_t* th2, string fname="", string id="");

        /// Destructor.
        virtual ~SurfKinetics(){delete m_kdata; delete m_xml;}

        virtual int ID() { return 10; }
        void import(string fname, string id);

        /**
         * The surface phase for which this is a kinetics manager.
         */
        SurfacePhase& sphase() { return *m_surfphase; }

        /**
         * Total number of species on the surface and in both phases.
         */
        int nTotal() { return m_ktot; }

        /**
         * Return a reference to one of the bulk phases.
         */
        Phase* bulkPhase(int n) {
            return m_phase[n];
        }

        /**
         * Get the forward rates of progress for all surface reactions.
         */
        virtual void getFwdRatesOfProgress(doublereal* fwdROP) { 
            updateROP(); 
            copy(m_kdata->m_ropf.begin(), m_kdata->m_ropf.end(), fwdROP);
        }

        /**
         * Get the reverse rates of progress for all surface reactions.
         * All reactions are currently modeled as irreversible, so this
         * returns all zeros.
         */
        virtual void getRevRatesOfProgress(doublereal* revROP) {
            int i;
            for (i = 0; i < m_ii; i++) 
                revROP[i] = 0.0;
        }

        virtual void getNetRatesOfProgress(doublereal* netROP) { 
            getFwdRatesOfProgress(netROP);
        }

        virtual void getNetProductionRates(doublereal* net);
        virtual void getCreationRates(doublereal* cdot);
        virtual void getDestructionRates(doublereal* ddot);
        virtual void getChemRates(doublereal* rtau);

        virtual void init();

        virtual void integrate(doublereal dt);

        ///  Add a reaction to the mechanism.
        void addReaction(const vector_int& r, const vector_int& rstoich, 
            const vector_int& order, const vector_int& p, 
            const vector_int& pstoich, 
            const vector_fp& rateParams);

        void saveReactionData(const vector_int& r, const vector_int& rstoich, 
            const vector_int& order, const vector_int& p, 
            const vector_int& pstoich, 
            const vector_fp& rateParams);

        virtual void finalize();
        virtual bool ready() const;

        void updateROP();

        virtual int reactionType(int i) const {return SURFACE_RXN;}
        virtual bool isReversible(int i) {return false;}
        void save(string fname, string idtag, string comment);

    protected:

        SurfacePhase* m_surfphase;
        SurfKineticsData* m_kdata;

        // objects for bulk phase 2. Those for bulk phase 1
        // are declared in Kinetics
        phase_t* m_phase2;
        thermo_t* m_thermo2;

        int m_kk;   // number of surface species
        int m_kk1;  // number of bulk phase 1 species
        int m_kk2;  // number of bulk phase 2 species
        int m_ktot; // total number of species

        Rate1<Arrhenius>                    m_rates;
        
        vector_int m_irrev;
        int m_nirrev;

        //mutable vector<map<int, doublereal> >     m_rrxn;
        //mutable vector<map<int, doublereal> >     m_prxn;

        //map<int, map<int, doublereal> >  m_rstoich;
        //map<int, map<int, doublereal> >  m_pstoich;

        vector<vector_int>               m_order;
        vector<vector_int>               m_rst;
        vector<vector_int>               m_pst;
        vector_int                       m_nr, m_np;

        vector_fp m_conc;

        ImplicitSurfChem* m_integrator;
        map<string, int>                 m_bsp1, m_bsp2;

    private:

        void _update_rates_T();
        void _update_rates_C();
        XML_Node* m_xml;
        bool m_twobulk;

        bool m_finalized;
    };
}

#endif
