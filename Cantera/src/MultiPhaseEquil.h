#ifndef CT_MULTIPHASE_EQUIL
#define CT_MULTIPHASE_EQUIL

#include "ct_defs.h"
#include "MultiPhase.h"

namespace Cantera {

    int _equilflag(const char* xy);

    class MultiPhaseEquil {

    public:

        typedef MultiPhase       mix_t;
        typedef size_t           index_t;
        typedef DenseMatrix      matrix_t;

        MultiPhaseEquil(mix_t* mix, bool start=true);

        virtual ~MultiPhaseEquil() {}

        int constituent(index_t m) { 
            if (m < m_nel) return m_order[m]; 
            else return -1;
        }

        void getStoichVector(index_t rxn, vector_fp& nu) {
            index_t k;
            nu.resize(m_nsp, 0.0);
            if (rxn > m_nsp - m_nel) return;
            for (k = 0; k < m_nsp; k++) {
                nu[m_order[k]] = m_N(k, rxn);
            }
        }

        int iterations() { return m_iter; }

        doublereal equilibrate(int XY, doublereal err = 1.0e-9, 
            int maxsteps = 1000, int loglevel=0);

        string reactionString(index_t j);
        doublereal error();
        void printInfo();
        void finish();

    protected:

        void getComponents(const vector_int& order);
        int setInitialMoles();
        void computeN();
        doublereal stepComposition(int loglevel);
        void sort(vector_fp& x);
        void unsort(vector_fp& x);
        void step(doublereal omega, vector_fp& deltaN);
        doublereal computeReactionSteps(vector_fp& dxi);
        void setMoles();

        // moles of the species with sorted index ns
        double moles(int ns) const { return m_moles[m_order[ns]]; }
        double& moles(int ns) { return m_moles[m_order[ns]]; }
        int solutionSpecies(int n) const { return m_dsoln[m_order[n]]; }
        bool isStoichPhase(int n) const { return (m_dsoln[m_order[n]] == 0); }
        
        index_t m_nel_mix, m_nsp_mix, m_np;
        index_t m_nel, m_nsp;
        index_t m_eloc;
        int m_iter;
        mix_t* m_mix;
        doublereal m_press, m_temp;
        vector_int m_order;
        matrix_t m_N, m_A;
        vector_fp m_work, m_work2, m_work3;
        vector_fp m_moles, m_lastmoles, m_dxi;
        vector_fp m_deltaG_RT, m_mu;
        vector<bool> m_majorsp;
        vector_int m_sortindex;
        vector_int m_lastsort;
        vector_int m_dsoln;
        vector_int m_incl_element, m_incl_species;

        // Vector of indices for species that are included in the
        // calculation.  This is used to exclude pure-phase species
        // with invalid thermo data
        vector_int m_species;
        vector_int m_element;
        vector<bool> m_solnrxn;
        bool m_force;
    };

    //-----------------------------------------------------------
    //              convenience functions
    //-----------------------------------------------------------

    /**
     * Set a mixture to a state of chemical equilibrium. The flag 'XY'
     * determines the two properties that will be held fixed in the
     * calculation.
     */
    inline doublereal equilibrate(MultiPhase& s, int XY, 
        doublereal tol = 1.0e-9, int maxsteps = 1000, int loglevel = 0) {
        s.init();
        //MultiPhaseEquil e(&s);
        if (XY == TP || XY == HP) {
            double err = s.equilibrate(XY, tol, maxsteps, maxsteps, loglevel);
            if (loglevel > 0) {
                write_logfile("equilibrate.html");
            }
            return err;
        }
        else {
            throw CanteraError("equilibrate","only fixed T, P supported");
            return -1.0;
        }
    }

    /**
     * Set a mixture to a state of chemical equilibrium. The flag 'XY'
     * determines the two properties that will be held fixed in the
     * calculation.
     */
    inline doublereal equilibrate(MultiPhase& s, const char* XY,
        doublereal tol = 1.0e-9, int maxsteps = 1000, int loglevel = 0) {
        return equilibrate(s,_equilflag(XY), tol, maxsteps, loglevel);
    }
}


#endif
