/**
 * @file OneDim.h
 */

#ifndef CT_ONEDIM_H
#define CT_ONEDIM_H

#include "Resid1D.h"

namespace Cantera {

    class MultiJac;
    class MultiNewton;

    /**
     * Container class for multiple-domain 1D problems. Each domain is
     * represented by an instance of Resid1D.
     */
    class OneDim {
    public: 

        // Default constructor.
        OneDim();

        // Constructor.
        OneDim(vector<Resid1D*> domains);

        /// Destructor.
        virtual ~OneDim();

        /// Add a domain.
        void addDomain(Resid1D* d);

        /// Return a reference to the Jacobian.
        MultiJac& jacobian();

        /// Return a reference to the Newton iterator.
        MultiNewton& newton();

        /// Solve F(x) = 0, where F(x) is the multi-domain residual.
        int solve(doublereal* x, doublereal* xnew, int loglevel);

        /// Number of domains.
        int nDomains() const { return m_nd; }

        /// Return a reference to domain i.
        Resid1D& domain(int i) const { return *m_dom[i]; }

        /// The index of the start of domain i in the solution vector.
        int start(int i) const { return m_dom[i]->loc(); }

        /// Total solution vector length;
        int size() const { return m_size; }

        /// Pointer to left-most domain (first added).
        Resid1D* left() { return m_dom[0]; }

        /// Pointer to right-most domain (last added).
        Resid1D* right() { return m_dom.back(); }

        /// Number of solution components at global point jg.
        int nVars(int jg) { return m_nvars[jg]; }

        /** Location in the solution vector of the first component of
            global point jg. */
        int loc(int jg) { return m_loc[jg]; }

        /// Jacobian bandwidth.
        int bandwidth() const { return m_bw; }

        /// Initialize.
        void init();

        /// Total number of points.
        int points() { return m_pts; }

        /// Staedy-state max norm of the residual.
        doublereal ssnorm(doublereal* x, doublereal* r);

        /// Reciprocal of the time step. 
        doublereal rdt() const { return m_rdt; }
        
        /// Prepare for time stepping.
        void initTimeInteg(doublereal dt, doublereal* x);

        /// True if transient mode.
        bool transient() const { return (m_rdt != 0.0);}

        /// True if steady mode.
        bool steady() const { return (m_rdt == 0.0); }

        /// Set steady mode
        void setSteadyMode();

        /// Evaluate the multi-domain residual function 
        void eval(int j, double* x, double* r, doublereal rdt=-1.0, 
            int count = 1);

        /// Pointer to the domain global point i belongs to.
        Resid1D* pointDomain(int i);

        void resize();
        doublereal solveTime() { return m_solve_time; }

        //void setTransientMask();
        vector_int& transientMask() { return m_mask; }

        double timeStep(int nsteps, double dt, double* x, 
            double* r, int loglevel);

        void writeStats();
        void saveStats();

        void save(string fname, string id, string desc, doublereal* sol);

    protected:

        MultiJac* m_jac;
        MultiNewton* m_newt;
        doublereal m_rdt;
        bool m_jac_ok;
        int m_nd;
        int m_bw;
        int m_size;
        vector_int m_states;
        vector_int m_start;
        vector_int m_comp, m_points;
        vector<Resid1D*> m_dom, m_connect, m_bulk;
        vector_int m_flow, m_bdry;
        int m_nflow, m_nbdry;
        bool m_init;
        vector_int m_nvars;
        vector_int m_loc;
        vector_int m_mask;
        int m_pts;
        doublereal m_solve_time;

        int m_ss_jac_age, m_ts_jac_age;

        // stats
        int m_nevals;
        doublereal m_evaltime;
        vector_int m_gridpts;
        vector_int m_jacEvals;
        vector_fp m_jacElapsed;
        vector_int m_funcEvals;
        vector_fp m_funcElapsed;

    private:
    };

}

#endif


