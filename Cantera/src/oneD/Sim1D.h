/**
 * @file Sim1D.h
 */

#ifndef CT_SIM1D_H
#define CT_SIM1D_H

#include "OneDim.h"
#include "../funcs.h"

namespace Cantera {

    /**
     * One-dimensional simulations.
     */
    class Sim1D : public OneDim {

    public: 

        /**
         * Default constructor. This constructor can be used to create
         * a dummy object if necessary, but is not usually called in
         * user programs. Use the next constructor instead.
         */
        Sim1D();


        /**
         * Standard constructor. 
         * @param domains A vector of pointers to the domains to be linked together.
         * The domains must appear in left-to-right order.
         */ 
        Sim1D(vector<Resid1D*>& domains);

        /// Destructor. Does nothing.
        virtual ~Sim1D(){}


        /// Set one entry in the solution vector.
        void setValue(int dom, int comp, int localPoint,  doublereal value);

        /// Get one entry in the solution vector.
        doublereal value(int dom, int comp, int localPoint) const;

        /// Specify a profile for one component of one domain.
        void setProfile(int dom, int comp, const vector_fp& pos, 
            const vector_fp& values);

        /// Set component 'comp' of domain 'dom' to value 'v' at all points.
        void setFlatProfile(int dom, int comp, doublereal v);

        /// Print to stream s the current solution for all domains.     
        void showSolution(ostream& s);

        /// Calls method _finalize in each domain.
        void finalize();

        void setTimeStep(doublereal stepsize, int n, integer* tsteps);

        //void setMaxTimeStep(doublereal tmax) { m_maxtimestep = tmax; }

        void solve(int loglevel = 0, bool refine_grid = true);

        int refine(int loglevel=0);

        void setRefineCriteria(int dom = -1, doublereal ratio = 10.0,
            doublereal slope = 0.8, doublereal curve = 0.8);

    protected:

        vector_fp m_x;           // the solution vector
        vector_fp m_xnew;       // a work array used to hold the residual 
                                // or the new solution
        doublereal m_tstep;     // timestep
        vector_int m_steps;     // array of number of steps to take before 
                                // re-attempting the steady-state solution

    private:

        void newtonSolve(int loglevel);


    };

}

#endif


