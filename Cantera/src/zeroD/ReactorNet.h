/**
 *  @file ReactorNet.h
 */

/*
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2004  California Institute of Technology

#ifndef CT_REACTORNET_H
#define CT_REACTORNET_H

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "Reactor.h"
#include "../FuncEval.h"
#include "../CVode.h"

namespace CanteraZeroD {


    class ReactorNet : public FuncEval {

    public:

        ReactorNet();
        virtual ~ReactorNet(){}

        //-----------------------------------------------------

        /** @name Methods to set up a simulation. */
        //@{

        /**
         * Set initial time. Default = 0.0 s. Restarts integration
         * from this time using the current mixture state as the
         * initial condition.
         */
        void setInitialTime(doublereal time) {
            m_time = time;
            m_init = false;
        }

        /// Set the maximum time step.
        void setMaxTimeStep(double maxstep) {
            m_maxstep = maxstep;
            m_init = false;
        }
            
        void setTolerances(doublereal rtol, doublereal atol) {
            if (rtol >= 0.0) m_rtol = rtol;
            if (atol >= 0.0) m_atols = atol;
            m_init = false;
        }

        /// Current value of the simulation time.
        doublereal time() { return m_time; }

        /// Relative tolerance.
        doublereal rtol() { return m_rtol; }
        doublereal atol() { return m_atols; }
        
        /**
         * Initialize the reactor network. 
         */
        void initialize(doublereal t0 = 0.0);

        /**
         * Advance the state of all reactors in time.
         * @param time Time to advance to (s). 
         */
        void advance(doublereal time);

        double step(doublereal time);

        //@}

        void addReactor(ReactorBase* r) {
            if (r->type() == ReactorType) {
                m_r.push_back(r);
                m_nr++;
            }
        }

        ReactorBase& reactor(int n) {
            return *m_r[n];
        }


        /// Return a reference to the integrator.
        Integrator& integrator() { return *m_integ; }

        void updateState(doublereal* y);

        //-----------------------------------------------------

        // overloaded methods of class FuncEval
        virtual int neq() { return m_nv; }
	virtual void eval(doublereal t, doublereal* y, doublereal* ydot);
        virtual void getInitialConditions(doublereal t0, size_t leny, 
            doublereal* y);


    protected:

        vector<ReactorBase*> m_r;
        vector<Reactor*> m_reactors;
        int m_nr;
        int m_nreactors;
        Integrator* m_integ;
        doublereal m_time;
        bool m_init;
        int m_nv;
        vector_int m_size;
        vector_fp m_atol;
        doublereal m_rtol;
        doublereal m_atols;
        doublereal m_maxstep;
        bool m_verbose;

    private:

    };
}

#endif

