/**
 *  @file Reactor.h
 *
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_REACTOR_H
#define CT_REACTOR_H

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ReactorBase.h"
#include "../FuncEval.h"
#include "../CVode.h"
#include "../Kinetics.h"

#define INCL_REACTOR_INTEG

namespace Cantera {

    /**
     * Class Reactor is a general-purpose class for stirred
     * reactors. The reactor may have an arbitrary number of inlets
     * and outlets, each of which may be connected to a "flow device"
     * such as a mass flow controller, a pressure regulator,
     * etc. Additional reactors may be connected to the other end of
     * the flow device, allowing construction of arbitrary reactor
     * networks.
     *
     * The reactor class integrates the same governing equations no
     * mattter what type of reactor is simulated. The differences
     * among reactor types are completely specified by the attached
     * flow devices and the time-dependent user-specified boundary
     * conditions. 
     *
     * If an instance of class Reactor is used directly, it will
     * simulate an adiabatic, constant volume reactor with gas-phase
     * chemistry but no surface chemistry. Other reactor types may be
     * simulated by deriving a class from Reactor and overloading
     * method getParams.  This method allows specifying the following
     * in terms of the instantaneous reactor state:
     *
     *  - rate of change of the total volume (m^3/s) 
     *  - surface heat loss rate (W) 
     *  - species surface production rates (kmol/s)
     * 
     * class Reactor inherits from both ReactorBase and
     * FuncEval. ReactorBase provides the basic reactor-like methods
     * that FlowDevice instances can access to determine their mass
     * flow rate. Class FuncEval is the class used to define a system
     * of ODE's to be integrated.
     */

    class Reactor : public ReactorBase, public FuncEval {

    public:

        /**
         * Default constructor.
         */
        Reactor();


        /**
         * Destructor. Deletes the integrator.
         */
        virtual ~Reactor(){ 
#ifdef INCL_REACTOR_INTEG
            delete m_integ; 
#endif
}
        
        virtual int type() const { return ReactorType; }

        /** 
         * Advance the state of the reactor in time. On the first
         * call, internal method 'initialize' is called, and the maximum
         * integrator step size is set. By default, this is set to
         * 'time'. To specify a different maximum step size, precede the
         * call to advance with a call to setMaxStep. Note that this
         * cannot be reset after advance has been called.
         * 
         * @param time Final time (s).
         */
        virtual void advance(doublereal time) {
#ifdef INCL_REACTOR_INTEG
            if (!m_init) {
                setMaxStep(time);
                initialize();
            }
            m_integ->integrate(time);
            m_time = time;
            updateState(m_integ->solution());
            m_mix->saveState(m_state);
#else
            throw CanteraError("Reactor::advance",
                "Reactor::advance is deprecated. Use ReactorNet::advance");
#endif
        }

        virtual double step(doublereal time) {
#ifdef INCL_REACTOR_INTEG
            if (!m_init) {
                setMaxStep(time);
                initialize();
            }
            m_time = m_integ->step(time);
            updateState(m_integ->solution());
            m_mix->saveState(m_state);
            return m_time;
#else
            throw CanteraError("Reactor::step",
                "Reactor::step is deprecated. Use ReactorNet::step");
#endif
        }

        /**
         * Insert something into the reactor. The 'something' must
         * belong to a class that is a subclass of both ThermoPhase
         * and Kinetics.
         */
        template<class G>
        void insert(G& contents) {
            setThermoMgr(contents);
            setKineticsMgr(contents);
        }

        void setKineticsMgr(Kinetics& kin) {
            m_kin = &kin;
            if (m_kin->nReactions() == 0) disableChemistry();
        }

        /**
         * Set the maximum step size for integration.
         */
        void setMaxStep(doublereal maxstep) {
            m_maxstep = maxstep;
        }

        void disableChemistry() { m_chem = false; }
        void enableChemistry() { m_chem = true; }

        /// Set the energy equation on or off.
        void setEnergy(int eflag = 1) { 
            if (eflag > 0) m_energy = true;
            else m_energy = false;
        } 

        //-----------------------------------------------------

        /** @name References to internal objects */
        //@{

        /// Return a reference to the integrator.
        Integrator& integrator() { return *m_integ; }

        //@}


        //-----------------------------------------------------

        // overloaded methods of class FuncEval
        virtual int neq() { return m_nv; }
	virtual void eval(doublereal t, doublereal* y, doublereal* ydot);
        virtual void getInitialConditions(doublereal t0, size_t leny, 
            doublereal* y);

        //-----------------------------------------------------

        virtual void initialize(doublereal t0 = 0.0);
	void evalEqs(doublereal t, doublereal* y, doublereal* ydot);

        /**
         * Set the mixture to a state consistent with solution
         * vector y.
         */

        virtual void updateState(doublereal* y);

    protected:
        
        Kinetics*   m_kin;

        Integrator* m_integ;         // pointer to integrator
        doublereal m_temp_atol;      // tolerance on T
        doublereal m_maxstep;        // max step size
        doublereal m_vdot, m_Q;
        vector_fp m_atol;
        doublereal m_rtol;
        vector_fp m_work;
        vector_fp m_sdot;            // surface production rates
        bool m_chem;
        bool m_energy;
        int m_nv;

    private:
    };
}

#endif

