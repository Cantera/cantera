/**
 *  @file FlowReactor.h
 *
 * $Author: dggoodwin $
 * $Revision: 1.3 $
 * $Date: 2006/11/27 21:43:34 $
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_FLOWREACTOR_H
#define CT_FLOWREACTOR_H

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "Reactor.h"

namespace CanteraZeroD {

    /**
     * Adiabatic, reversible flow in a constant-area duct.
     */
    class FlowReactor : public Reactor {

    public:

        /**
         * Default constructor.
         */
        FlowReactor();

        /**
         * Destructor. 
         */
        virtual ~FlowReactor(){}
        
        virtual int type() const { return FlowReactorType; }

        //-----------------------------------------------------

        virtual void getInitialConditions(doublereal t0, size_t leny, 
            doublereal* y);


        //-----------------------------------------------------
        virtual int neq() { return m_nv; }

        virtual void initialize(doublereal t0 = 0.0);
	virtual void evalEqs(doublereal t, doublereal* y, 
            doublereal* ydot, doublereal* params);
        virtual void updateState(doublereal* y);

        void setMassFlowRate(doublereal mdot) {
            m_rho0 = m_thermo->density();
            m_speed = mdot/m_rho0;
            m_speed0 = m_speed;
            m_T = m_thermo->temperature();
            m_P0 = m_thermo->pressure() + m_rho0*m_speed*m_speed;
            m_h0 = m_thermo->enthalpy_mass() + 0.5*m_speed*m_speed;
        }

        void setTimeConstant(doublereal tau) {
            m_fctr = 1.0/tau;
        }

        double speed() const { return m_speed; }
        double distance() const { return m_dist; }
        virtual int componentIndex(std::string nm) const;

    protected:
        
        doublereal m_speed, m_dist, m_T;
        doublereal m_fctr;
        doublereal m_rho0, m_speed0, m_P0, m_h0;

    private:
    };
}

#endif

