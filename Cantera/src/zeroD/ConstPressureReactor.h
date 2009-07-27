/**
 *  @file Reactor.h
 *
 * $Author: dggoodwin $
 * $Revision: 1.3 $
 * $Date: 2006/11/27 21:43:34 $
 */

// Copyright 2001  California Institute of Technology
 
#ifndef CT_CONSTP_REACTOR_H
#define CT_CONSTP_REACTOR_H

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "Reactor.h"

namespace CanteraZeroD {

    /**
     * Class ConstPressureReactor is a class for constant-pressure 
     * reactors. The reactor may have an arbitrary number of inlets
     * and outlets, each of which may be connected to a "flow device"
     * such as a mass flow controller, a pressure regulator,
     * etc. Additional reactors may be connected to the other end of
     * the flow device, allowing construction of arbitrary reactor
     * networks.
     *
     */
    class ConstPressureReactor : public Reactor {

    public:

        /**
         * Default constructor.
         */
        ConstPressureReactor();

        /**
         * Destructor. Deletes the integrator.
         */
        virtual ~ConstPressureReactor(){}
        

        virtual int type() const { return ConstPressureReactorType; }

        //-----------------------------------------------------

        //virtual int neq() { return m_nv; }

        virtual void getInitialConditions(doublereal t0, size_t leny, 
            doublereal* y);

        virtual void initialize(doublereal t0 = 0.0);
	virtual void evalEqs(doublereal t, doublereal* y, 
            doublereal* ydot, doublereal* params);

        virtual void updateState(doublereal* y);

        virtual int componentIndex(std::string nm) const;

    protected:
        
    private:

    };
}

#endif

