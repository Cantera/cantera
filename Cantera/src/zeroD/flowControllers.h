/**
 *  @file FlowDevice.h
 *
 *  Some flow devices derived from class FlowDevice.
 *
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_FLOWCONTR_H
#define CT_FLOWCONTR_H

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "FlowDevice.h"
#include "ReactorBase.h"
#include "PID_Controller.h"
#include "../Func1.h"

namespace Cantera {


    ///////////////////////////////////////////////////////////////


    /** 
     * A base class for devices that do not use closed-loop control.
     * This is defined only for convenience, in order to overload
     * virtual methods of FlowDevice that print warnings with ones
     * that do nothing.
     */
    class NoController : public FlowDevice {
    public:

        NoController() {}
        virtual ~NoController() {}

        NoController(const NoController& a) 
            : FlowDevice(a) {}

        NoController& operator=(const NoController& a) {
            return *this;
        }

        // unneeded methods 
        virtual void update() {}
        virtual void reset() {}
        virtual bool setGains(int n, const doublereal* gains) {return true;}
        virtual bool getGains(int n, doublereal* gains) {return true;}
        virtual doublereal maxError() { return 0.0; }
        virtual doublereal setpoint() { return 0.0; }
        virtual void setSetpoint(doublereal mdot) { }
        virtual bool ready() {
            return FlowDevice::ready();
        }

    protected:
    private:
    };


    //////////////////////////////////////////////////////////


    /**
     * A class for mass flow controllers. The mass flow rate is constant,
     * independent of any other parameters.
     */
    class MassFlowController : public NoController {
    public:

        MassFlowController() {}
        virtual ~MassFlowController() {}

        MassFlowController(const MassFlowController& a) 
            : NoController(a) {}

        MassFlowController& operator=(const MassFlowController& a) {
            if (this == &a) return *this;
            m_mdot = a.m_mdot;
            return *this;
        }

        virtual doublereal setpoint() { return m_mdot; }
        virtual void setSetpoint(doublereal mdot) { m_mdot = mdot; }
        virtual bool ready() {
            return FlowDevice::ready() && m_mdot >= 0.0;
        }

    protected:
    private:
    };


    class UserValve : public NoController {
    public:

        UserValve() : m_func(0) {}
        virtual ~UserValve() {}

        UserValve(const UserValve& a) : NoController(a) {}

        UserValve& operator=(const UserValve& a) {
            if (this == &a) return *this;
            m_func = a.m_func;
            return *this;
        }

        virtual bool ready() {
            return FlowDevice::ready() && m_func != 0;
        }

        virtual void setFunction(Func1* f) { m_func = f; }

        virtual doublereal massFlowRate() {
            return m_func->eval(in().pressure() - out().pressure());
        }

    protected:
        Func1* m_func;

    private:
    };


    /**
     * A class for mass flow controllers. The mass flow rate is constant,
     * independent of any other parameters.
     */
    class Valve : public NoController {
    public:

        Valve() {}
        virtual ~Valve() {}

        Valve(const Valve& a) : NoController(a) {}

        Valve& operator=(const Valve& a) {
            if (this == &a) return *this;
            m_mdot = a.m_mdot;
            return *this;
        }

        virtual bool ready() {
            return FlowDevice::ready() && m_coeffs.size() >= 1;
        }

        virtual doublereal massFlowRate() {
            m_mdot = m_coeffs[0]* (in().pressure() - out().pressure());
            return (m_mdot > 0.0 ? m_mdot : 0.0);
        }

    protected:

    private:
    };



    /**
     * A PressureRegulator is a device that controls the pressure 
     * of the upstream reactor by regulating the mass flow rate.
     */ 
    class PressureRegulator : public FlowDevice  {

    public:

        PressureRegulator() {}
        virtual ~PressureRegulator() {}
        PressureRegulator(const PressureRegulator& p) : m_pid(p.m_pid) {}
        PressureRegulator& operator=(const PressureRegulator& p) {
            if (this == &p) return *this;
            m_pid = p.m_pid;
            return *this;
        }

        // overloaded virtual methods
        virtual void setSetpoint(doublereal p0) { m_pid.setpoint(-p0); }
        virtual doublereal setpoint() { return -m_pid.setpoint(); }
        virtual bool ready() { 
            return FlowDevice::ready() && m_pid.ready(); }
        virtual void reset() {
            m_pid.reset(in().time()-1.e-12, -in().pressure());
        }
        virtual void update() {
            m_pid.update(in().time(), -in().pressure());
        }

        virtual bool setGains(int n, const doublereal* gains) {
            return m_pid.setGains(n, gains);
        }

        virtual bool getGains(int n, doublereal* gains) {
            return m_pid.getGains(n, gains);
        }

        virtual doublereal maxError() { return m_pid.maxError(); }

        virtual doublereal massFlowRate() {
            m_mdot =  m_pid.output(-in().pressure());
            return m_mdot;
        }

    protected:

    private:
        PID_Controller m_pid;
    };
}
#endif

