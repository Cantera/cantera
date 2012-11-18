/**
 *  @file FlowDevice.h
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_FLOWDEVICE_H
#define CT_FLOWDEVICE_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/global.h"
#include "cantera/base/stringUtils.h"

namespace Cantera
{
class Func1;
class ReactorBase;  // forward reference

const int MFC_Type = 1;
const int PressureController_Type = 2;
const int Valve_Type = 3;

/**
 * Base class for 'flow devices' (valves, pressure regulators,
 * etc.)  connecting reactors. Allowance is made for devices that
 * are closed-loop controllers. Several methods for these are
 * defined here that do nothing but may be overloaded to set or
 * get the setpoint, gains, etc. The behavior of overloaded
 * methods should be consistent with the behavior described
 * here. The base-class versions of these methods print a warning
 * if called.
 * @ingroup reactor0
 */
class FlowDevice
{

public:

    /// Constructor
    FlowDevice() : m_mdot(0.0), m_func(0), m_type(0),
        m_nspin(0), m_nspout(0),
        m_in(0), m_out(0) {}

    /// Destructor (does nothing)
    virtual ~FlowDevice() {}

    //        /// Copy constructor.
    //         FlowDevice(const FlowDevice& a) : m_in(a.m_in), m_out(a.m_out) {}

    //         /// Assignment operator
    //         FlowDevice& operator=(const FlowDevice& a) {
    //             if (this == &a) return *this;
    //             m_in = a.m_in;
    //             m_out = a.m_out;
    //             return *this;
    //         }

    int type() {
        return m_type;
    }

    /**
     * Mass flow rate (kg/s).
     */
    doublereal massFlowRate(double time = -999.0) {
        if (time != -999.0) {
            updateMassFlowRate(time);
        }
        return m_mdot;
    }

    // Update the mass flow rate at time 'time'. This must be
    // overloaded in subclassess to update m_mdot.
    virtual void updateMassFlowRate(doublereal time) {}

    // mass flow rate of outlet species k
    doublereal outletSpeciesMassFlowRate(size_t k);

    // specific enthalpy
    doublereal enthalpy_mass();

    //         /**
    //          * Setpoint. Default = 0.0.
    //          */
    //         virtual doublereal setpoint() { warn("setpoint"); return 0.0; }


    //         /* Update the internal state, if necessary. By default this method
    //          * does nothing, but may be overloaded for devices that have a
    //          * state.
    //          */
    //         virtual void update() {warn("update");}


    //         /* Reset the device. By default this method does nothing, but
    //          * may be overloaded for devices that have a state that depends on
    //          * past history.
    //          */
    //         virtual void reset() {warn("reset");}

    //         /**
    //          * Set the setpoint. May be changed at any time. By default,
    //          * this does nothing.
    //          */
    //         virtual void setSetpoint(doublereal value) {warn("setSetpoint");}

    //         /**
    //          * Set the controller gains. Returns false if the number of
    //          * gains is too small, or if an illegal value is specified.
    //          */
    //         virtual bool setGains(int n, const doublereal* gains) {
    //             warn("setGains");
    //             return true;
    //         }

    //         /**
    //          * Get the controller gains. Returns false if the 'gains'
    //          * array is too small.
    //          */
    //         virtual bool getGains(int n, doublereal* gains) {
    //             warn("getGains");
    //             return true;
    //         }

    //         /**
    //          * Maximum difference between input and setpoint since
    //          * last call to 'reset'.
    //          */
    //         virtual doublereal maxError() {warn("maxError"); return 0.0;}

    /**
     * Install a flow device between two reactors.
     * @param in Upstream reactor.
     * @param out Downstream reactor.
     */
    bool install(ReactorBase& in, ReactorBase& out);

    virtual bool ready() {
        return (m_in != 0 && m_out != 0);
    }

    /// Return a reference to the upstream reactor.
    ReactorBase& in() const {
        return *m_in;
    }

    /// Return a const reference to the downstream reactor.
    const ReactorBase& out() const {
        return *m_out;
    }

    /// set parameters
    virtual void setParameters(int n, doublereal* coeffs) {
        m_coeffs.resize(n);
        std::copy(coeffs, coeffs + n, m_coeffs.begin());
    }

    void setFunction(Cantera::Func1* f);
    void setMassFlowRate(doublereal mdot) {
        m_mdot = mdot;
    }


protected:

    doublereal m_mdot;
    Cantera::Func1* m_func;
    vector_fp m_coeffs;
    int m_type;

private:

    size_t m_nspin, m_nspout;
    ReactorBase* m_in;
    ReactorBase* m_out;
    std::vector<size_t> m_in2out, m_out2in;

    void warn(std::string meth) {
        writelog(std::string("Warning: method ") + meth + " of base class "
                 + " FlowDevice called. Nothing done.\n");
    }
};

}

#endif
