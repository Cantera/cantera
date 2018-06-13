//! @file FlowDevice.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_FLOWDEVICE_H
#define CT_FLOWDEVICE_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/global.h"
#include "cantera/base/stringUtils.h"

namespace Cantera
{
class Func1;
class ReactorBase;

const int MFC_Type = 1;
const int PressureController_Type = 2;
const int Valve_Type = 3;

/**
 * Base class for 'flow devices' (valves, pressure regulators, etc.)
 * connecting reactors.
 * @ingroup reactor0
 */
class FlowDevice
{
public:
    FlowDevice() : m_mdot(0.0), m_func(0), m_type(0),
        m_nspin(0), m_nspout(0),
        m_in(0), m_out(0) {}

    virtual ~FlowDevice() {}
    FlowDevice(const FlowDevice&) = delete;
    FlowDevice& operator=(const FlowDevice&) = delete;

    //! Return an integer indicating the type of flow device
    int type() {
        return m_type;
    }

    //! Mass flow rate (kg/s).
    doublereal massFlowRate(double time = -999.0) {
        if (time != -999.0) {
            updateMassFlowRate(time);
        }
        return m_mdot;
    }

    //! Update the mass flow rate at time 'time'. This must be overloaded in
    //! subclassess to update m_mdot.
    virtual void updateMassFlowRate(doublereal time) {}

    //! Mass flow rate (kg/s) of outlet species k. Returns zero if this species
    //! is not present in the upstream mixture.
    doublereal outletSpeciesMassFlowRate(size_t k);

    //! specific enthalpy
    doublereal enthalpy_mass();

    //! Install a flow device between two reactors.
    /*!
     * @param in Upstream reactor.
     * @param out Downstream reactor.
     */
    bool install(ReactorBase& in, ReactorBase& out);

    virtual bool ready() {
        return (m_in != 0 && m_out != 0);
    }

    //! Return a reference to the upstream reactor.
    ReactorBase& in() const {
        return *m_in;
    }

    //! Return a const reference to the downstream reactor.
    const ReactorBase& out() const {
        return *m_out;
    }

    //! set parameters. Generic function used only in the Matlab interface. From
    //! Python or C++, device-specific functions like Valve::setPressureCoeff
    //! should be used instead.
    virtual void setParameters(int n, const double* coeffs) {
        m_coeffs.resize(n);
        std::copy(coeffs, coeffs + n, m_coeffs.begin());
    }

    //! Set a function of a single variable that is used in determining the
    //! mass flow rate through the device. The meaning of this function
    //! depends on the parameterization of the derived type.
    void setFunction(Func1* f);

    //! Set the fixed mass flow rate (kg/s) through the flow device.
    void setMassFlowRate(doublereal mdot) {
        m_mdot = mdot;
    }

protected:
    doublereal m_mdot;
    Func1* m_func;
    vector_fp m_coeffs;
    int m_type;

private:
    size_t m_nspin, m_nspout;
    ReactorBase* m_in;
    ReactorBase* m_out;
    std::vector<size_t> m_in2out, m_out2in;
};

}

#endif
