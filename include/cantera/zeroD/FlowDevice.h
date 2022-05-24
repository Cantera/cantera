//! @file FlowDevice.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_FLOWDEVICE_H
#define CT_FLOWDEVICE_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/global.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{
class Func1;
class ReactorBase;

/**
 * Base class for 'flow devices' (valves, pressure regulators, etc.)
 * connecting reactors.
 * @ingroup ZeroD
 */
class FlowDevice
{
public:
    FlowDevice();

    virtual ~FlowDevice() {}
    FlowDevice(const FlowDevice&) = delete;
    FlowDevice& operator=(const FlowDevice&) = delete;

    //! String indicating the flow device implemented. Usually
    //! corresponds to the name of the derived class.
    virtual std::string type() const {
        return "FlowDevice";
    }

    //! Mass flow rate (kg/s).
    double massFlowRate() {
        if (m_mdot == Undef) {
            throw CanteraError("FlowDevice::massFlowRate",
                               "Flow device is not ready. Try initializing the reactor network.");
        } else {
            return m_mdot;
        }
    }

    //! Update the mass flow rate at time 'time'. This must be overloaded in
    //! subclassess to update m_mdot.
    virtual void updateMassFlowRate(double time) {}

    //! Mass flow rate (kg/s) of outlet species k. Returns zero if this species
    //! is not present in the upstream mixture.
    double outletSpeciesMassFlowRate(size_t k);

    //! specific enthalpy
    double enthalpy_mass();

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

    //! Set a function of pressure that is used in determining the
    //! mass flow rate through the device. The evaluation of mass flow
    //! depends on the derived flow device class.
    virtual void setPressureFunction(Func1* f);

    //! Set a function of time that is used in determining
    //! the mass flow rate through the device. The evaluation of mass flow
    //! depends on the derived flow device class.
    virtual void setTimeFunction(Func1* g);

protected:
    double m_mdot;

    //! Function set by setPressureFunction; used by updateMassFlowRate
    Func1* m_pfunc;

    //! Function set by setTimeFunction; used by updateMassFlowRate
    Func1* m_tfunc;

    //! Coefficient set by derived classes; used by updateMassFlowRate
    double m_coeff;

private:
    size_t m_nspin, m_nspout;
    ReactorBase* m_in;
    ReactorBase* m_out;
    std::vector<size_t> m_in2out, m_out2in;
};

}

#endif
