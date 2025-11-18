//! @file FlowDevice.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_FLOWDEVICE_H
#define CT_FLOWDEVICE_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/global.h"
#include "cantera/base/ctexceptions.h"
#include "ConnectorNode.h"

namespace Cantera
{
class Func1;
class ReactorBase;

/**
 * Base class for 'flow devices' (valves, pressure regulators, etc.)
 * connecting reactors.
 * @ingroup connectorGroup
 */
class FlowDevice : public ConnectorNode
{
public:
    FlowDevice(shared_ptr<ReactorBase> r0, shared_ptr<ReactorBase> r1,
               const string& name="(none)");
    using ConnectorNode::ConnectorNode;  // inherit constructors

    string type() const override {
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

    //! Update the mass flow rate at time 'time'.
    //! This must be overloaded in derived classes to update m_mdot.
    virtual void updateMassFlowRate(double time) {}

    //! Set the fixed mass flow rate (kg/s) through a flow device.
    virtual void setMassFlowRate(double mdot) {
        throw NotImplementedError("FlowDevice::setMassFlowRate");
    }

    //! Get the device coefficient (defined by derived class).
    //! @since  New in %Cantera 3.2.
    double deviceCoefficient() const {
        return m_coeff;
    }

    //! Set the device coefficient (defined by derived class).
    //! @since  New in %Cantera 3.2.
    void setDeviceCoefficient(double c) {
        m_coeff = c;
    }

    //! Set the primary mass flow controller.
    //! @since New in %Cantera 3.2.
    virtual void setPrimary(shared_ptr<ConnectorNode> primary) {
        throw NotImplementedError("FlowDevice::setPrimary");
    }

    //! Mass flow rate (kg/s) of outlet species k. Returns zero if this species
    //! is not present in the upstream mixture.
    double outletSpeciesMassFlowRate(size_t k);

    //! specific enthalpy
    double enthalpy_mass();

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

    //! Return a mutable reference to the downstream reactor.
    ReactorBase& out() {
        return *m_out;
    }

    //! Return current value of the pressure function.
    /*!
     * The mass flow rate [kg/s] is calculated given the pressure drop [Pa] and a
     * coefficient set by a flow device specific function; unless a user-defined
     * pressure function is set, this is the pressure difference across the device.
     * The calculation of mass flow rate depends to the flow device.
     * @since New in %Cantera 3.0.
     */
    double evalPressureFunction();

    //! Set a function of pressure to modify the pressure response.
    //! Set a function of pressure that is used in determining the
    //! mass flow rate through the device. The evaluation of mass flow
    //! depends on the derived flow device class.
    //! @since  Changed in %Cantera 3.2. Previous version used a raw pointer.
    virtual void setPressureFunction(shared_ptr<Func1> f) {
        m_pfunc = f.get();
    }

    //! Return current value of the time function.
    /*!
     * The mass flow rate [kg/s] is calculated for a Flow device, and multiplied by a
     * function of time, which returns 1.0 unless a user-defined function is provided.
     * The calculation of mass flow rate depends on the flow device.
     * @since New in %Cantera 3.0.
     */
    double evalTimeFunction();

    //! Set a function of time to modulate the mass flow rate.
    //! Set a function of time that is used in determining
    //! the mass flow rate through the device. The evaluation of mass flow
    //! depends on the derived flow device class.
    //! @since  Changed in %Cantera 3.2. Previous version used a raw pointer.
    virtual void setTimeFunction(shared_ptr<Func1> g) {
        m_tfunc = g.get();
    }

    //! Set current reactor network time
    /*!
     * @since New in %Cantera 3.0.
     */
    void setSimTime(double time) {
        m_time = time;
    }

protected:
    double m_mdot = Undef;

    //! Function set by setPressureFunction; used by updateMassFlowRate
    Func1* m_pfunc = nullptr;

    //! Function set by setTimeFunction; used by updateMassFlowRate
    Func1* m_tfunc = nullptr;

    //! Coefficient set by derived classes; used by updateMassFlowRate
    double m_coeff = 1.0;

    //! Current reactor network time
    double m_time = 0.;

private:
    size_t m_nspin = 0;
    size_t m_nspout = 0;
    ReactorBase* m_in = nullptr;
    ReactorBase* m_out = nullptr;
    vector<size_t> m_in2out, m_out2in;
};

}

#endif
