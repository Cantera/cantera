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
 * @ingroup flowDeviceGroup
 */
class FlowDevice
{
public:
    FlowDevice(const string& name="(none)") : m_name(name) {}

    virtual ~FlowDevice() = default;
    FlowDevice(const FlowDevice&) = delete;
    FlowDevice& operator=(const FlowDevice&) = delete;

    //! String indicating the flow device implemented. Usually
    //! corresponds to the name of the derived class.
    virtual string type() const {
        return "FlowDevice";
    }

    //! Retrieve flow device name.
    string name() const {
        return m_name;
    }

    //! Set flow device name.
    void setName(const string& name) {
        m_name = name;
    }

    //! Set the default name of a flow device. Returns `false` if it was previously set.
    bool setDefaultName(map<string, int>& counts);

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
    //! subclasses to update m_mdot.
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

    //! Set a function of pressure that is used in determining the
    //! mass flow rate through the device. The evaluation of mass flow
    //! depends on the derived flow device class.
    virtual void setPressureFunction(Func1* f);

    //! Return current value of the time function.
    /*!
     * The mass flow rate [kg/s] is calculated for a Flow device, and multiplied by a
     * function of time, which returns 1.0 unless a user-defined function is provided.
     * The calculation of mass flow rate depends on the flow device.
     * @since New in %Cantera 3.0.
     */
    double evalTimeFunction();

    //! Set a function of time that is used in determining
    //! the mass flow rate through the device. The evaluation of mass flow
    //! depends on the derived flow device class.
    virtual void setTimeFunction(Func1* g);

    //! Set current reactor network time
    /*!
     * @since New in %Cantera 3.0.
     */
    void setSimTime(double time) {
        m_time = time;
    }

protected:
    string m_name;  //!< Flow device name.
    bool m_defaultNameSet = false;  //!< `true` if default name has been previously set.

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
