//! @file FlowDevice.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_FLOWDEVICE_H
#define CT_FLOWDEVICE_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/global.h"
#include "cantera/base/stringUtils.h"

namespace Cantera
{
class Func1;
class ReactorBase;

//! Magic numbers
//! @deprecated To be removed after Cantera 2.5.
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
    FlowDevice();

    virtual ~FlowDevice() {}
    FlowDevice(const FlowDevice&) = delete;
    FlowDevice& operator=(const FlowDevice&) = delete;

    //! String indicating the flow device implemented. Usually
    //! corresponds to the name of the derived class.
    virtual std::string typeStr() const {
        return "FlowDevice";
    }

    //! Return an integer indicating the type of flow device
    /*!
     * @deprecated To be changed after Cantera 2.5.
     */
    virtual int type() const {
        warn_deprecated("FlowDevice::type",
                        "To be changed after Cantera 2.5. "
                        "Return string instead of magic number; use "
                        "FlowDevice::typeStr during transition.");
        return m_type;
    }

    //! Mass flow rate (kg/s).
    double massFlowRate(double time = -999.0) {
        if (time != -999.0) {
            updateMassFlowRate(time);
        }
        return m_mdot;
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

    //! Set parameters. Generic function formerly used in the Matlab interface.
    //! @deprecated To be removed after Cantera 2.5.
    virtual void setParameters(int n, const double* coeffs) {
        warn_deprecated("FlowDevice::setParameters",
                        "To be removed after Cantera 2.5. "
                        "Use device-specific functions (e.g. "
                        "Valve::setValveCoeff) instead.");
        m_coeff = coeffs[0]; // vectorized coefficients are not used
    }

    //! Set a function of a single variable that is used in determining the
    //! mass flow rate through the device. The meaning of this function
    //! depends on the parameterization of the derived type.
    //! @deprecated To be removed after Cantera 2.5.
    void setFunction(Func1* f) {
        warn_deprecated("FlowDevice::setFunction",
                        "To be removed after Cantera 2.5. "
                        "Use FlowDevice::setTimeFunction or "
                        "FlowDevice::setPressureFunction instead.");
        if (typeStr()=="MassFlowController") {
            setTimeFunction(f);
        } else if (typeStr()=="Valve") {
            setPressureFunction(f);
        }
    }

    //! Set a function of pressure that is used in determining the
    //! mass flow rate through the device. The evaluation of mass flow
    //! depends on the derived flow device class.
    virtual void setPressureFunction(Func1* f);

    //! Set a function of time that is used in determining
    //! the mass flow rate through the device. The evaluation of mass flow
    //! depends on the derived flow device class.
    virtual void setTimeFunction(Func1* g);

    //! Set the fixed mass flow rate (kg/s) through the flow device.
    //! @deprecated To be removed after Cantera 2.5.
    void setMassFlowRate(double mdot) {
        warn_deprecated("FlowDevice::setMassFlowRate",
                        "To be removed after Cantera 2.5. "
                        "Use device-specific functions (e.g. "
                        "MassFlowController::setMassFlowRate or "
                        "Valve::setValveCoeff) instead.");
        m_mdot = mdot;
    }

protected:
    double m_mdot;

    //! Function set by setPressureFunction; used by updateMassFlowRate
    Func1* m_pfunc;

    //! Function set by setTimeFunction; used by updateMassFlowRate
    Func1* m_tfunc;

    //! Coefficient set by derived classes; used by updateMassFlowRate
    double m_coeff;

    int m_type; //!< @deprecated To be removed after Cantera 2.5.

private:
    size_t m_nspin, m_nspout;
    ReactorBase* m_in;
    ReactorBase* m_out;
    std::vector<size_t> m_in2out, m_out2in;
};

}

#endif
