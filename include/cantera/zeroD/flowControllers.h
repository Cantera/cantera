//! @file flowControllers.h Some flow devices derived from class FlowDevice.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_FLOWCONTR_H
#define CT_FLOWCONTR_H

#include "FlowDevice.h"

namespace Cantera
{

/**
 * A class for mass flow controllers. The mass flow rate is constant or
 * specified as a function of time..
 */
class MassFlowController : public FlowDevice
{
public:
    MassFlowController();

    virtual bool ready() {
        return FlowDevice::ready() && m_mdot >= 0.0;
    }

    /// If a function of time has been specified for mdot, then update the
    /// stored mass flow rate. Otherwise, mdot is a constant, and does not
    /// need updating.
    virtual void updateMassFlowRate(double time);
};

/**
 * A class for flow controllers where the flow rate is equal to the flow rate
 * of a "master" mass flow controller plus a correction proportional to the
 * pressure difference between the inlet and outlet.
 */
class PressureController : public FlowDevice
{
public:
    PressureController();

    virtual bool ready() {
        return FlowDevice::ready() && m_master != 0 && m_coeffs.size() == 1;
    }

    void setMaster(FlowDevice* master) {
        m_master = master;
    }

    //! Set the proportionality constant between pressure drop and mass flow
    //! rate
    /*!
     * *c* has units of kg/s/Pa. The mass flow rate is computed as:
     * \f[\dot{m} = \dot{m}_{master} + c \Delta P \f]
     */
    void setPressureCoeff(double c) {
        m_coeffs = {c};
    }

    virtual void updateMassFlowRate(double time);

protected:
    FlowDevice* m_master;
};

//! Supply a mass flow rate that is a function of the pressure drop across the
//! valve.
/*!
 * The default behavior is a linearly proportional to the pressure difference.
 * Note that real valves do not have this behavior, so this class does not
 * model real, physical valves.
 */
class Valve : public FlowDevice
{
public:
    Valve();

    virtual bool ready() {
        return FlowDevice::ready() && (m_coeffs.size() == 1 || m_func);
    }

    //! Set the proportionality constant between pressure drop and mass flow
    //! rate
    /*!
     * *c* has units of kg/s/Pa. The mass flow rate is computed as:
     * \f[\dot{m} = c \Delta P \f]
     */
    void setPressureCoeff(double c) {
        m_coeffs = {c};
    }

    /// Compute the currrent mass flow rate, based on the pressure difference.
    virtual void updateMassFlowRate(double time);
};

}
#endif
