//! @file flowControllers.h Some flow devices derived from class FlowDevice.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_FLOWCONTR_H
#define CT_FLOWCONTR_H

#include "FlowDevice.h"
#include "cantera/base/ctexceptions.h"

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

    virtual std::string typeStr() const {
        return "MassFlowController";
    }

    virtual bool ready() {
        return FlowDevice::ready() && m_mdot >= 0.0;
    }

    //! Set the fixed mass flow rate (kg/s) through the mass flow controller.
    void setMassFlowRate(double mdot);

    //! Set the mass flow coefficient.
    /*!
     * *m* has units of kg/s. The mass flow rate is computed as:
     * \f[\dot{m} = m g(t) \f]
     * where *g* is a function of time that is set by `setTimeFunction`.
     * If no function is specified, the mass flow rate defaults to:
     * \f[\dot{m} = m \f]
     */
    void setMassFlowCoeff(double m) {
        m_coeff = m;
    }

    //! Get the mass flow coefficient.
    double getMassFlowCoeff() {
        return m_coeff;
    }

    virtual void setPressureFunction(Func1* f) {
        throw NotImplementedError("MassFlowController::setPressureFunction");
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

    virtual std::string typeStr() const {
        return "PressureController";
    }

    virtual bool ready() {
        return FlowDevice::ready() && m_master != 0;
    }

    void setMaster(FlowDevice* master) {
        m_master = master;
    }

    virtual void setTimeFunction(Func1* g) {
        throw NotImplementedError("PressureController::setTimeFunction");
    }

    //! Set the proportionality constant between pressure drop and mass flow
    //! rate
    /*!
     * *c* has units of kg/s/Pa. The mass flow rate is computed as:
     * \f[\dot{m} = \dot{m}_{master} + c f(\Delta P) \f]
     * where *f* is a functions of pressure drop that is set by
     * `setPressureFunction`. If no functions is specified, the mass flow
     * rate defaults to:
     * \f[\dot{m} = \dot{m}_{master} + c \Delta P \f]
     */
    void setPressureCoeff(double c) {
        m_coeff = c;
    }

    //! Get the pressure coefficient.
    double getPressureCoeff() {
        return m_coeff;
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

    virtual std::string typeStr() const {
        return "Valve";
    }

    //! Set the proportionality constant between pressure drop and mass flow
    //! rate
    /*!
     * *c* has units of kg/s/Pa. The mass flow rate is computed as:
     * \f[\dot{m} = c \Delta P \f]
     */
    //! @deprecated To be removed after Cantera 2.5.
    void setPressureCoeff(double c) {
        warn_deprecated("Valve::setPressureCoeff",
                        "To be removed after Cantera 2.5. "
                        "Use Valve::setValveCoeff instead.");
        m_coeff = c;
    }

    //! Set the proportionality constant between pressure drop and mass flow
    //! rate
    /*!
     * *c* has units of kg/s/Pa. The mass flow rate is computed as:
     * \f[\dot{m} = c g(t) f(\Delta P) \f]
     * where *g* and *f* are functions of time and pressure drop that are set
     * by `setTimeFunction` and `setPressureFunction`, respectively. If no functions are
     * specified, the mass flow rate defaults to:
     * \f[\dot{m} = c \Delta P \f]
     */
    void setValveCoeff(double c) {
        m_coeff = c;
    }

    //! Get the valve coefficient.
    double getValveCoeff() {
        return m_coeff;
    }

    /// Compute the currrent mass flow rate, based on the pressure difference.
    virtual void updateMassFlowRate(double time);
};

}
#endif
