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
 * specified as a function of time.
 *
 * The device coefficient *c* is specifies the mass flow coefficient with units of kg/s.
 * The mass flow rate is computed as:
 * @f[\dot{m} = m g(t) @f]
 * where *g* is a function of time that is set by `setTimeFunction`.
 * If no function is specified, the mass flow rate defaults to:
 * @f[\dot{m} = m @f]
 *
 * @ingroup connectorGroup
 */
class MassFlowController : public FlowDevice
{
public:
    using FlowDevice::FlowDevice;  // inherit constructors

    string type() const override {
        return "MassFlowController";
    }

    void setMassFlowRate(double mdot) override;

    //! Set the mass flow coefficient.
    /*!
     * *m* has units of kg/s. The mass flow rate is computed as:
     * @f[\dot{m} = m g(t) @f]
     * where *g* is a function of time that is set by `setTimeFunction`.
     * If no function is specified, the mass flow rate defaults to:
     * @f[\dot{m} = m @f]
     */
    void setMassFlowCoeff(double m) {
        m_coeff = m;
    }

    //! Get the mass flow coefficient.
    double getMassFlowCoeff() {
        return m_coeff;
    }

    void setPressureFunction(Func1* f) override {
        throw NotImplementedError("MassFlowController::setPressureFunction");
    }

    void setPressureFunction(shared_ptr<Func1> f) override {
        throw NotImplementedError("MassFlowController::setPressureFunction");
    }

    //! If a function of time has been specified for mdot, then update the
    //! stored mass flow rate. Otherwise, mdot is a constant, and does not
    //! need updating.
    void updateMassFlowRate(double time) override;
};

/**
 * A class for flow controllers where the flow rate is equal to the flow rate
 * of a primary mass flow controller plus a correction proportional to the
 * pressure difference between the inlet and outlet.
 *
 * The device coefficient *c* sets the proportionality constant between pressure drop
 * and mass flow rate and has units of kg/s/Pa. The mass flow rate is computed as:
 * @f[\dot{m} = \dot{m}_{primary} + c f(\Delta P) @f]
 * where *f* is a functions of pressure drop that is set by
 * `setPressureFunction`. If no functions is specified, the mass flow
 * rate defaults to:
 * @f[\dot{m} = \dot{m}_{primary} + c \Delta P @f]
 *
 * @ingroup connectorGroup
 */
class PressureController : public FlowDevice
{
public:
    using FlowDevice::FlowDevice;  // inherit constructors

    string type() const override {
        return "PressureController";
    }

    bool ready() override {
        return FlowDevice::ready() && m_primary != 0;
    }

    void setPrimary(shared_ptr<ConnectorNode> primary) override {
        if (!std::dynamic_pointer_cast<MassFlowController>(primary)) {
            throw CanteraError("PressureController::setPrimary",
                               "Invalid primary mass flow controller with type {}.",
                               primary->type());
        }
        auto dev = std::dynamic_pointer_cast<FlowDevice>(primary);
        m_primary = dev.get();
    }

    //! Set the primary mass flow controller.
    //! @since New in %Cantera 3.0.
    //! @deprecated  To be removed after %Cantera 3.2. Replaceable by version using
    //!     shared pointer.
    void setPrimary(FlowDevice* primary) {
        warn_deprecated("PressureController::setPrimary",
                        "To be removed after %Cantera 3.2. Replaceable by version "
                        "using shared pointer.");
        m_primary = primary;
    }

    void setTimeFunction(Func1* g) override {
        throw NotImplementedError("PressureController::setTimeFunction");
    }

    void setTimeFunction(shared_ptr<Func1> f) override {
        throw NotImplementedError("PressureController::setTimeFunction");
    }

    //! Set the proportionality constant between pressure drop and mass flow rate.
    /*!
     * *c* has units of kg/s/Pa. The mass flow rate is computed as:
     * @f[\dot{m} = \dot{m}_{primary} + c f(\Delta P) @f]
     * where *f* is a functions of pressure drop that is set by
     * `setPressureFunction`. If no functions is specified, the mass flow
     * rate defaults to:
     * @f[\dot{m} = \dot{m}_{primary} + c \Delta P @f]
     */
    void setPressureCoeff(double c) {
        m_coeff = c;
    }

    //! Get the pressure coefficient.
    double getPressureCoeff() {
        return m_coeff;
    }

    void updateMassFlowRate(double time) override;

protected:
    FlowDevice* m_primary = nullptr;
};

//! Supply a mass flow rate that is a function of the pressure drop across the valve.
/*!
 * The default behavior is a linearly proportional to the pressure difference.
 * Note that real valves do not have this behavior, so this class does not
 * model real, physical valves.
 *
 * The device coefficient *c* sets the proportionality constant between pressure drop
 * and mass flow rate and has units of kg/s/Pa. The mass flow rate is computed as:
 * @f[\dot{m} = c g(t) f(\Delta P) @f]
 * where *g* and *f* are functions of time and pressure drop that are set
 * by `setTimeFunction` and `setPressureFunction`, respectively. If no functions are
 * specified, the mass flow rate defaults to:
 * @f[\dot{m} = c \Delta P @f]
 *
 * @ingroup connectorGroup
 */
class Valve : public FlowDevice
{
public:
    using FlowDevice::FlowDevice;  // inherit constructors

    string type() const override {
        return "Valve";
    }

    //! Set the proportionality constant between pressure drop and mass flow rate.
    /*!
     * *c* has units of kg/s/Pa. The mass flow rate is computed as:
     * @f[\dot{m} = c g(t) f(\Delta P) @f]
     * where *g* and *f* are functions of time and pressure drop that are set
     * by `setTimeFunction` and `setPressureFunction`, respectively. If no functions are
     * specified, the mass flow rate defaults to:
     * @f[\dot{m} = c \Delta P @f]
     */
    void setValveCoeff(double c) {
        m_coeff = c;
    }

    //! Get the valve coefficient.
    double getValveCoeff() {
        return m_coeff;
    }

    //! Compute the current mass flow rate, based on the pressure difference.
    void updateMassFlowRate(double time) override;
};

}
#endif
