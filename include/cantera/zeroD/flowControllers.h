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
 * @ingroup connectorGroup
 */
class MassFlowController : public FlowDevice
{
public:
    using FlowDevice::FlowDevice;  // inherit constructors

    //! Create a new MassFlowController.
    //! @param r0  Reactor upstream of the mass flow controller.
    //! @param r1  Reactor downstream of the mass flow controller.
    //! @param name  Name of the mass flow controller. Optional; if left empty, a
    //!     default name will be assigned when the object is integrated into a
    //!     ReactorNet.
    static shared_ptr<MassFlowController> create(
        shared_ptr<ReactorNode> r0, shared_ptr<ReactorNode> r1, const string& name="")
    {
        return shared_ptr<MassFlowController>( new MassFlowController(r0, r1, name) );
    }

    string type() const override {
        return "MassFlowController";
    }

    //! Set the fixed mass flow rate (kg/s) through the mass flow controller.
    void setMassFlowRate(double mdot);

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

    //! If a function of time has been specified for mdot, then update the
    //! stored mass flow rate. Otherwise, mdot is a constant, and does not
    //! need updating.
    void updateMassFlowRate(double time) override;
};

/**
 * A class for flow controllers where the flow rate is equal to the flow rate
 * of a primary mass flow controller plus a correction proportional to the
 * pressure difference between the inlet and outlet.
 * @ingroup connectorGroup
 */
class PressureController : public FlowDevice
{
public:
    using FlowDevice::FlowDevice;  // inherit constructors

    //! Create a new PressureController.
    //! @param r0  Reactor upstream of the pressure controller.
    //! @param r1  Reactor downstream of the pressure controller.
    //! @param name  Name of the pressure controller. Optional; if left empty, a default
    //!     name will be assigned when the object is integrated into a ReactorNet.
    static shared_ptr<PressureController> create(
        shared_ptr<ReactorNode> r0, shared_ptr<ReactorNode> r1, const string& name="")
    {
        return shared_ptr<PressureController>( new PressureController(r0, r1, name) );
    }

    string type() const override {
        return "PressureController";
    }

    bool ready() override {
        return FlowDevice::ready() && m_primary != 0;
    }

    //! Set the primary mass flow controller.
    /*!
     * @since New in %Cantera 3.0.
     */
    void setPrimary(FlowDevice* primary) {
        m_primary = primary;
    }

    void setTimeFunction(Func1* g) override {
        throw NotImplementedError("PressureController::setTimeFunction");
    }

    //! Set the proportionality constant between pressure drop and mass flow
    //! rate
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

//! Supply a mass flow rate that is a function of the pressure drop across the
//! valve.
/*!
 * The default behavior is a linearly proportional to the pressure difference.
 * Note that real valves do not have this behavior, so this class does not
 * model real, physical valves.
 * @ingroup connectorGroup
 */
class Valve : public FlowDevice
{
public:
    using FlowDevice::FlowDevice;  // inherit constructors

    //! Create a new Valve.
    //! @param r0  Reactor left of the valve.
    //! @param r1  Reactor right of the valve.
    //! @param name  Name of the valve. Optional; if left empty, a default name will be
    //!     assigned when the object is integrated into a ReactorNet.
    static shared_ptr<Valve> create(
        shared_ptr<ReactorNode> r0, shared_ptr<ReactorNode> r1, const string& name="")
    {
        return shared_ptr<Valve>( new Valve(r0, r1, name) );
    }

    string type() const override {
        return "Valve";
    }

    //! Set the proportionality constant between pressure drop and mass flow
    //! rate
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
