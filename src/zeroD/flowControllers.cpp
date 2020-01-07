//! @file flowControllers.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/flowControllers.h"
#include "cantera/zeroD/ReactorBase.h"
#include "cantera/numerics/Func1.h"

namespace Cantera
{

MassFlowController::MassFlowController() : FlowDevice() {
    m_type = MFC_Type;
}

void MassFlowController::setMassFlowRate(double mdot)
{
    if (m_tfunc) {
        delete m_tfunc;
    }
    m_coeff = mdot;
}

void MassFlowController::updateMassFlowRate(double time)
{
    if (!ready()) {
        throw CanteraError("MassFlowController::updateMassFlowRate",
                           "Device is not ready; some parameters have not been set.");
    }
    double mdot = m_coeff;
    if (m_tfunc) {
        mdot *= m_tfunc->eval(time);
    }
    m_mdot = std::max(mdot, 0.0);
}

PressureController::PressureController() : FlowDevice(), m_master(0) {
    m_type = PressureController_Type;
}

void PressureController::updateMassFlowRate(double time)
{
    if (!ready()) {
        throw CanteraError("PressureController::updateMassFlowRate",
                           "Device is not ready; some parameters have not been set.");
    }
    double mdot = m_coeff;
    double delta_P = in().pressure() - out().pressure();
    if (m_pfunc) {
        mdot *= m_pfunc->eval(delta_P);
    } else {
        mdot *= delta_P;
    }
    mdot += m_master->massFlowRate(time);
    m_mdot = std::max(mdot, 0.0);
}

Valve::Valve() : FlowDevice() {
    m_type = Valve_Type;
}

void Valve::updateMassFlowRate(double time)
{
    if (!ready()) {
        throw CanteraError("Valve::updateMassFlowRate",
                           "Device is not ready; some parameters have not been set.");
    }
    double mdot = m_coeff;
    if (m_tfunc) {
        mdot *= m_tfunc->eval(time);
    }
    double delta_P = in().pressure() - out().pressure();
    if (m_pfunc) {
        mdot *= m_pfunc->eval(delta_P);
    } else {
        mdot *= delta_P;
    }
    m_mdot = std::max(mdot, 0.0);
}

}
