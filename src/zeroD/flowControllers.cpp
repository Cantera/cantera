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

void MassFlowController::updateMassFlowRate(double time)
{
    if (m_func) {
        m_mdot = m_func->eval(time);
    }
    m_mdot = std::max(m_mdot, 0.0);
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
    m_mdot = m_master->massFlowRate(time)
      + m_coeffs[0]*(in().pressure() - out().pressure());
    m_mdot = std::max(m_mdot, 0.0);
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
    double delta_P = in().pressure() - out().pressure();
    if (m_func) {
        m_mdot = m_func->eval(delta_P);
    } else {
        m_mdot = m_coeffs[0]*delta_P;
    }
    m_mdot = std::max(m_mdot, 0.0);
}
  
}
