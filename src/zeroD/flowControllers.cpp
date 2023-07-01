//! @file flowControllers.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/flowControllers.h"
#include "cantera/zeroD/ReactorBase.h"
#include "cantera/numerics/Func1.h"

namespace Cantera
{

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

void PressureController::setMaster(FlowDevice* master)
{
    warn_deprecated("PressureController::setMaster",
        "To be removed after Cantera 3.0; replaced by setPrimary.");
    m_primary = master;
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
    m_primary->updateMassFlowRate(time);
    mdot += m_primary->massFlowRate();
    m_mdot = std::max(mdot, 0.0);
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
