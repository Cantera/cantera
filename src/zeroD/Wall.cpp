//! @file Wall.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/stringUtils.h"
#include "cantera/numerics/Func1.h"
#include "cantera/zeroD/Wall.h"

namespace Cantera
{

WallBase::WallBase(shared_ptr<ReactorNode> r0, shared_ptr<ReactorNode> r1,
                   const string& name) : Connector(r0, r1, name)
{
    if (!m_nodes.first || !m_nodes.second) {
        warn_deprecated("FlowDevice::FlowDevice",
            "After Cantera 3.1, Reactors must be provided to a FlowDevice "
            "constructor.");
        return;
    }
    // todo: switch to shared pointers after Cantera 3.1.
    m_left = std::dynamic_pointer_cast<ReactorBase>(r0).get();
    m_right = std::dynamic_pointer_cast<ReactorBase>(r1).get();
    m_left->addWall(*this, 0);
    m_right->addWall(*this, 1);
}

bool WallBase::install(ReactorBase& rleft, ReactorBase& rright)
{
    warn_deprecated("WallBase::install",
        "To be removed after Cantera 3.1. Reactors should be provided to constructor "
        "instead.");
    // check if wall is already installed
    if (m_left || m_right) {
        return false;
    }
    m_left =  &rleft;
    m_right = &rright;
    m_left->addWall(*this, 0);
    m_right->addWall(*this, 1);
    return true;
}

void WallBase::setArea(double a) {
    m_area = a;
}

double Wall::velocity() const {
    if (m_vf) {
        return m_vf->eval(m_time);
    }
    return 0.;
}

double Wall::expansionRate()
{
    if (!ready()) {
        throw CanteraError("Wall::expansionRate",
                           "Wall is not ready; some parameters have not been set.");
    }
    double rate = m_k * m_area * (m_left->pressure() - m_right->pressure());

    if (m_vf) {
        rate += m_area * m_vf->eval(m_time);
    }
    return rate;
}

double Wall::heatFlux() const {
    if (m_qf) {
        return m_qf->eval(m_time);
    }
    return 0.;
}

double Wall::heatRate()
{
    if (!ready()) {
        throw CanteraError("Wall::heatRate",
                           "Wall is not ready; some parameters have not been set.");
    }
    double q1 = (m_area * m_rrth) *
                (m_left->temperature() - m_right->temperature());
    if (m_emiss > 0.0) {
        double tl = m_left->temperature();
        double tr = m_right->temperature();
        q1 += m_emiss * m_area * StefanBoltz * (tl*tl*tl*tl - tr*tr*tr*tr);
    }

    if (m_qf) {
        q1 += m_area * m_qf->eval(m_time);
    }
    return q1;
}

}
