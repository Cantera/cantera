//! @file Wall.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/Wall.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/numerics/Func1.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/base/stringUtils.h"

namespace Cantera
{
Wall::Wall() : m_left(0), m_right(0),
    m_surf(2),
    m_area(1.0), m_k(0.0), m_rrth(0.0), m_emiss(0.0),
    m_vf(0), m_qf(0)
{
}

bool Wall::install(ReactorBase& rleft, ReactorBase& rright)
{
    // check if wall is already installed
    if (m_left || m_right) {
        return false;
    }
    m_left =  &rleft;
    m_right = &rright;
    m_left->addWall(*this, 0);
    m_right->addWall(*this, 1);
    m_surf[0].setReactor(&rleft);
    m_surf[1].setReactor(&rright);
    return true;
}

doublereal Wall::vdot(doublereal t)
{
    double rate1 = m_k * m_area * (m_left->pressure() - m_right->pressure());
    if (m_vf) {
        rate1 += m_area * m_vf->eval(t);
    }
    return rate1;
}

doublereal Wall::Q(doublereal t)
{
    double q1 = (m_area * m_rrth) *
                (m_left->temperature() - m_right->temperature());
    if (m_emiss > 0.0) {
        double tl = m_left->temperature();
        double tr = m_right->temperature();
        q1 += m_emiss * m_area * StefanBoltz * (tl*tl*tl*tl - tr*tr*tr*tr);
    }
    if (m_qf) {
        q1 += m_area * m_qf->eval(t);
    }
    return q1;
}

}
