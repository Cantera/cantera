
#include "Wall.h"
#include "ReactorBase.h"
#include "Func1.h"

namespace Cantera {

        
    Wall::Wall() : m_left(0), m_right(0), m_area(0.0), m_k(0.0), m_rrth(0.0), 
                   m_vf(0), m_qf(0) {}

    bool Wall::install(ReactorBase& rleft, ReactorBase& rright) {
        // check if wall is already installed
        if (m_left || m_right) return false;
        m_left =  &rleft;
        m_right = &rright;
        m_left->addWall(*this, 1);
        m_right->addWall(*this, -1);
        return true; 
    }

    doublereal Wall::vdot(doublereal t) {
        double rate1 = m_k * m_area * 
                       (m_left->pressure() - m_right->pressure()); 
        if (m_vf) rate1 += m_vf->eval(t);
        return rate1;
    }

    doublereal Wall::Q(doublereal t) {
        double q1 = (m_area * m_rrth) *
                    (m_left->temperature() - m_right->temperature());
        if (m_qf) q1 += m_area * m_qf->eval(t);
        return q1;
    }

}        
