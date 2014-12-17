/**
 * @file Group.cpp
 *
 *  Implementation file for the Group class used in reaction path analysis.
 */

// Copyright 2001  California Institute of Technology

// reaction path analysis support

#include "cantera/kinetics/Group.h"
#include <iostream>

namespace Cantera
{

void Group::validate()
{

    size_t n = m_comp.size();

    // if already checked and not valid, return
    if (m_sign == -999) {
        return;
    }

    m_sign = 0;
    bool ok = true;
    for (size_t m = 0; m < n; m++) {
        if (m_comp[m] != 0) {
            if (m_sign == 0) {
                m_sign = m_comp[m]/abs(m_comp[m]);
            } else if (m_sign * m_comp[m] < 0) {
                ok = false;
                break;
            }
        }
    }
    if (!ok) {
        m_sign = -999;
        m_comp.resize(n,0);
    }
}

std::ostream& Group::fmt(std::ostream& s,
                         const std::vector<std::string>& esymbols) const
{
    s << "(";
    int nm;
    bool first = true;
    size_t n = m_comp.size();
    for (size_t m = 0; m < n; m++) {
        nm = m_comp[m];
        if (nm != 0) {
            if (!first) {
                s << "-";
            }
            s << esymbols[m];
            if (nm != 1) {
                s << nm;
            }
            first = false;
        }
    }
    s << ")";
    return s;
}

std::ostream& operator<<(std::ostream& s, const Cantera::Group& g)
{
    if (g.valid()) {
        s << g.m_comp;
    } else {
        s << "<none>";
    }
    return s;
}

}
