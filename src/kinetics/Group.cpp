/**
 * @file Group.cpp Implementation file for the Group class used in reaction path
 *  analysis.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/Group.h"
#include <iostream>

namespace Cantera
{

void Group::validate()
{
    // if already checked and not valid, return
    if (m_sign == -999) {
        return;
    }

    m_sign = 0;
    bool ok = true;
    for (size_t m = 0; m < m_comp.size(); m++) {
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
    }
}

std::ostream& Group::fmt(std::ostream& s,
                         const std::vector<std::string>& esymbols) const
{
    s << "(";
    bool first = true;
    for (size_t m = 0; m < m_comp.size(); m++) {
        int nm = m_comp[m];
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

std::ostream& operator<<(std::ostream& s, const Group& g)
{
    if (g.valid()) {
        s << g.m_comp;
    } else {
        s << "<none>";
    }
    return s;
}

}
