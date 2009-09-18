/**
 * @file Group.cpp
 *
 *  Implementation file for the Group class used in reaction path analysis.
 *
 * $Author: hkmoffa $
 * $Revision: 1.2 $
 * $Date: 2008/12/29 21:34:08 $
 */

// Copyright 2001  California Institute of Technology


// reaction path analysis support

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "Group.h"

#include <algorithm>
#include <cmath>

namespace Cantera {
     
    /**
     * A group is 'valid' if all of its nonzero atom numbers have
     * the same sign, either positive or negative. This method
     * checks for this, and if the group is not valid it sets
     * m_sign to -999, and sets all atom numbers to zero.
     */
    void Group::validate() {

        int n = m_comp.size();

        // if already checked and not valid, return
        if (m_sign == -999) return;

        m_sign = 0;
        bool ok = true;
        for (int m = 0; m < n; m++) 
        {
            if (m_comp[m] != 0) 
            {
                if (m_sign == 0) { 
                    m_sign = m_comp[m]/abs(m_comp[m]);
                }
                else if (m_sign * m_comp[m] < 0) {
                    ok = false; break;
                }
            }
        }
        if (!ok) { m_sign = -999; m_comp.resize(n,0); }
    }

    std::ostream& Group::fmt(std::ostream& s, 
        const std::vector<std::string>& esymbols) const {
        s << "(";
        int nm;
        bool first = true;
        int n = m_comp.size();
        for (int m = 0; m < n; m++) {
            nm = m_comp[m];
            if (nm != 0) {
                if (!first) s << "-";
                s << esymbols[m];
                if (nm != 1) s << nm;
                first = false;
            }
        }
        s << ")";
        return s;
    }

    std::ostream& operator<<(std::ostream& s, const Cantera::Group& g) {
        if (g.valid()) {
 	   s << g.m_comp;
	} else {
	   s << "<none>";
        }
        return s;
    }

}
