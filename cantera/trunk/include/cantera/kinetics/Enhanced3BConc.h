/**
 *  @file Enhanced3BConc.h
 */
// Copyright 2001  California Institute of Technology


#ifndef CT_ENH_CONC_H
#define CT_ENH_CONC_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/global.h"

namespace Cantera
{

/**
 * Computes enhanced third-body concentrations.
 * @deprecated Replaced by ThirdBodyCalc. To be removed after Cantera 2.2.
 * @see GasKinetics
 */
class Enhanced3BConc
{

public:

    Enhanced3BConc() : m_deflt(1.0) {
        warn_deprecated("class Enhanced3BConc",
                        "To be removed after Cantera 2.2.");
    }

    Enhanced3BConc(const std::map<size_t, doublereal>& enhanced,
                   doublereal deflt = 1.0) {
        warn_deprecated("class Enhanced3BConc",
                        "To be removed after Cantera 2.2.");
        std::map<size_t, doublereal>::const_iterator iter;
        for (iter = enhanced.begin(); iter != enhanced.end(); ++iter) {
            m_index.push_back(iter->first);
            m_eff.push_back(iter->second - deflt);
        }
        m_deflt = deflt;
    }

    doublereal update(const vector_fp& c, doublereal ctot) const {
        doublereal sum = 0.0;
        for (size_t i = 0; i < m_eff.size(); i++) {
            sum += m_eff[i] * c[m_index[i]];
        }
        return m_deflt * ctot  + sum;
    }

    void getEfficiencies(vector_fp& eff) const {
        for (size_t i = 0; i < m_eff.size(); i++) {
            eff[m_index[i]] = m_eff[i] + m_deflt;
        }
    }

private:
    std::vector<size_t> m_index;
    vector_fp   m_eff;
    doublereal  m_deflt;
};

}

#endif
