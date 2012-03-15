/**
 *  @file ThirdBodyMgr.h
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_THIRDBODY_MGR_H
#define CT_THIRDBODY_MGR_H

#include <algorithm>

#include "cantera/base/ct_defs.h"
#include "cantera/base/utilities.h"
#include "Enhanced3BConc.h"


namespace Cantera
{

template<class _E>
class ThirdBodyMgr
{

public:

    ThirdBodyMgr<_E>() : m_n(0) {}

    void install(size_t rxnNumber, const std::map<size_t, doublereal>& enhanced,
                 doublereal dflt=1.0) {
        m_n++;
        m_reaction_index.push_back(rxnNumber);
        m_concm.push_back(_E(static_cast<int>(enhanced.size()),
                             enhanced, dflt));
    }

    void update(const vector_fp& conc, doublereal ctot, doublereal* work) {
        typename std::vector<_E>::const_iterator b = m_concm.begin();
        //doublereal* v = m_values.begin();
        for (; b != m_concm.end(); ++b, ++work) {
            *work = b->update(conc, ctot);
        }
    }

    void multiply(doublereal* output, const doublereal* work) {
        scatter_mult(work, work + m_n,
                     output, m_reaction_index.begin());
    }

    size_t workSize() {
        return m_concm.size();
    }
    bool contains(int rxnNumber) {
        return (find(m_reaction_index.begin(),
                     m_reaction_index.end(), rxnNumber)
                != m_reaction_index.end());
    }

protected:

    int m_n;
    std::vector<size_t> m_reaction_index;
    std::vector<_E>      m_concm;
};

}

#endif
























































































































