/**
 *  @file ThirdBodyMgr.h
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_THIRDBODY_MGR_H
#define CT_THIRDBODY_MGR_H

#include "cantera/base/utilities.h"
#include "Enhanced3BConc.h"


namespace Cantera
{

//! @deprecated Replaced by ThirdBodyCalc. To be removed after Cantera 2.2.
template<class _E>
class ThirdBodyMgr
{

public:
    ThirdBodyMgr() {
        warn_deprecated("class ThirdBodyMgr", "To be removed after Cantera 2.2.");
    }

    void install(size_t rxnNumber, const std::map<size_t, doublereal>& enhanced,
                 doublereal dflt=1.0) {
        m_reaction_index.push_back(rxnNumber);
        m_concm.push_back(_E(enhanced, dflt));
    }

    void update(const vector_fp& conc, doublereal ctot, doublereal* work) {
        typename std::vector<_E>::const_iterator b = m_concm.begin();
        for (; b != m_concm.end(); ++b, ++work) {
            *work = b->update(conc, ctot);
        }
    }

    void multiply(doublereal* output, const doublereal* work) {
        scatter_mult(work, work + m_reaction_index.size(),
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
    std::vector<size_t> m_reaction_index;
    std::vector<_E>      m_concm;
};

}

#endif
