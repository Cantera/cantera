/**
 *  @file ThirdBodyMgr.h
 *
 * $Author: dggoodwin $
 * $Revision: 1.3 $
 * $Date: 2005/11/22 17:59:04 $
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_THIRDBODY_MGR_H
#define CT_THIRDBODY_MGR_H

#include <algorithm>

#include "ct_defs.h"
#include "utilities.h"
#include "Enhanced3BConc.h"


namespace Cantera {

    template<class _E>
    class ThirdBodyMgr{

    public:

        ThirdBodyMgr<_E>() : m_n(0) {}

        void install( int rxnNumber, const map<int, doublereal>& enhanced,
            doublereal dflt=1.0) {
            m_n++;
            m_reaction_index.push_back( rxnNumber );
            m_concm.push_back( _E(static_cast<int>(enhanced.size()),
			                      enhanced, dflt ) );
        }
        
        void update(const vector_fp& conc, doublereal ctot, workPtr work) {
            TYPENAME_KEYWORD vector<_E>::const_iterator b = m_concm.begin();
            //doublereal* v = m_values.begin();
            for (; b != m_concm.end(); ++b, ++work) 
                *work = b->update(conc, ctot);
        }

        void multiply(doublereal* output, const_workPtr work) {
            scatter_mult(work, work + m_n, 
                output, m_reaction_index.begin());
        }

        size_t workSize() { return m_concm.size(); }
        bool contains(int rxnNumber) {
            return (find(m_reaction_index.begin(), 
                        m_reaction_index.end(), rxnNumber) 
                != m_reaction_index.end());
        }
        
    protected:
        
        int m_n;
        vector_int      m_reaction_index;
        vector<_E>      m_concm;
    };
    
}
 
#endif
























































































































