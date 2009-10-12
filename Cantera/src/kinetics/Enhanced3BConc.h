/**
 *  @file Enhanced3BConc.h
 */

/*  $Author: dggoodwin $
 *  $Date: 2007/05/04 14:27:23 $
 *  $Revision: 1.1 $
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_ENH_CONC_H
#define CT_ENH_CONC_H

#include <map>

namespace Cantera {

    /**
     * Computes enhanced third-body concentrations. 
     * @see GasKinetics
     */
    class Enhanced3BConc {

    public:

        Enhanced3BConc() : m_n (0), m_deflt (1.0) {}

	Enhanced3BConc(int n, const std::map<int, doublereal>& enhanced, 
            doublereal deflt = 1.0) {
            std::map<int, doublereal>::const_iterator iter;
            for (iter = enhanced.begin(); iter != enhanced.end(); ++iter) {
                m_index.push_back( iter->first );
                m_eff.push_back( iter->second - deflt);
            }
            m_deflt = deflt;
            m_n = n;
	}

	Enhanced3BConc(int n, const vector_int& e_index, 
            const vector_fp& efficiencies, doublereal deflt = 1.0)
            : m_index (e_index), m_eff (efficiencies) {
            int i;
            m_n = n;
            m_deflt = deflt;
            for (i = 0; i < m_n; i++) {
                m_eff[i] -= m_deflt;
            }
	}

	doublereal update(const vector_fp& c, doublereal ctot) const {
            int i;
            doublereal sum = 0.0;
            for (i = 0; i < m_n; i++) {
                sum += m_eff[i] * c[m_index[i]];
            }
            return m_deflt * ctot  + sum;
	}
       
        void getEfficiencies(vector_fp& eff) const {
            int i;
            for (i = 0; i < m_n; i++) {
                eff[m_index[i]] = m_eff[i] + m_deflt;
            }	    
	}
        
    private:
	int         m_n;
	vector_int  m_index;
	vector_fp   m_eff;
	doublereal  m_deflt;
    };

}

#endif




