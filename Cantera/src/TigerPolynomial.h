/**
 *  @file TigerPolynomial.h
 *
 * $Author$
 * $Revision$
 * $Date$
 */

WARNING: May be out of date. Test before using!


// Copyright 2001  California Institute of Technology


#ifndef CT_TigerPOLY_H
#define CT_TigerPOLY_H

#include "TempFuncMgr.h"

namespace Cantera {

    /**
     * The Tiger polynomial parameterization. 
     */

    
    class TigerPolynomial {

    public:
      
        TigerPolynomial() : m_lowT(0.0), m_highT(0.0), m_Pref(0.0), 
            m_coeff(vector_fp(11)) {}

        /**
         * construct a Tiger polynomial parameterization from the
         * values in vector coeffs.
         * @param coeffs(0) minimum temperature (K)
         * @param coeffs(1) maximum temperature (K)
         * @param coeffs(2) standard-state pressure (Pa)
         * @param coeffs(3)-coeffs(11) polynomial coefficients
         */

        TigerPolynomial(int n, const vector_fp& coeffs) :

            m_lowT        (coeffs[0]),
            m_highT       (coeffs[1]),
            m_Pref        (coeffs[2]),
            m_coeff       (vector_fp(8)) {

            m_coeff[0]    = coeffs[3];
            m_coeff[1]    = 1.e-3 * coeffs[4];
            m_coeff[2]    = 1.e-6 * coeffs[5];
            m_coeff[3]    = 1.e-9 * coeffs[6];
            m_coeff[4]    = 1.e3 * coeffs[7];
            m_coeff[5]    = 1.e6 * coeffs[8];
            m_coeff[6]    = 1.e9 * coeffs[9]; 
            m_coeff[7]    = 1.e9 * coeffs[10] / 
                            (GasConst_cal_mol_K * 298.15);

            if (Debug::on) {
                Debug::require(coeffs.size(), 11, Debug::NotLessThan);
            }
        }

        doublereal minTemp() const { return m_lowT; }
        doublereal maxTemp() const { return m_highT; }
        doublereal refPressure() const { return m_Pref; }

        void updateProperties(const TempFuncMgr& tt, 
            vector_fp& cp_R, vector_fp& h_RT, vector_fp& s_R,
            vector_fp& g_RT) const  {

            doublereal ct0 = m_coeff[0];
            doublereal ct1 = m_coeff[1]*tt.value(0);
            doublereal ct2 = m_coeff[2]*tt.value(T2_INDEX);
            doublereal ct3 = m_coeff[3]*tt.value(T3_INDEX);
            doublereal ctm1 = m_coeff[4]*tt.value(TRECIP_INDEX);
            doublereal ctm2 = m_coeff[5]/tt.value(T2_INDEX);
            doublereal ctm3 = m_coeff[6]/tt.value(T3_INDEX);
                        
            cp_R[m_index] = ct0 + ct1 + ct2 + ct3 
                            + ctm1 + ctm2 + ctm3;

            h_RT[m_index] = ct0 + 0.5 * ct1 +  ct2/3.0 
                            + 0.25 * ct3 
                            + m_coeff[4] * tt.value(TLOG_INDEX) / tt.value(0)
                            - ctm2 - 0.5 * ctm3
                            + m_coeff[7]*tt.value(TRECIP_INDEX);

            s_R[m_index] = m_coeff[0] * tt.value(TLOG_INDEX) + ct1 
                           + 0.5 * ct2 + ct3/3.0
                           - ctm1 - 0.5 * ctm2 - ctm3/3.0
                           + m_coeff[8] / GasConst_cal_mol_K ;

            g_RT[m_index] = h_RT[m_index] - s_R[m_index];
        }

    protected:

        doublereal m_lowT, m_highT, m_Pref;
        vector_fp m_coeff;
        int m_index;

    private:

    };

}

#endif






