/**
 *  @file PolyThermo.h
 *
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#include "utilities.h"

#ifndef CT_POLYTHERMO_H
#define CT_POLYTHERMO_H

namespace Cantera {

    /**
     * A polynomial parameterization for one temperature range. 
     * Seven coefficients \f$(a_0,\dots,a_6)\f$ are used to represent
     * \f$ c_p^0(T)\f$, \f$ h^0(T)\f$, and \f$ s^0(T) \f$ as 
     * polynomials in \f$ T \f$ :  
     * \f[
     * \frac{c_p(T)}{R} = a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4
     * \f]
     * \f[
     * \frac{h^0(T)}{RT} = a_0 + \frac{a_1}{2} T + \frac{a_2}{3} T^2 
     + \frac{a_3}{4} T^3 + \frac{a_4}{5} T^4  + a_5.
     * \f]
     * \f[
     * \frac{s^0(T)}{R} = a_0\ln T + a_1 T + \frac{a_2}{2} T^2 
     + \frac{a_3}{3} T^3 + \frac{a_4}{4} T^4  + a_6.
     * \f]
     * 
     * This class is designed specifically for use by class NasaThermo.
     */

    template <int N>
    class PolyThermo {

    public:

        PolyThermo() {}

        void setCoefficients(const vector_fp& coeffs) {
            m_c.resize(N + 1);
            m_h.resize(N + 1);
            m_s.resize(N + 1);
            //copy(coeffs.begin(), coeffs.begin() + N + 3, m_c.begin());
            m_c[0] = coeffs[2];
            m_h[0] = m_c[0];
            m_s[0] = coeffs[1];
            for (int i = 1; i <= N; i++) {
                m_c[i] = coeffs[i+2];
                m_h[i] = m_c[i]/(i+1);
                m_s[i] = m_c[i]/i;
            }
            m_h0 = coeffs[0];
            m_s0 = coeffs[2];
        }

        virtual ~PolyThermo(){}
      
        // these inline functions are defined only for readability
        static doublereal tlog(const doublereal* tt) { return tt[0]; }
        static doublereal trecip(const doublereal* tt) { return tt[1]; }

        doublereal h0_RT(const doublereal* tt) const { 
            return m_c[0] * trecip(tt) + m_c[2]; }
        doublereal s0_R(const doublereal* tt) const { 
            return m_c[2] * tlog(tt) + m_c[1]; }

        void updateProperties(const doublereal* tt, 
            doublereal* cp_R, doublereal* h_RT, doublereal* s_R) const {
            doublereal cpt, ht, st;
            cpt = m_c[0] + dot4(tt+2,m_c.begin()+1);
            ht = m_h0*tt[1] + m_h[0] + dot4(tt+2,m_h.begin()+1);
            st  = m_s0*tt[0] + m_s[0] + dot4(tt+2,m_s.begin()+1); //s0_R(tt);
            *cp_R = cpt;
            *h_RT = ht;
            *s_R  = st;
        }

    protected:

        doublereal m_h0, m_s0;
        vector_fp m_c, m_h, m_s;
    };

}

#endif






