/**
 *  @file ShomatePoly.h
 *
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_SHOMATEPOLY1_H
#define CT_SHOMATEPOLY1_H

namespace Cantera {

    /**
     * The Shomate polynomial parameterization for one temperature range. 
     * Seven coefficients \f$(a_0,\dots,a_6)\f$ are used to represent
     * \f$ c_p^0(T)\f$, \f$ h^0(T)\f$, and \f$ s^0(T) \f$ as 
     * polynomials in \f$ T \f$ :  
     * \f[
     * \hat c_p(T) = A + B t + C t^2 + D t^3 + \frac{E}{t^2}
     * \f]
     * \f[
     * \hat h^0(T)} = A t + \frac{B t^2}{2} + \frac{C t^3}{3} 
     + \frac{D t^4}{4}  - \frac{E}{t}  + F.
     * \f]
     * \f[
     * s^0(T) = A\ln t + B t + \frac{C t^2}{2}  
     + \frac{D t^3}{3} - \frac{E}{2t^2}  + G.
     * \f]
     */

    class ShomatePoly {

    public:

        ShomatePoly() 
            : m_lowT(0.0), m_highT (0.0),
              m_Pref(0.0), m_index (0) {m_coeff.resize(7);}

        ShomatePoly(int n, doublereal tlow, doublereal thigh, doublereal pref,
            const doublereal* coeffs) :
            m_lowT      (tlow),
            m_highT     (thigh),
            m_Pref      (pref),
            m_index     (n) {
            m_coeff.resize(7);
            copy(coeffs, coeffs + 7, m_coeff.begin());
        }

        virtual ~ShomatePoly(){}
      
        doublereal minTemp() const { return m_lowT;}
        doublereal maxTemp() const { return m_highT;}
        doublereal refPressure() const { return m_Pref; }
      
        /**
         * t is T/1000.
         *
         * tt[0]      t
         * tt[1]      t*t
         * tt[2]      t*t*t
         * tt[3]      t^4
         * tt[4]      ln t
         */
        void updateProperties(const doublereal* tt, 
            doublereal* cp_R, doublereal* h_RT, doublereal* s_R) const {

            doublereal A      = m_coeff[0]; 
            doublereal Bt     = m_coeff[1]*tt[0];    
            doublereal Ct2    = m_coeff[2]*tt[1];    
            doublereal Dt3    = m_coeff[3]*tt[2]; 
            doublereal Etm2   = m_coeff[4]*tt[3]; 
            doublereal F      = m_coeff[5];
            doublereal G      = m_coeff[6];

            doublereal cp, h, s;
            cp = A + Bt + Ct2 + Dt3 + Etm2;
            h = tt[0]*(A + 0.5*Bt + OneThird*Ct2 + 0.25*Dt3 - Etm2) + F;
            s = A*tt[4] + Bt + 0.5*Ct2 + OneThird*Dt3 - 0.5*Etm2 + G; 
            h *= 1.e6;

            cp_R[m_index] = 1.e3 * cp * tt[5];
            h_RT[m_index] = h * tt[6];
            s_R[m_index] = 1.e3 * s * tt[5];
        }

    protected:
        
        doublereal m_lowT, m_highT, m_Pref;
        array_fp m_coeff;
        int m_index;
    
    private:

    };

}

#endif






