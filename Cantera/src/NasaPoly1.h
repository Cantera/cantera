/**
 *  @file NasaPoly1.h
 */

/* $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_NASAPOLY1_H
#define CT_NASAPOLY1_H

namespace Cantera {

    /**
     * The NASA polynomial parameterization for one temperature range.
     * This parameterization expresses the heat capacity as a
     * fourth-order polynomial. Note that this is the form used in the
     * 1971 NASA equilibrium program and by the Chemkin software
     * package, but differs from the form used in the more recent NASA
     * equilibrium program.
     *
     * Seven coefficients \f$(a_0,\dots,a_6)\f$ are used to represent
     * \f$ c_p^0(T)\f$, \f$ h^0(T)\f$, and \f$ s^0(T) \f$ as 
     * polynomials in \f$ T \f$ :  
     * \f[
     * \frac{c_p(T)}{R} = a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4
     * \f]
     * \f[
     * \frac{h^0(T)}{RT} = a_0 + \frac{a_1}{2} T + \frac{a_2}{3} T^2 
     * + \frac{a_3}{4} T^3 + \frac{a_4}{5} T^4  + \frac{a_5}{T}.
     * \f]
     * \f[
     * \frac{s^0(T)}{R} = a_0\ln T + a_1 T + \frac{a_2}{2} T^2 
     + \frac{a_3}{3} T^3 + \frac{a_4}{4} T^4  + a_6.
     * \f]
     * 
     * This class is designed specifically for use by class NasaThermo.
     * @ingroup spthermo
     */
    class NasaPoly1 {

    public:

        NasaPoly1() 
            : m_lowT(0.0), m_highT (0.0),
              m_Pref(0.0), m_index (0), m_coeff(array_fp(7)) {}

        NasaPoly1(int n, doublereal tlow, doublereal thigh, doublereal pref,
            const doublereal* coeffs) :
            m_lowT      (tlow),
            m_highT     (thigh),
            m_Pref      (pref),
            m_index     (n),
            m_coeff     (array_fp(7)) {
            copy(coeffs, coeffs + 7, m_coeff.begin());
        }

        virtual ~NasaPoly1(){}
      
        doublereal minTemp() const     { return m_lowT;}
        doublereal maxTemp() const     { return m_highT;}
        doublereal refPressure() const { return m_Pref; }
      
        /**
         * Update the properties for this species. This method is called 
         * with a pointer to an array containing the functions of temperature needed by this 
         * parameterization, and three pointers to arrays where the computed property values
         * should be written. This method updates only one value in each array.  
         */
        void updateProperties(const doublereal* tt, 
            doublereal* cp_R, doublereal* h_RT, doublereal* s_R) const {
          
            doublereal ct0 = m_coeff[2];          // a0 
            doublereal ct1 = m_coeff[3]*tt[0];    // a1 * T
            doublereal ct2 = m_coeff[4]*tt[1];    // a2 * T^2
            doublereal ct3 = m_coeff[5]*tt[2];    // a3 * T^3
            doublereal ct4 = m_coeff[6]*tt[3];    // a4 * T^4
 
            doublereal cp, h, s;
            cp = ct0 + ct1 + ct2 + ct3 + ct4;
            h = ct0 + 0.5*ct1 + OneThird*ct2 + 0.25*ct3 + 0.2*ct4
                + m_coeff[0]*tt[4];               // last term is a5/T
            s = ct0*tt[5] + ct1 + 0.5*ct2 + OneThird*ct3
                +0.25*ct4 + m_coeff[1];           // last term is a6

            // return the computed properties in the location in the output 
            // arrays for this species
            cp_R[m_index] = cp;
            h_RT[m_index] = h;
            s_R[m_index] = s;
        }

    protected:
        
        doublereal m_lowT;     // lowest valid temperature
        doublereal m_highT;    // highest valid temperature
        doublereal m_Pref;     // standard-state pressure
        int m_index;           // species index
        array_fp m_coeff;      // array of polynomial coefficients
    
    private:

    };

}

#endif






