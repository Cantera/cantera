/**
 *  @file NasaPoly1.h
 */

/* $Author: dggoodwin $
 * $Revision: 1.2 $
 * $Date: 2006/04/28 17:22:23 $
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_NASAPOLY2_H
#define CT_NASAPOLY2_H

#include "SpeciesThermoInterpType.h"

namespace Cantera {

    /**
     *
     * 
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
     * This class is designed specifically for use by class 
     * GeneralSpeciesThermo.
     * @ingroup spthermo
     */
    class NasaPoly2 : public SpeciesThermoInterpType {

    public:

        NasaPoly2() 
            : m_lowT(0.0),
	      m_midT(0.0), 
	      m_highT (0.0),
              m_Pref(0.0),
	      mnp_low(0),
	      mnp_high(0), 
	      m_index(0),
	      m_coeff(array_fp(15)) {
	}

        NasaPoly2(int n, doublereal tlow, doublereal thigh, doublereal pref,
            const doublereal* coeffs) :
            m_lowT(tlow),
            m_highT(thigh),
            m_Pref(pref),
	    mnp_low(0),
	    mnp_high(0),
            m_index(n),
            m_coeff(array_fp(15)) {

            copy(coeffs, coeffs + 15, m_coeff.begin());
	    m_midT = coeffs[0];
	    mnp_low  = new NasaPoly1(m_index, m_lowT, m_midT,
				     m_Pref, &m_coeff[1]);
	    mnp_high = new NasaPoly1(m_index, m_midT, m_highT,
				     m_Pref, &m_coeff[8]);
        }

	NasaPoly2(const NasaPoly2& b) :
	    m_lowT(b.m_lowT),
	    m_midT(b.m_midT),
            m_highT(b.m_highT),
            m_Pref(b.m_Pref),
	    mnp_low(0),
	    mnp_high(0),
            m_index(b.m_index),
            m_coeff(array_fp(15)) {

            copy(b.m_coeff.begin(),
		 b.m_coeff.begin() + 15,
		 m_coeff.begin());
	    mnp_low  = new NasaPoly1(m_index, m_lowT, m_midT,
				     m_Pref, &m_coeff[1]);
	    mnp_high = new NasaPoly1(m_index, m_midT, m_highT,
				     m_Pref, &m_coeff[8]);
	}

	NasaPoly2& operator=(const NasaPoly2& b) {
	    if (&b != this) {
	      m_lowT   = b.m_lowT;
	      m_midT   = b.m_midT;
	      m_highT  = b.m_highT;
	      m_Pref   = b.m_Pref;
	      m_index  = b.m_index;
	      copy(b.m_coeff.begin(),
		   b.m_coeff.begin() + 15,
		   m_coeff.begin());
	      if (mnp_low) delete mnp_low;
	      if (mnp_high) delete mnp_high;
	      mnp_low  = new NasaPoly1(m_index, m_lowT, m_midT, 
				       m_Pref, &m_coeff[1]);
	      mnp_high = new NasaPoly1(m_index, m_midT, m_highT, 
				       m_Pref, &m_coeff[8]);
	    }
	    return *this;
	}

        virtual ~NasaPoly2(){
	    delete mnp_low;
	    delete mnp_high;
	}

	virtual SpeciesThermoInterpType *
	duplMyselfAsSpeciesThermoInterpType() const {
	    NasaPoly2* np = new NasaPoly2(*this);
	    return (SpeciesThermoInterpType *) np;
	}
      
        doublereal minTemp() const     { return m_lowT;}
        doublereal maxTemp() const     { return m_highT;}
        doublereal refPressure() const { return m_Pref; }
        virtual int reportType() const { return NASA2; }
      
        /**
         * Update the properties for this species. This method is called 
         * with a pointer to an array containing the functions of
	 * temperature needed by this 
         * parameterization, and three pointers to arrays where the
	 * computed property values
         * should be written. This method updates only one value in 
	 * each array.  
	 *
	 * Temperature Polynomial:
	 *  tt[0] = t;
	 *  tt[1] = t*t;
	 *  tt[2] = m_t[1]*t;
	 *  tt[3] = m_t[2]*t;
	 *  tt[4] = 1.0/t;
	 *  tt[5] = log(t);
         */
        void updateProperties(const doublereal* tt, 
            doublereal* cp_R, doublereal* h_RT, doublereal* s_R) const {
          
       	    double T = tt[0];
	    if (T <= m_midT) {
	      mnp_low->updateProperties(tt, cp_R, h_RT, s_R);
	    } else {
	      mnp_high->updateProperties(tt, cp_R, h_RT, s_R);
	    }
        }

	/**
	 *   updatePropertiesTemp():
	 *     This formulation creates its own temperature 
	 *     polynomial. Then, it calls updateProperties();
	 *
	 *     (note: this is slow, but it is general)
	 */
	void updatePropertiesTemp(const doublereal temp, 
				  doublereal* cp_R,
				  doublereal* h_RT, 
				  doublereal* s_R) const {
	    if (temp <= m_midT) {
	      mnp_low->updatePropertiesTemp(temp, cp_R, h_RT, s_R);
	    } else {
	      mnp_high->updatePropertiesTemp(temp, cp_R, h_RT, s_R);
	    }
	}

	void reportParameters(int &n, int &type,
			      doublereal &tlow, doublereal &thigh,
			      doublereal &pref,
			      doublereal* const coeffs) const {
	    n = m_index;
	    type = NASA2;
	    tlow = m_lowT;
	    thigh = m_highT;
	    pref = m_Pref;
	    for (int i = 0; i < 15; i++) {
	      coeffs[i] = m_coeff[i];
	    }
	}

    protected:
        
        doublereal m_lowT;     // lowest valid temperature
	doublereal m_midT;
        doublereal m_highT;    // highest valid temperature
        doublereal m_Pref;     // standard-state pressure
	NasaPoly1 *mnp_low;
	NasaPoly1 *mnp_high;
        int m_index;           // species index
        array_fp m_coeff;      // array of polynomial coefficients
    
    private:

    };

}
#endif




