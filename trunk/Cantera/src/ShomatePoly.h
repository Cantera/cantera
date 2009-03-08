/**
 *  @file ShomatePoly.h
 *
 * $Author: dggoodwin $
 * $Revision: 1.5 $
 * $Date: 2006/05/06 15:34:18 $
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_SHOMATEPOLY1_H
#define CT_SHOMATEPOLY1_H

#include "SpeciesThermoInterpType.h"

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
     * \hat h^0(T) = A t + \frac{B t^2}{2} + \frac{C t^3}{3} 
     + \frac{D t^4}{4}  - \frac{E}{t}  + F.
     * \f]
     * \f[
     * s^0(T) = A\ln t + B t + \frac{C t^2}{2}  
     + \frac{D t^3}{3} - \frac{E}{2t^2}  + G.
     * \f]
     */

    class ShomatePoly : public SpeciesThermoInterpType {

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

	ShomatePoly(const ShomatePoly& b) :
	    m_lowT      (b.m_lowT),
            m_highT     (b.m_highT),
            m_Pref      (b.m_Pref),
            m_coeff     (array_fp(7)),
	    m_index     (b.m_index) {
            copy(b.m_coeff.begin(),
		 b.m_coeff.begin() + 7,
		 m_coeff.begin());
	}

	ShomatePoly& operator=(const ShomatePoly& b) {
	    if (&b != this) {
	      m_lowT   = b.m_lowT;
	      m_highT  = b.m_highT;
	      m_Pref   = b.m_Pref;
	      m_index  = b.m_index;
	      copy(b.m_coeff.begin(),
		   b.m_coeff.begin() + 7,
		   m_coeff.begin());
	    }
	    return *this;
	}
        virtual ~ShomatePoly(){}

	virtual SpeciesThermoInterpType *
	duplMyselfAsSpeciesThermoInterpType() const {
	    ShomatePoly* sp = new ShomatePoly(*this);
	    return (SpeciesThermoInterpType *) sp;
	}
      
        doublereal minTemp() const { return m_lowT;}
        doublereal maxTemp() const { return m_highT;}
        doublereal refPressure() const { return m_Pref; }
        virtual int reportType() const { return SHOMATE; }
      
        /**
	 *  This formulation calculates the thermo functions
	 *  given the native formulation of the temperature
	 *  polynomial 
	 *
         *  tt is T/1000. 
	 *  m_t[0] = tt;
	 *  m_t[1] = tt*tt;
	 *  m_t[2] = m_t[1]*tt;
	 *  m_t[3] = 1.0/m_t[1];
	 *  m_t[4] = log(tt);
	 *  m_t[5] = 1.0/GasConstant;
	 *  m_t[6] = 1.0/(GasConstant * T);
         */
        void updateProperties(const doublereal* tt, 
			      doublereal* cp_R, doublereal* h_RT, 
			      doublereal* s_R) const {

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

            /*
             *  Shomate polynomials parameterizes assuming units of
             *  J/(gmol*K) for cp_r and s_R and kJ/(gmol) for h.
             *  However, Cantera assumes default MKS units of 
             *  J/(kmol*K). This requires us to multiply cp and s
             *  by 1.e3 and h by 1.e6, before we then nondimensionlize
             *  the results by dividing by (GasConstant * T),
             *  where GasConstant has units of J/(kmol * K).
             */
            cp_R[m_index] = 1.e3 * cp * tt[5];
            h_RT[m_index] = 1.e6 * h  * tt[6];
            s_R[m_index]  = 1.e3 * s  * tt[5];
        }

	/**
	 *   updatePropertiesTemp():
	 *     This formulation creates its own temperature 
	 *     polynomial. Then, it calls updateProperties();
	 *     -> general, but slow.
	 */
	void updatePropertiesTemp(const doublereal temp, 
				  doublereal* cp_R, doublereal* h_RT, 
				  doublereal* s_R) const {
	    double tPoly[7];
	    doublereal tt = 1.e-3*temp;
	    tPoly[0] = tt;
	    tPoly[1] = tt * tt;
	    tPoly[2] = tPoly[1] * tt;
	    tPoly[3] = 1.0/tPoly[1];
	    tPoly[4] = log(tt);
	    tPoly[5] = 1.0/GasConstant;
	    tPoly[6] = 1.0/(GasConstant * temp);
	    updateProperties(tPoly, cp_R, h_RT, s_R);
	}


	void reportParameters(int &n, int &type,
			      doublereal &tlow, doublereal &thigh,
			      doublereal &pref,
			      doublereal* const coeffs) const {
	    n = m_index;
	    type = SHOMATE;
	    tlow = m_lowT;
	    thigh = m_highT;
	    pref = m_Pref;
	    for (int i = 0; i < 7; i++) {
	      coeffs[i] = m_coeff[i];
	    }
	}
    protected:
        
        doublereal m_lowT, m_highT, m_Pref;
        array_fp m_coeff;
        int m_index;
    
    private:

    };



    class ShomatePoly2 : public SpeciesThermoInterpType {
    public:

	ShomatePoly2() 
            : m_lowT(0.0),
	      m_midT(0.0), 
	      m_highT (0.0),
	      m_Pref(0.0),
	      msp_low(0),
	      msp_high(0),
              m_index(0) {
	    m_coeff.resize(15);
	}

	ShomatePoly2(int n, doublereal tlow, doublereal thigh, doublereal pref,
		    const doublereal* coeffs) :
            m_lowT      (tlow),
	    m_midT(0.0),
            m_highT     (thigh),
	    m_Pref      (pref),
	    msp_low(0),
	    msp_high(0),
            m_index     (n)  {
            m_coeff.resize(15);
            copy(coeffs, coeffs + 15, m_coeff.begin());
	    m_midT = coeffs[0];
	    msp_low  = new ShomatePoly(n, tlow, m_midT, pref, coeffs+1);
	    msp_high = new ShomatePoly(n, m_midT, thigh, pref, coeffs+8);
        }

	ShomatePoly2(const ShomatePoly2& b) :
	    m_lowT      (b.m_lowT),
	    m_midT      (b.m_midT),
            m_highT     (b.m_highT),
            m_Pref      (b.m_Pref),
	    msp_low(0),
	    msp_high(0),
            m_coeff     (array_fp(15)),
	    m_index     (b.m_index) {
            copy(b.m_coeff.begin(),
		 b.m_coeff.begin() + 15,
		 m_coeff.begin());
	    msp_low  = new ShomatePoly(m_index, m_lowT, m_midT, 
				       m_Pref, &m_coeff[1]);
	    msp_high = new ShomatePoly(m_index, m_midT, m_highT, 
				       m_Pref, &m_coeff[8]);
	}

	ShomatePoly2& operator=(const ShomatePoly2& b) {
	    if (&b != this) {
	      m_lowT   = b.m_lowT;
	      m_midT   = b.m_midT;
	      m_highT  = b.m_highT;
	      m_Pref   = b.m_Pref;
	      m_index  = b.m_index;
	      copy(b.m_coeff.begin(),
		   b.m_coeff.begin() + 15,
		   m_coeff.begin());
	      if (msp_low) delete msp_low;
	      if (msp_high) delete msp_high;
	      msp_low  = new ShomatePoly(m_index, m_lowT, m_midT, 
					 m_Pref, &m_coeff[1]);
	      msp_high = new ShomatePoly(m_index, m_midT, m_highT, 
					 m_Pref, &m_coeff[8]);
	    }
	    return *this;
	}

	virtual ~ShomatePoly2(){
	    delete msp_low;
	    delete msp_high;
	}
	
	virtual SpeciesThermoInterpType *
	duplMyselfAsSpeciesThermoInterpType() const {
	    ShomatePoly2* sp = new ShomatePoly2(*this);
	    return (SpeciesThermoInterpType *) sp;
	}

	doublereal minTemp() const { return m_lowT;}
        doublereal maxTemp() const { return m_highT;}
        doublereal refPressure() const { return m_Pref; }
        virtual int reportType() const { return SHOMATE2; }

        /**
	 *  This formulation calculates the thermo functions
	 *  given the native formulation of the temperature
	 *  polynomial 
	 *
         *  tt is T/1000. 
	 *  m_t[0] = tt;
	 *  m_t[1] = tt*tt;
	 *  m_t[2] = m_t[1]*tt;
	 *  m_t[3] = 1.0/m_t[1];
	 *  m_t[4] = log(tt);
	 *  m_t[5] = 1.0/GasConstant;
	 *  m_t[6] = 1.0/(GasConstant * T);
         */
	void updateProperties(const doublereal* tt, 
			      doublereal* cp_R, doublereal* h_RT, 
			      doublereal* s_R) const {
	    double T = 1000 * tt[0];
	    if (T <= m_midT) {
	      msp_low->updateProperties(tt, cp_R, h_RT, s_R);
	    } else {
	      msp_high->updateProperties(tt, cp_R, h_RT, s_R);
	    }
	    
	}

	/**
	 *   updatePropertiesTemp():
	 *     This formulation creates its own temperature 
	 *     polynomial. Then, it calls updateProperties();
	 *     -> general, but slow.
	 */
	void updatePropertiesTemp(const doublereal temp, 
				  doublereal* cp_R,
				  doublereal* h_RT, 
				  doublereal* s_R) const {
	    if (temp <= m_midT) {
	      msp_low->updatePropertiesTemp(temp, cp_R, h_RT, s_R);
	    } else {
	      msp_high->updatePropertiesTemp(temp, cp_R, h_RT, s_R);
	    }
	}

	void reportParameters(int &n, int &type,
			      doublereal &tlow, doublereal &thigh,
			      doublereal &pref,
			      doublereal* const coeffs) const {
	    n = m_index;
	    type = SHOMATE2;
	    tlow = m_lowT;
	    thigh = m_highT;
	    pref = m_Pref;
	    for (int i = 0; i < 15; i++) {
	      coeffs[i] = m_coeff[i];
	    }
	}

    protected:
        
        doublereal m_lowT, m_midT;
	doublereal m_highT;
	doublereal m_Pref;
	ShomatePoly *msp_low;
	ShomatePoly *msp_high;
        array_fp m_coeff;
        int m_index;
    };
}

#endif
