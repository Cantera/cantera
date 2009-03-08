/**
 *  @file ConstCpPoly.h
 *
 * $Author: hkmoffa $
 * $Revision: 1.1 $
 * $Date: 2005/10/22 00:18:16 $
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_CONSTCPPOLY_H
#define CT_CONSTCPPOLY_H

#include "SpeciesThermoInterpType.h"

namespace Cantera {


    class ConstCpPoly: public SpeciesThermoInterpType {

    public:

        ConstCpPoly();
        ConstCpPoly(int n, doublereal tlow, doublereal thigh, 
		    doublereal pref,
		    const doublereal* coeffs);
	ConstCpPoly(const ConstCpPoly&);
	ConstCpPoly& operator=(const ConstCpPoly&);
        virtual ~ConstCpPoly();
	virtual SpeciesThermoInterpType *
    	   duplMyselfAsSpeciesThermoInterpType() const;

        doublereal minTemp() const;
        doublereal maxTemp() const;
        doublereal refPressure() const;
        virtual int reportType() const { return CONSTANT_CP; }
      
        void updateProperties(const doublereal* tt, 
			      doublereal* cp_R, doublereal* h_RT, 
			      doublereal* s_R) const;

        void updatePropertiesTemp(const doublereal temp, 
				  doublereal* cp_R, doublereal* h_RT,
				  doublereal* s_R) const;

	void reportParameters(int &n, int &type,
			      doublereal &tlow, doublereal &thigh,
			      doublereal &pref,
			      doublereal* const coeffs) const;
    protected:
	doublereal m_t0;
        doublereal m_cp0_R;
	doublereal m_h0_R;
	doublereal m_s0_R;
	doublereal m_logt0;
        doublereal m_lowT, m_highT, m_Pref;
        int m_index;
    
    private:

    };

}

#endif
