/**
 *  @file Mu0Poly.h
 */

/* $Author: dggoodwin $
 * $Revision: 1.2 $
 * $Date: 2006/05/06 15:34:18 $
 */



#ifndef CT_MU0POLY_H
#define CT_MU0POLY_H

#include "SpeciesThermoInterpType.h"

namespace Cantera {
    class SpeciesThermo;
    class XML_Node;
    /**
     *   The Mu0Poly class implements a linear interpolation
     *   of the standard state chemical potential of one
     *   species at a single reference pressure.
     *   The chemical potential is input as a series of (T, mu0)
     *   values. The first temperature is assumed to be equal
     *   to 298.15 K; however, this may be relaxed in the future. 
     *   This information, and an assumption of a constant
     *   heat capacity within each interval is enough to 
     *   calculate all thermodynamic functions.
     *   
     *   The basic equation for going from point 1 to point 2
     *  are as follows for T, T1 <= T <= T2
     *
     *      mu1 = H1 - T1 * S1
     * 
     *      mu2 - mu1 = Cp1(T2 - T1) - Cp1(ln(T2/T1)) - S1(T2 - T1)
     *
     *      S2 = S1 + Cp1(ln(T2/T1))
     * 
     *      H2 = H1 + Cp1(T2 - T1)
     *
     *   In the future, a better assumption about the heat
     *   capacity may be employed, so that it can be continuous.
     *
     *   Notes about temperature interpolation for T < T1 and T > Tn
     *     These are achieved by assuming a constant heat capacity
     *     equal to the value in the closest temperature interval.
     *     No error is thrown.
     */
    class Mu0Poly: public SpeciesThermoInterpType {

    public:

        Mu0Poly();

        Mu0Poly(int n, doublereal tlow, doublereal thigh, 
		doublereal pref, const doublereal* coeffs);
	Mu0Poly(const Mu0Poly &);
	Mu0Poly& operator=(const Mu0Poly&);
        virtual ~Mu0Poly();
	SpeciesThermoInterpType *
	duplMyselfAsSpeciesThermoInterpType() const;
	
        doublereal minTemp() const;
        doublereal maxTemp() const;
        doublereal refPressure() const;
        virtual int reportType() const { return MU0_INTERP; }

	/**
	 * Update all of the properties, using the polynomial
	 * tPoly[]
	 *
	 * tPoly[0] = temp (Kelvin)
	 */
        void updateProperties(const doublereal* tPoly, 
			      doublereal* cp_R, doublereal* h_RT, 
			      doublereal* s_R) const ;

	void updatePropertiesTemp(const doublereal temp, 
				  doublereal* cp_R,
				  doublereal* h_RT, 
				  doublereal* s_R) const ;

	/**
	 * report all of the parameters that make up this 
	 * interpolation.
	 */
	void reportParameters(int &n, int &type,
			      doublereal &tlow, doublereal &thigh,
			      doublereal &pref,
			      doublereal* const coeffs) const;
    protected:
	/**
	 * Number of intervals in the interpolating linear
	 * approximation. Number of points is one more than the
	 * number of intervals.
	 */
	int m_numIntervals;
	/**
	 * Value of the enthalpy at T = 298.15.
	 *  This value is tied to the Heat of formation of 
	 *  the species at 298.15.
	 */
	doublereal m_H298;
	/**
	 * Points at which the standard state chemical potential
	 * are given. 
	 */
	vector_fp m_t0_int;

	/*
	 * Mu0's are primary input data. They aren't strictly
	 * needed, but are kept here for convenience. 
	 */
	vector_fp m_mu0_R_int;
	vector_fp m_h0_R_int;
	vector_fp m_s0_R_int;
	vector_fp m_cp0_R_int;
        doublereal m_lowT, m_highT, m_Pref;
        int m_index;
    
    private:

    };

    void installMu0ThermoFromXML(string speciesName,
				 SpeciesThermo& sp, int k, 
				 const XML_Node* Mu0Node_ptr);
}

#endif






