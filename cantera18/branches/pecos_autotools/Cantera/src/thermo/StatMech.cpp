/**
 *  @file StatMech.cpp
 *  \link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType\endlink 
 */

/* $Author: hkmoffa $
 * $Revision: 279 $
 * $Date: 2009-12-05 13:08:43 -0600 (Sat, 05 Dec 2009) $
 */
// Copyright 2007  Sandia National Laboratories

#include "StatMech.h"

namespace Cantera {

 
  // The NASA 9 polynomial parameterization for one temperature range.
  /*
   * @ingroup spthermo
   */
 

  //! Empty constructor
  StatMech::StatMech() 
    : m_lowT(0.0), m_highT (0.0),
      m_Pref(1.0E5), m_index (0), m_coeff(array_fp(9)) {}


  // constructor used in templated instantiations
  /*
   * @param n            Species index
   * @param tlow         Minimum temperature
   * @param thigh        Maximum temperature
   * @param pref         reference pressure (Pa).
   * @param coeffs       Vector of coefficients used to set the
   *                     parameters for the standard state.
   */
  StatMech::StatMech(int n, doublereal tlow, doublereal thigh, 
			 doublereal pref,
			 const doublereal* coeffs) :
    m_lowT      (tlow),
    m_highT     (thigh),
    m_Pref      (pref),
    m_index     (n),
    m_coeff     (array_fp(9)) {
    std::copy(coeffs, coeffs + 9, m_coeff.begin());
  }

  // copy constructor
  /*
   * @param b object to be copied
   */
  StatMech::StatMech(const StatMech& b) :
    m_lowT      (b.m_lowT),
    m_highT     (b.m_highT),
    m_Pref      (b.m_Pref),
    m_index     (b.m_index),
    m_coeff     (array_fp(9)) {
    std::copy(b.m_coeff.begin(),
	      b.m_coeff.begin() + 9,
	      m_coeff.begin());
  }

  // assignment operator
  /*
   * @param b object to be copied
   */
  StatMech& StatMech::operator=(const StatMech& b) {
    if (&b != this) {
      m_lowT   = b.m_lowT;
      m_highT  = b.m_highT;
      m_Pref   = b.m_Pref;
      m_index  = b.m_index;
      std::copy(b.m_coeff.begin(),
		b.m_coeff.begin() + 9,
		m_coeff.begin());
    }
    return *this;
  }

  // Destructor
  StatMech::~StatMech() {
  }

  // duplicator
  SpeciesThermoInterpType *
  StatMech::duplMyselfAsSpeciesThermoInterpType() const {
    StatMech* np = new StatMech(*this);
    return (SpeciesThermoInterpType *) np;
  }

  // Returns the minimum temperature that the thermo
  // parameterization is valid
  doublereal StatMech::minTemp() const     { return m_lowT;}

  // Returns the maximum temperature that the thermo
  // parameterization is valid
  doublereal StatMech::maxTemp() const  { 
    return m_highT;
  }

  // Returns the reference pressure (Pa)
  doublereal StatMech::refPressure() const { return m_Pref; }

  // Returns an integer representing the type of parameterization
  int StatMech::reportType() const {
    return NASA9;
  }
      
  // Returns an integer representing the species index
  int StatMech::speciesIndex() const { 
    return m_index;
  }

  // Update the properties for this species, given a temperature polynomial
  /*
   * This method is called with a pointer to an array containing the functions of
   * temperature needed by this  parameterization, and three pointers to arrays where the
   * computed property values should be written. This method updates only one value in
   * each array.
   *
   * Temperature Polynomial:
   *  tt[0] = t;
   *  tt[1] = t*t;
   *  tt[2] = t*t*t;
   *  tt[3] = t*t*t*t;
   *  tt[4] = 1.0/t;
   *  tt[5] = 1.0/(t*t);
   *  tt[6] = std::log(t);
   *
   * @param tt      vector of temperature polynomials
   * @param cp_R    Vector of Dimensionless heat capacities.
   *                (length m_kk).
   * @param h_RT    Vector of Dimensionless enthalpies.
   *                (length m_kk).
   * @param s_R     Vector of Dimensionless entropies.
   *                (length m_kk).
   */
  void StatMech::updateProperties(const doublereal* tt, 
				    doublereal* cp_R, doublereal* h_RT,
				    doublereal* s_R) const {
          
    doublereal ct0 = m_coeff[0] * tt[5];   // a0 / (T^2)
    doublereal ct1 = m_coeff[1] * tt[4];   // a1 / T
    doublereal ct2 = m_coeff[2];           // a2
    doublereal ct3 = m_coeff[3] * tt[0];   // a3 * T
    doublereal ct4 = m_coeff[4] * tt[1];   // a4 * T^2
    doublereal ct5 = m_coeff[5] * tt[2];   // a5 * T^3
    doublereal ct6 = m_coeff[6] * tt[3];   // a6 * T^4
 

    doublereal cpdivR = ct0 + ct1 + ct2 + ct3 + ct4 + ct5 + ct6;
    doublereal hdivRT = -ct0 + tt[6]*ct1  + ct2 + 0.5*ct3 + OneThird*ct4  
      + 0.25*ct5  + 0.2*ct6 + m_coeff[7] * tt[4]; 
    doublereal sdivR  = -0.5*ct0  - ct1 + tt[6]*ct2  + ct3  + 0.5*ct4 
      + OneThird*ct5 + 0.25*ct6 + m_coeff[8]; 

    // return the computed properties in the location in the output 
    // arrays for this species
    cp_R[m_index] = cpdivR;
    h_RT[m_index] = hdivRT;
    s_R[m_index] = sdivR;
    //writelog("NASA9poly1: for species "+int2str(m_index)+", h_RT = "+
    //    fp2str(h)+"\n");
  }

 
  // Compute the reference-state property of one species
  /*
   * Given temperature T in K, this method updates the values of
   * the non-dimensional heat capacity at constant pressure,
   * enthalpy, and entropy, at the reference pressure, Pref
   * of one of the species. The species index is used
   * to reference into the cp_R, h_RT, and s_R arrays.
   *
   * Temperature Polynomial:
   *  tt[0] = t;
   *  tt[1] = t*t;
   *  tt[2] = t*t*t;
   *  tt[3] = t*t*t*t;
   *  tt[4] = 1.0/t;
   *  tt[5] = 1.0/(t*t);
   *  tt[6] = std::log(t);
   *
   * @param temp    Temperature (Kelvin)
   * @param cp_R    Vector of Dimensionless heat capacities.
   *                (length m_kk).
   * @param h_RT    Vector of Dimensionless enthalpies.
   *                (length m_kk).
   * @param s_R     Vector of Dimensionless entropies.
   *                (length m_kk).
   */
  void StatMech::updatePropertiesTemp(const doublereal temp, 
					doublereal* cp_R, doublereal* h_RT, 
					doublereal* s_R) const {
    double tPoly[7];
    tPoly[0]  = temp;
    tPoly[1]  = temp * temp;
    tPoly[2]  = tPoly[1] * temp;
    tPoly[3]  = tPoly[2] * temp;
    tPoly[4]  = 1.0 / temp;
    tPoly[5]  = tPoly[4] / temp;
    tPoly[6]  = std::log(temp);
    updateProperties(tPoly, cp_R, h_RT, s_R);
  }

  //This utility function reports back the type of 
  // parameterization and all of the parameters for the 
  // species, index.
  /*
   * All parameters are output variables
   *
   * @param n         Species index
   * @param type      Integer type of the standard type
   * @param tlow      output - Minimum temperature
   * @param thigh     output - Maximum temperature
   * @param pref      output - reference pressure (Pa).
   * @param coeffs    Vector of coefficients used to set the
   *                  parameters for the standard state.
   */
  void StatMech::reportParameters(int &n, int &type,
				    doublereal &tlow, doublereal &thigh,
				    doublereal &pref,
				    doublereal* const coeffs) const {
    n = m_index;
    type = NASA9;
    tlow = m_lowT;
    thigh = m_highT;
    pref = m_Pref;
    coeffs[0] = 1;
    coeffs[1] = m_lowT;
    coeffs[2] = m_highT;
    for (int i = 0; i < 9; i++) {
      coeffs[i+3] = m_coeff[i];
    }

  }

  // Modify parameters for the standard state
  /*
   * @param coeffs   Vector of coefficients used to set the
   *                 parameters for the standard state.
   */
  void StatMech::modifyParameters(doublereal* coeffs) {
    for (int i = 0; i < 9; i++) {
      m_coeff[i] = coeffs[i];
    }            
  }


}

