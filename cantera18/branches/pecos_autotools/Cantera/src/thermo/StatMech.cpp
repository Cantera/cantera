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

 
  // Statistical mechanics 
  /*
   * @ingroup spthermo
   */
 

  //! Empty constructor
  StatMech::StatMech() 
    : m_lowT(0.0), m_highT (0.0),
      m_Pref(1.0E5), m_index (0) {}


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
    m_index     (n)
  {

  }

  // copy constructor
  /*
   * @param b object to be copied
   */
  StatMech::StatMech(const StatMech& b) :
    m_lowT      (b.m_lowT),
    m_highT     (b.m_highT),
    m_Pref      (b.m_Pref),
    m_index     (b.m_index)
  {

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

  // Update the properties for this species
  /**
   *
   * \f[
   * \frac{C_p^0(T)}{R} = \frac{C_v^0(T)}{R} + 1 
   * \f]
   *
   * Where,
   * \f[
   * \frac{C_v^0(T)}{R} = \frac{C_v^{tr}(T)}{R} + \frac{C_v^{vib}(T)}{R}
   * \f]
   *
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
    
    // translational + rotational specific heat
    doublereal ctr = 0.0;
    double atomicity = 0.0;

    // 5/2 * R for molecules, 3/2 * R for atoms
    if(atomicity < 1.0) // atom
      {
	ctr += 3/2 * GasConstant; 
      }
    else  // molecule
      {
	ctr += 5/2 * GasConstant; 
	doublereal theta = 1.0;
	ctr += GasConstant * theta* ( theta* exp(theta/tt[0])/(tt[0]*tt[0]))/((exp(theta/tt[0])-1) * (exp(theta/tt[0])-1));
      }   
    doublereal cpdivR = ctr/GasConstant + 1;

    // ACTUNG: fix enthalpy and entropy 
    doublereal hdivRT = 0.0;
    doublereal sdivR  = 0.0;

    // return the computed properties in the location in the output 
    // arrays for this species
    cp_R[m_index] = cpdivR;
    h_RT[m_index] = hdivRT;
    s_R[m_index]  = sdivR;
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
      coeffs[i+3] = 0.0;
    }

  }

  // Modify parameters for the standard state
  /*
   * @param coeffs   Vector of coefficients used to set the
   *                 parameters for the standard state.
   */
  void StatMech::modifyParameters(doublereal* coeffs) 
  {

    

  }


}

