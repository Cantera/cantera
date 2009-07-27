/**
 *  @file State.cpp
 * Definitions for the class State, that manages the independent variables of temperature, mass density,
 * and species mass/mole fraction that define the thermodynamic state (see \ref phases and
 * class \link Cantera::State State\endlink).
 */

/*
 *  
 *  $Date: 2009/03/03 19:53:34 $
 *  $Revision: 1.8 $
 *
 *  Copyright 2003-2004 California Institute of Technology
 *  See file License.txt for licensing information
 *
 */

#include "utilities.h"
#include "ctexceptions.h"
#include "stringUtils.h"
#include "State.h"

//#ifdef DARWIN
//#include <Accelerate.h>
//#endif

using namespace std;

namespace Cantera {

  inline void State::stateMFChangeCalc(bool forcerChange) {
    // Right now we assume that the mole fractions have changed every time
    // the function is called
    m_stateNum++;
    if (m_stateNum > 1000000) m_stateNum = -10000000;
  }

  State::State() : 
    m_kk(0),
    m_temp(0.0), 
    m_dens(0.001), 
    m_mmw(0.0),
    m_stateNum(-1)
  {
  }

  State::~State() 
  {
  }

  State::State(const State& right) : 
    m_kk(0),
    m_temp(0.0),
    m_dens(0.001), 
    m_mmw(0.0),
    m_stateNum(-1)
  {
    /*
     * Call the assignment operator.
     */
    *this = operator=(right);
  }

  /*
   * Assignment operator for the State Class
   */
  State& State::operator=(const State& right) {
    /*
     * Check for self assignment.
     */
    if (this == &right) return *this;
    /*
     * We do a straight assignment operator on all of the
     * data. The vectors are copied.
     */
    m_kk             = right.m_kk;
    m_temp           = right.m_temp;
    m_dens           = right.m_dens;
    m_mmw            = right.m_mmw;
    m_ym             = right.m_ym;
    m_y              = right.m_y;
    m_molwts         = right.m_molwts;
    m_rmolwts        = right.m_rmolwts;
    m_stateNum = -1;
    /*
     * Return the reference to the current object
     */
    return *this;
  }

  doublereal State::moleFraction(const int k) const {
    if (k >= 0 && k < m_kk) {
      return m_ym[k] * m_mmw;
    }
    else {
      throw CanteraError("State:moleFraction",
			 "illegal species index number");
    }
  }

  void State::setMoleFractions(const doublereal* const x) {
    doublereal sum = dot(x, x + m_kk, m_molwts.begin());
    doublereal rsum = 1.0/sum;
    transform(x, x + m_kk, m_ym.begin(), timesConstant<double>(rsum));
    transform(m_ym.begin(), m_ym.begin() + m_kk, m_molwts.begin(), 
	      m_y.begin(), multiplies<double>());
    doublereal norm = accumulate(x, x + m_kk, 0.0);
    m_mmw = sum/norm;
   
    //! Call a routine to determin whether state has changed.
    stateMFChangeCalc();
  }

  void State::setMoleFractions_NoNorm(const doublereal* const x) {
    m_mmw = dot(x, x + m_kk, m_molwts.begin());
    doublereal rmmw = 1.0/m_mmw;
    transform(x, x + m_kk, m_ym.begin(), timesConstant<double>(rmmw));
    transform(m_ym.begin(), m_ym.begin() + m_kk, m_molwts.begin(), 
	      m_y.begin(), multiplies<double>());

    //! Call a routine to determin whether state has changed.
    stateMFChangeCalc();
  }

  doublereal State::massFraction(const int k) const {
    if (k >= 0 && k < m_kk) {
      return m_y[k];
    }
    else {
      throw CanteraError("State:massFraction",
			 "illegal species index number");
    }
  }

  doublereal State::concentration(const int k) const {
    if (k >= 0 && k < m_kk) {
      return m_y[k] * m_dens * m_rmolwts[k] ;
    }
    else {
      throw CanteraError("State:massFraction",
			 "illegal species index number");
    }
  }

  void State::setMassFractions(const doublereal* const y) {
    doublereal norm = 0.0, sum = 0.0;
    //cblas_dcopy(m_kk, y, 1, m_y.begin(), 1);
    norm = accumulate(y, y + m_kk, 0.0);
    copy(y, y + m_kk, m_y.begin());
    scale(y, y + m_kk, m_y.begin(), 1.0/norm);
    //        for (k = 0; k != m_kk; ++k) {
    //    norm += y[k];
    //    m_y[k] = y[k];
    //}

    //scale(m_kk, 1.0/norm, m_y.begin());
    transform(m_y.begin(), m_y.begin() + m_kk, m_rmolwts.begin(),
	      m_ym.begin(), multiplies<double>());
    sum = accumulate(m_ym.begin(), m_ym.begin() + m_kk, 0.0);
    //    for (k = 0; k != m_kk; ++k) {
    //      m_ym[k] = m_y[k] * m_rmolwts[k];
    //      sum += m_ym[k];
    //    }
    m_mmw = 1.0/sum;

    //! Call a routine to determin whether state has changed.
    stateMFChangeCalc();
  }

  void State::setMassFractions_NoNorm(const doublereal* const y) {
    doublereal sum = 0.0;
    copy(y, y + m_kk, m_y.begin());
    transform(m_y.begin(), m_y.end(), m_rmolwts.begin(), m_ym.begin(),
	      multiplies<double>());
    sum = accumulate(m_ym.begin(), m_ym.end(), 0.0);
    //for (k = 0; k != m_kk; ++k) {
    //    m_y[k] = y[k];
    //    m_ym[k] = m_y[k] * m_rmolwts[k];
    //    sum += m_ym[k];
    //}
    m_mmw = 1.0/sum;

    //! Call a routine to determine whether state has changed.
    stateMFChangeCalc();
  }

  doublereal State::sum_xlogx() const {
    return m_mmw* Cantera::sum_xlogx(m_ym.begin(), m_ym.end()) + log(m_mmw);
  }

  doublereal State::sum_xlogQ(doublereal* Q) const {
    return m_mmw * Cantera::sum_xlogQ(m_ym.begin(), m_ym.end(), Q);
  }

  doublereal State::molarDensity() const { 
    return density()/meanMolecularWeight(); 
  }

  void State::setConcentrations(const doublereal* const c) {
    int k;
    doublereal sum = 0.0, norm = 0.0;
    for (k = 0; k != m_kk; ++k) {
      sum += c[k]*m_molwts[k];
      norm += c[k];
    }
    m_mmw = sum/norm;
    setDensity(sum);
    doublereal rsum = 1.0/sum;
    for (k = 0; k != m_kk; ++k) {
      m_ym[k] = c[k] * rsum;
      m_y[k] =  m_ym[k] * m_molwts[k];
    }

    //! Call a routine to determin whether state has changed.
    stateMFChangeCalc();
  }

  const doublereal* State::moleFractdivMMW() const { 
    return &m_ym[0];
  }

  void State::getConcentrations(doublereal* const c) const {
    scale(m_ym.begin(), m_ym.end(), c, m_dens);
  }

  doublereal State::mean_X(const doublereal* const Q) const {
    return m_mmw*std::inner_product(m_ym.begin(), m_ym.end(), Q, 0.0);
  }

  doublereal State::mean_Y(const doublereal* const Q) const {
    return dot(m_y.begin(), m_y.end(), Q);
  }

  void State::getMoleFractions(doublereal* const x) const {
    scale(m_ym.begin(), m_ym.end(), x, m_mmw);
  }

  void State::getMassFractions(doublereal* const y) const {
    copy(m_y.begin(), m_y.end(), y);
  }

  void State::setMolarDensity(const doublereal molarDensity) {
    m_dens = molarDensity*meanMolecularWeight();
  }


  void State::init(const array_fp& mw) {
    m_kk = mw.size();
    m_molwts.resize(m_kk);
    m_rmolwts.resize(m_kk);
    m_y.resize(m_kk, 0.0);
    m_ym.resize(m_kk, 0.0);
    copy(mw.begin(), mw.end(), m_molwts.begin());
    for (int k = 0; k < m_kk; k++) {
      if (m_molwts[k] < 0.0) {
	throw CanteraError("State::init",
			   "negative molecular weight for species number "+int2str(k));
      }
      /*
       * Some surface phases may define species representing
       * empty sites that have zero molecular weight. Give them
       * a very small molecular weight to avoid dividing by
       * zero.
       */
      if (m_molwts[k] < Tiny) m_molwts[k] = Tiny;
      m_rmolwts[k] = 1.0/m_molwts[k];
    }

    /*
     * Now that we have resized the State object, let's fill it with
     * a valid mass fraction vector that sums to one. The State object
     * should never have a mass fraction vector that doesn't sum to one.
     * We will assume that species 0 has a mass fraction of 1.0 and
     * mass fraction of all other species is 0.0.
     */
    m_y[0] = 1.0;
    m_ym[0] = m_y[0] * m_rmolwts[0];
    m_mmw = 1.0 / m_ym[0];
  }

  // True if the number of species has been set and fixed
  bool State::ready() const { 
    return (m_kk > 0); 
  }


}
