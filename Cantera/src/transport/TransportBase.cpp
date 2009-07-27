/**
 *  @file TransportBase.cpp
 *  Mixture-averaged transport properties for ideal gas mixtures.
 */
/* 
 * $Revision: 1.1 $
 * $Date: 2009/02/15 19:41:33 $
 */

#include "ThermoPhase.h"
#include "LiquidTransport.h"

#include "utilities.h"
#include "LiquidTransportParams.h"
#include "TransportFactory.h"

#include "ctlapack.h"

#include <iostream>
using namespace std;

/** 
 * Mole fractions below MIN_X will be set to MIN_X when computing
 * transport properties.
 */
#define MIN_X 1.e-20


namespace Cantera {

  //////////////////// class LiquidTransport methods //////////////




  Transport::Transport(thermo_t* thermo, int ndim) :
    m_thermo(thermo),
    m_ready(false),
    m_nmin(0),
    m_index(-1),
    m_nDim(ndim)
  {
  }

  Transport::Transport(const Transport &right) 
  {
    m_thermo        = right.m_thermo;
    m_ready         = right.m_ready;
    m_nmin          = right.m_nmin;
    m_index         = right.m_index;
    m_nDim          = right.m_nDim;
  }


  Transport& Transport::operator=(const Transport& right) {
    if (&right != this) {
      return *this; 
    }
    m_thermo        = right.m_thermo;
    m_ready         = right.m_ready;
    m_nmin          = right.m_nmin;
    m_index         = right.m_index;
    m_nDim          = right.m_nDim;
    return *this;
  }

  Transport *Transport::duplMyselfAsTransport() const {
    Transport* tr = new Transport(*this);
    return tr;
  }


  Transport::~Transport() { 
  }          
    
  bool Transport::ready() {
    return m_ready; 
  }
  
  int Transport::index() const {
    return m_index; 
  }

  /*
   * Set an integer index number. This is for internal use of
   * Cantera, and may be removed in the future.
   */
  void Transport::setIndex(int i) {
    m_index = i; 
  }

  //! Set the number of dimensions to be expected in flux expressions
  /*!
   * Internal memory will be set with this value
   */
  void Transport::setNDim(const int ndim) {
    m_nDim = ndim;
  }




  /*
   * Set transport model parameters. This method may be
   * overloaded in subclasses to set model-specific parameters.
   */
  void Transport::setParameters(const int type, const int k, 
				const doublereal* const p) 
  {
    err("setParameters"); 
  }
  
}
