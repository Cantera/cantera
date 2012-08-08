/**
 *  @file TransportBase.cpp
 *  Mixture-averaged transport properties for ideal gas mixtures.
 */
/* 
 * $Revision$
 * $Date$
 */

#include "ThermoPhase.h"
#include "LiquidTransport.h"
#include "ctexceptions.h"

#include "utilities.h"
#include "LiquidTransportParams.h"
#include "TransportFactory.h"
#include "stringUtils.h"

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
    m_nDim(ndim),
    m_velocityBasis(VB_MASSAVG)
  {
  }

  Transport::Transport(const Transport &right) 
  {
    m_thermo        = right.m_thermo;
    m_ready         = right.m_ready;
    m_nmin          = right.m_nmin;
    m_index         = right.m_index;
    m_nDim          = right.m_nDim;
    m_velocityBasis = right.m_velocityBasis;
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
    m_velocityBasis = right.m_velocityBasis;
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

  /* Set an integer index number. This is for internal use of
   * Cantera, and may be removed in the future.
   */
  void Transport::setIndex(int i) {
    m_index = i; 
  }

  // Set the number of dimensions to be expected in flux expressions
  /* Internal memory will be set with this value
   */
  void Transport::setNDim(const int ndim) {
    m_nDim = ndim;
  }



  /* Set transport model parameters. This method may be
   * overloaded in subclasses to set model-specific parameters.
   */
  void Transport::setParameters(const int type, const int k, 
				const doublereal* const p) 
  {
    err("setParameters"); 
  }


  void Transport::setThermo(thermo_t& thermo) { 
    if (!ready()) { 
      m_thermo = &thermo;
      m_nmin = m_thermo->nSpecies();
    }
    else  {
      int newNum = thermo.nSpecies();
      int oldNum = m_thermo->nSpecies();
      if (newNum != oldNum) { 
        throw CanteraError("Transport::setThermo",
                           "base object cannot be changed after "
			   "the transport manager has been constructed because num species isn't the same.");
      }
      for (int i = 0; i < newNum; i++) {
        std::string newS0 = thermo.speciesName(i);
        std::string oldS0 = m_thermo->speciesName(i);
        if (newNum != oldNum) {
          throw CanteraError("Transport::setThermo",
                           "base object cannot be changed after "
                           "the transport manager has been constructed because species names are not the same");
        }
      }
      m_thermo = &thermo;
    }
  }


  doublereal Transport::err(std::string msg) const {

    throw CanteraError("Transport Base Class",
		       "\n\n\n**** Method "+ msg +" not implemented in model "
		       + int2str(model()) + " ****\n"
		       "(Did you forget to specify a transport model?)\n\n\n");
   	 
    return 0.0;
  }

  
  void Transport::finalize() {
    if (!ready()) 
      m_ready = true;
    else 
      throw CanteraError("Transport::finalize",
			 "finalize has already been called.");
  }

  //====================================================================================================================
  void Transport::getSpeciesFluxes(int ndim, const doublereal * const grad_T, 
				   int ldx, const doublereal * const grad_X,
				   int ldf, doublereal * const fluxes) { 
    err("getSpeciesFluxes"); 
  }
  //====================================================================================================================
}
