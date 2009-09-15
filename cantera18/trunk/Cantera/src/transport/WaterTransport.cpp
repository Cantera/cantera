#include "ct_defs.h"
#include "WaterPropsIAPWS.h"
#include "TransportBase.h"
#include "DenseMatrix.h"
#include "LiquidTransportParams.h"
#include "VPStandardStateTP.h"

#include "WaterTransport.h"
#include "PDSS_Water.h"
#include "WaterSSTP.h"
#include "WaterProps.h"

#include <iostream>
using namespace std;


namespace Cantera {

  //! default constructor
  WaterTransport::WaterTransport(thermo_t* thermo, int ndim) :
    Transport(thermo, ndim)
  {
    initTP();
  }

  //  Copy Constructor for the %WaterThermo object.
  /* 
   *    @param right  ThermoPhase to be copied
   */
  WaterTransport::WaterTransport(const WaterTransport &right) :
    Transport(right.m_thermo, right.m_nDim)
  {
    *this = right;
  }

  // Assignment operator
  /*
   *
   * @param right    Reference to %WaterTransport object to be copied into the
   *                 current one.
   */
  WaterTransport&  WaterTransport::operator=(const  WaterTransport& right)
  {
    if (&right != this) {
      return *this;
    }
    Transport::operator=(right);
    
    // All pointers in this routine are shallow pointers. Therefore, it's
    // ok just to reinitialize them
    initTP();

    return *this;
  }

  // Duplication routine for objects which inherit from %Transport
  /*
   *  This virtual routine can be used to duplicate %Transport objects
   *  inherited from %Transport even if the application only has
   *  a pointer to %Transport to work with.
   *
   *  These routines are basically wrappers around the derived copy
   *  constructor.
   */
  Transport *  WaterTransport::duplMyselfAsTransport() const {
    WaterTransport* tr = new WaterTransport(*this);
    return dynamic_cast<Transport *> (tr);
  }


  // virtual destructor
  WaterTransport::~WaterTransport() {
  }

  // Routine to do some common initializations at the start of using
  // this routine.
  void WaterTransport::initTP() {
    // The expectation is that we have a VPStandardStateTP derived object
    VPStandardStateTP *vpthermo = dynamic_cast<VPStandardStateTP *>(m_thermo);
    if (!vpthermo) {

      WaterSSTP *wsstp = dynamic_cast<WaterSSTP *>(m_thermo);
      if (!wsstp) {
	throw CanteraError("WaterTransport::initTP()",
			   "Expectation is that ThermoPhase be a VPStandardStateTP");
      } else {

	m_sub = wsstp->getWater();
	AssertTrace(m_sub != 0);
	// Get a pointer to a changeable WaterProps object
	m_waterProps = wsstp->getWaterProps();
	AssertTrace(m_waterProps != 0);
      }
    } else {
      m_waterPDSS = dynamic_cast<PDSS_Water *>(vpthermo->providePDSS(0));
      if (!m_waterPDSS) {
	throw CanteraError("WaterTransport::initTP()",
			   "Expectation is that first species be water with a PDSS_Water object");
      }
      // Get a pointer to a changeable WaterPropsIAPWS object
      m_sub = m_waterPDSS->getWater();
      AssertTrace(m_sub != 0);
      // Get a pointer to a changeable WaterProps object
      m_waterProps = m_waterPDSS->getWaterProps();
      AssertTrace(m_waterProps != 0);
    }
  }


  double WaterTransport::viscosity() {
   double visc = m_waterProps->viscosityWater();
   return visc;
  }

}
