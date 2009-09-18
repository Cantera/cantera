/**
 *  @file phasereport.cpp
 *  Output routines for phases
 */

/* 
 *  $Date: 2009/01/14 22:46:33 $
 *  $Revision: 1.8 $
 *
 *  Copyright 2001 California Institute of Technology
 *
 */

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ThermoPhase.h"
#include "PureFluidPhase.h"
#include "mix_defs.h"

#include <cstdio>

using namespace std;

namespace Cantera {
  
  // Format a summary of the mixture state for output.
  /*
   * @param th  ThermoPhase object to create a report about
   * @param show_thermo Boolean indicating whether the thermo functions
   *                    of the phase should be written out
   *
   * @return  Returns a string containing the report
   */
  std::string report(const ThermoPhase& th, const bool show_thermo) {
    return th.report(show_thermo);
  }

  //! Write a phase report to the screen device
  /*!
   *  This routine is a wrapper around the report() function.
   *  It writes to the screen device.
   *
   * @param th          ThermoPhase object to create a report about
   * @param show_thermo Boolean indicating whether the thermo functions
   *                    of the phase should be written out
   */
  void writephase(const ThermoPhase& th, bool show_thermo) {
    std::string s = report(th, show_thermo);
    writelog(s+"\n");
  }
  


}


