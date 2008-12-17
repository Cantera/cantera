/**
 *  @file phasereport.cpp
 *  Output routines for phases
 */

/*  $Author$
 *  $Date$
 *  $Revision$
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
#include <cstdio>
#include "mix_defs.h"

using namespace std;

namespace Cantera {
  
  //! Format a summary of the mixture state for output.
  /*!
   * @param th  ThermoPhase object to create a report about
   * @param show_thermo Boolean indicating whether the thermo functions
   *                    of the phase should be written out
   *
   * @return  Returns a string containing the report
   */
  std::string report(const ThermoPhase& th, bool show_thermo) {
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
  
  //! Format a composition list for output.
  /*!
   *  Depending upon xyc, this creates a list of
   *  mole fractions, mass fractions, or concentrations.
   * 
   *  @param mix  Mixture to create a list from
   *  @param xyc  Integer 0  mole fractions
   *                      1  mass fractions
   *                      2 concentrations 
   *
   * @return Returns a string containing the list
   */
  std::string formatCompList(const Phase& mix, int xyc) {

    const doublereal Threshold = 1.e-20;

    char p[200];
    string s = "";
    int kk = mix.nSpecies();
    array_fp zz(kk);
    switch (xyc) {
    case 0: mix.getMoleFractions(&zz[0]); break;
    case 1: mix.getMassFractions(&zz[0]); break;
    case 2: mix.getConcentrations(&zz[0]); break;
    default: return "error: xyc must be 0, 1, or 2";
    }

    doublereal z;
    int k;
    for (k = 0; k < kk; k++) {
      z = fabs(zz[k]);
      if (z < Threshold) zz[k] = 0.0;
    }

    for (k = 0; k < kk; k++) {
      sprintf(p, "%18s\t %12.6e\n", mix.speciesName(k).c_str(), 
	      zz[k]);
      s += p;
    }
    return s;
  }

}


