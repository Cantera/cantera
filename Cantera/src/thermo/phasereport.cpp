
// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ThermoPhase.h"
#include "PureFluidPhase.h"
#include <stdio.h>
#include "mix_defs.h"

using namespace std;

namespace Cantera {

  /**
   * Format a summary of the mixture state for output.
   */           
  std::string report(const ThermoPhase& th, bool show_thermo) {
    return th.report(show_thermo);
  }

  void writephase(const ThermoPhase& th, bool show_thermo) {
    string s = report(th, show_thermo);
    writelog(s+"\n");
  }

  /**
   * Format a composition list for output.
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


