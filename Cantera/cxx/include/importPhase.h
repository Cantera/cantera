#ifndef CXX_IMPORTPHASE
#define CXX_IMPORTPHASE

#include <string>

#include "kernel/ThermoPhase.h"
#include "kernel/ThermoFactory.h"

/*
 * The Definitions for these functions are all located in 
 * libctxx.a
 */
namespace Cantera {
   ThermoPhase* importPhase(std::string infile, std::string id="");

   std::string report(const ThermoPhase& th, bool show_thermo);
   std::string formatCompList(const Phase& mix, int xyc);
}


#endif
