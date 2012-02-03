#ifndef CXX_IMPORTPHASE
#define CXX_IMPORTPHASE

#include <string>

#include "kernel/ThermoPhase.h"
#include "kernel/ThermoFactory.h"

/*
 * The Definitions for these functions are all located in 
 * libctxx.a
 */

namespace Cantera_CXX {
    Cantera::ThermoPhase* importPhase(std::string infile, std::string id="");

    // -> this is a duplicate of a src/thermo/phasereport function
    //     We'll leave it here so that these are available externally
    std::string report(const Cantera::ThermoPhase& th, bool show_thermo);

    std::string formatCompList(const Cantera::Phase& mix, int xyc);
}


#endif
