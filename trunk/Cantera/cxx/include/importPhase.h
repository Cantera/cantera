#ifndef CXX_IMPORTPHASE
#define CXX_IMPORTPHASE

#include <string>

#include "kernel/ThermoPhase.h"
#include "kernel/importCTML.h"

namespace Cantera {

    ThermoPhase* importPhase(string infile, string id="") {
        ThermoPhase* p = newPhase(infile, id);
        return p;
    }
}


#endif
