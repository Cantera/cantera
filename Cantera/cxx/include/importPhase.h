#ifndef CXX_IMPORTPHASE
#define CXX_IMPORTPHASE

#include <string>

#include "kernel/ThermoPhase.h"
#include "kernel/ThermoFactory.h"

namespace Cantera {

    ThermoPhase* importPhase(std::string infile, std::string id="") {
        ThermoPhase* p = newPhase(infile, id);
        return p;
    }
}


#endif
