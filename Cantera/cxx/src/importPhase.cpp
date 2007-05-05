/**
 *  @file cxx/src/importPhase.cpp
 */

#include <string>

#include "ThermoPhase.h"
#include "ThermoFactory.h"
//#include "../include/importPhase.h"

namespace Cantera {

    ThermoPhase* importPhase(std::string infile, std::string id) {
        ThermoPhase* p = newPhase(infile, id);
        return p;
    }
}


