/**
 *  @file cxx/src/importPhase.cpp
 */

#include <string>

#include "ThermoPhase.h"
#include "ThermoFactory.h"
//#include "../include/importPhase.h"

namespace Cantera {}

namespace Cantera_CXX {

    Cantera::ThermoPhase* importPhase(std::string infile, std::string id) {
        Cantera::ThermoPhase* p = Cantera::newPhase(infile, id);
        return p;
    }
}


