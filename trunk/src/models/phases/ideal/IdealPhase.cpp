/**
 *  @file IdealPhase.cpp
 *
 *  Implementation file for class IdealPhase
 */


#include "IdealPhase.h"

namespace Cantera {
    
    Real IdealPhase::activity(int k) {
        return moleFraction(k);
    }

}
