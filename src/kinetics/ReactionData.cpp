//! @file ReactionData.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ReactionData.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

void ArrheniusData::update(const ThermoPhase& bulk) {
    update(bulk.temperature());
}

void CustomFunc1Data::update(const ThermoPhase& bulk) {
    m_temperature = bulk.temperature();
}

}
