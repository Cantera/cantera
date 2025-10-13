/**
 *  @file Transport.cpp
 *  Mixture-averaged transport properties for ideal gas mixtures.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/transport/Transport.h"
#include "cantera/base/AnyMap.h"
#include "cantera/transport/TransportFactory.h"

namespace Cantera
{

shared_ptr<Transport> Transport::clone(shared_ptr<ThermoPhase> thermo) const
{
    return newTransport(thermo, transportModel());
}

size_t Transport::checkSpeciesIndex(size_t k) const
{
    if (k < m_nsp) {
        return k;
    }
    throw IndexError("Transport::checkSpeciesIndex", "species", k, m_nsp);
}

void Transport::checkSpeciesArraySize(size_t kk) const
{
    warn_deprecated("Transport::checkSpeciesArraySize",
        "To be removed after Cantera 3.2. Only used by legacy CLib.");
    if (m_nsp > kk) {
        throw ArraySizeError("Transport::checkSpeciesArraySize", kk, m_nsp);
    }
}

AnyMap Transport::parameters() const
{
    AnyMap out;
    string name = TransportFactory::factory()->canonicalize(transportModel());
    if (name != "") {
        out["transport"] = name;
    }
    return out;
}

}
