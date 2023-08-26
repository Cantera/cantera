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

void Transport::checkSpeciesIndex(size_t k) const
{
    if (k >= m_nsp) {
        throw IndexError("Transport::checkSpeciesIndex", "species", k, m_nsp-1);
    }
}

void Transport::checkSpeciesArraySize(size_t kk) const
{
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
