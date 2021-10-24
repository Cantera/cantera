/**
 *  @file TransportBase.cpp
 *  Mixture-averaged transport properties for ideal gas mixtures.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/transport/TransportBase.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/transport/TransportFactory.h"

using namespace std;

namespace Cantera
{
Transport::Transport(ThermoPhase* thermo, size_t ndim) :
    m_thermo(thermo),
    m_ready(false),
    m_nsp(0),
    m_nDim(ndim),
    m_velocityBasis(VB_MASSAVG)
{
}

bool Transport::ready()
{
    return m_ready;
}

void Transport::setNDim(const int ndim)
{
    m_nDim = ndim;
}

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
    string name = TransportFactory::factory()->canonicalize(transportType());
    if (name != "") {
        out["transport"] = name;
    }
    return out;
}

void Transport::setThermo(ThermoPhase& thermo)
{
    if (!ready()) {
        m_thermo = &thermo;
        m_nsp = m_thermo->nSpecies();
    } else  {
        size_t newNum = thermo.nSpecies();
        size_t oldNum = m_thermo->nSpecies();
        if (newNum != oldNum) {
            throw CanteraError("Transport::setThermo",
                               "base object cannot be changed after "
                               "the transport manager has been constructed because num species isn't the same.");
        }
        for (size_t i = 0; i < newNum; i++) {
            std::string newS0 = thermo.speciesName(i);
            std::string oldS0 = m_thermo->speciesName(i);
            if (newNum != oldNum) {
                throw CanteraError("Transport::setThermo",
                                   "base object cannot be changed after "
                                   "the transport manager has been constructed because species names are not the same");
            }
        }
        m_thermo = &thermo;
    }
}

void Transport::setRoot(std::shared_ptr<Solution> root)
{
    m_root = root;
}

void Transport::finalize()
{
    if (!ready()) {
        m_ready = true;
    } else {
        throw CanteraError("Transport::finalize",
                           "finalize has already been called.");
    }
}
}
