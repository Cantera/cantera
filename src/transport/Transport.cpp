/**
 *  @file Transport.cpp
 *  Mixture-averaged transport properties for ideal gas mixtures.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/transport/Transport.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/transport/TransportFactory.h"

namespace Cantera
{
Transport::Transport(ThermoPhase* thermo, size_t ndim) :
    m_thermo(thermo),
    m_nDim(ndim)
{
    if (thermo != nullptr) {
        warn_deprecated("Transport::Transport", "Specifying the ThermoPhase object "
            "in the Transport constructor is deprecated and will be removed after "
            "Cantera 3.0");
    }
    if (ndim != npos) {
        warn_deprecated("Transport::Transport", "The 'ndim' argument to the Transport "
            "constructor is deprecated and will be removed after Cantera 3.0");
    } else {
        m_nDim = 1;
    }
}

bool Transport::ready()
{
    warn_deprecated("Transport::ready", "To be removed after Cantera 3.0");
    return m_ready;
}

void Transport::setNDim(const int ndim)
{
    warn_deprecated("Transport::setNDim", "To be removed after Cantera 3.0");
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
    string name = TransportFactory::factory()->canonicalize(transportModel());
    if (name != "") {
        out["transport"] = name;
    }
    return out;
}

void Transport::setThermo(ThermoPhase& thermo)
{
    warn_deprecated("Transport::setThermo", "To be removed after Cantera 3.0");
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
            string newS0 = thermo.speciesName(i);
            string oldS0 = m_thermo->speciesName(i);
            if (newNum != oldNum) {
                throw CanteraError("Transport::setThermo",
                                   "base object cannot be changed after "
                                   "the transport manager has been constructed because species names are not the same");
            }
        }
        m_thermo = &thermo;
    }
}

void Transport::setRoot(shared_ptr<Solution> root)
{
    warn_deprecated("Transport::setRoot", "To be removed after Cantera 3.0");
    m_root = root;
}

void Transport::finalize()
{
    warn_deprecated("Transport::finalize", "To be removed after Cantera 3.0");
    if (!ready()) {
        m_ready = true;
    } else {
        throw CanteraError("Transport::finalize",
                           "finalize has already been called.");
    }
}
}
